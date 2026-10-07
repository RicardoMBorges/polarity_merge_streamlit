import io
import re
import math
import json
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
import streamlit as st
import plotly.graph_objects as go

H_MASS = 1.007276466621

st.set_page_config(page_title="Polarity Merge", page_icon="↔", layout="wide")

st.markdown("""
<style>
.block-container {padding-top: 1.4rem; padding-bottom: 2rem;}
.small-note {font-size: 0.86rem; color: #666;}
div[data-testid="stMetric"] {border: 1px solid rgba(128,128,128,.22); padding: 10px 14px; border-radius: 10px;}
</style>
""", unsafe_allow_html=True)


# ----------------------------- Parsers -----------------------------

def read_graphml(uploaded):
    raw = uploaded.getvalue()
    try:
        g = nx.read_graphml(io.BytesIO(raw))
    except Exception:
        # fallback for GraphML readers that prefer a path-like text stream
        g = nx.parse_graphml(raw.decode("utf-8", errors="replace"))
    return g


def node_table(g, polarity):
    rows = []
    for node_id, d in g.nodes(data=True):
        def first_num(*keys):
            for k in keys:
                if k in d:
                    try:
                        return float(d[k])
                    except Exception:
                        pass
            return np.nan
        rows.append({
            "node_id": str(node_id),
            "polarity": polarity,
            "mz": first_num("mz", "precursor_mz", "parent mass", "parent_mass"),
            "rt_min": first_num("rt_min", "rt", "retention_time", "RT"),
            **{f"attr__{k}": v for k, v in d.items()}
        })
    return pd.DataFrame(rows)


def parse_mgf(uploaded):
    text = uploaded.getvalue().decode("utf-8", errors="replace")
    spectra, block = [], None
    for line in text.splitlines():
        s = line.strip()
        if s == "BEGIN IONS":
            block = {"meta": {}, "mz": [], "intensity": []}
        elif s == "END IONS":
            if block is not None:
                spectra.append(block)
            block = None
        elif block is not None:
            if "=" in s:
                k, v = s.split("=", 1)
                block["meta"][k.strip()] = v.strip()
            elif s and not s.startswith("#"):
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        block["mz"].append(float(parts[0]))
                        block["intensity"].append(float(parts[1]))
                    except ValueError:
                        pass
    return spectra


def spectrum_index(spectra):
    idx = {}
    for sp in spectra:
        m = sp["meta"]
        keys = [m.get("FEATURE_ID"), m.get("SCANS"), m.get("FEATURELIST_FEATURE_ID")]
        for k in keys:
            if k is not None:
                idx[str(k)] = sp
        # Also extract a terminal integer from feature-list labels when possible
        fl = m.get("FEATURELIST_FEATURE_ID", "")
        mt = re.search(r"(\d+)\s*$", fl)
        if mt:
            idx[mt.group(1)] = sp
    return idx


# ----------------------------- Matching -----------------------------

def make_candidates(pos, neg, ppm_tol, rt_tol, use_abs_da=False, da_tol=0.005):
    pos = pos.dropna(subset=["mz", "rt_min"]).copy()
    neg = neg.dropna(subset=["mz", "rt_min"]).copy()
    rows = []

    # Neutral mass: [M+H]+ -> M = mz-H ; [M-H]- -> M = mz+H
    neg_nm = neg["mz"].to_numpy(float) + H_MASS
    neg_rt = neg["rt_min"].to_numpy(float)
    neg_ids = neg["node_id"].astype(str).to_numpy()

    for _, p in pos.iterrows():
        pnm = float(p.mz) - H_MASS
        drt = np.abs(neg_rt - float(p.rt_min))
        dm = neg_nm - pnm
        ppm = np.abs(dm) / max(abs(pnm), 1e-12) * 1e6
        mask_mass = np.abs(dm) <= da_tol if use_abs_da else ppm <= ppm_tol
        inds = np.where(mask_mass & (drt <= rt_tol))[0]
        for j in inds:
            rows.append({
                "source": f"POS::{p.node_id}",
                "target": f"NEG::{neg_ids[j]}",
                "pos_node": str(p.node_id),
                "neg_node": str(neg_ids[j]),
                "pos_mz": float(p.mz),
                "neg_mz": float(neg.iloc[j].mz),
                "pos_rt_min": float(p.rt_min),
                "neg_rt_min": float(neg.iloc[j].rt_min),
                "delta_mz_observed": float(p.mz - neg.iloc[j].mz),
                # Signed shift: positive means NEG elutes after POS.
                "delta_rt_min": float(neg.iloc[j].rt_min - p.rt_min),
                "abs_delta_rt_min": float(drt[j]),
                "delta_rt_sec": float((neg.iloc[j].rt_min - p.rt_min) * 60.0),
                "abs_delta_rt_sec": float(drt[j] * 60.0),
                "neutral_mass_pos": pnm,
                "neutral_mass_neg": float(neg_nm[j]),
                "delta_neutral_mass_da": float(dm[j]),
                "mass_error_ppm": float(dm[j] / max(abs(pnm), 1e-12) * 1e6),
                "edge_type": "PolarityMerge",
                "ion_pair": "[M+H]+ / [M-H]-",
            })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["abs_mass_error_ppm"] = out["mass_error_ppm"].abs()
        out = out.sort_values(["abs_mass_error_ppm", "abs_delta_rt_min"]).reset_index(drop=True)
        out["pair_id"] = [f"PM{i+1:05d}" for i in range(len(out))]
        out["accepted"] = True
    return out


# ----------------------------- Optional abundance correlation -----------------------------

def read_table(uploaded):
    name = uploaded.name.lower()
    raw = uploaded.getvalue()
    if name.endswith(".csv"):
        try:
            return pd.read_csv(io.BytesIO(raw))
        except Exception:
            return pd.read_csv(io.BytesIO(raw), sep=None, engine="python")
    if name.endswith((".tsv", ".txt")):
        return pd.read_csv(io.BytesIO(raw), sep=None, engine="python")
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(io.BytesIO(raw))
    raise ValueError("Unsupported abundance-table format.")


def normalize_sample_name(x):
    x = str(x).strip()
    x = re.sub(r"\.(mzxml|mzml|raw|cdf)$", "", x, flags=re.I)
    # conservative POS/NEG token normalization; user can edit mappings
    x = re.sub(r"(?i)(^|[_\-\s])(pos|positive|p|neg|negative|n)(?=([_\-\s]|$))", r"\1POL\3", x)
    return x


def feature_matrix(df, id_col):
    d = df.copy()
    d[id_col] = d[id_col].astype(str)
    numeric = []
    for c in d.columns:
        if c == id_col:
            continue
        s = pd.to_numeric(d[c], errors="coerce")
        if s.notna().sum() > 0:
            d[c] = s
            numeric.append(c)
    return d.set_index(id_col)[numeric]


def correlate_pairs(edges, pos_q, neg_q, pos_id_col, neg_id_col, mapping):
    P = feature_matrix(pos_q, pos_id_col)
    N = feature_matrix(neg_q, neg_id_col)
    valid_map = mapping.dropna(subset=["POS_sample", "NEG_sample"]).copy()
    valid_map = valid_map[
        valid_map["POS_sample"].isin(P.columns) & valid_map["NEG_sample"].isin(N.columns)
    ]
    rvals, nvals = [], []
    for _, e in edges.iterrows():
        pid, nid = str(e.pos_node), str(e.neg_node)
        if pid not in P.index or nid not in N.index or len(valid_map) < 3:
            rvals.append(np.nan); nvals.append(0); continue
        x = np.array([P.loc[pid, c] for c in valid_map.POS_sample], dtype=float)
        y = np.array([N.loc[nid, c] for c in valid_map.NEG_sample], dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 3 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0:
            rvals.append(np.nan); nvals.append(int(ok.sum())); continue
        rvals.append(float(np.corrcoef(x[ok], y[ok])[0, 1]))
        nvals.append(int(ok.sum()))
    z = edges.copy()
    z["intensity_pearson_r"] = rvals
    z["n_paired_samples"] = nvals
    return z



def rt_qc_summary(df):
    if df is None or df.empty:
        return {}
    signed = pd.to_numeric(df["delta_rt_min"], errors="coerce").dropna().to_numpy(float)
    absolute = np.abs(signed)
    if not len(signed):
        return {}
    median = float(np.median(signed))
    mad = float(np.median(np.abs(signed - median)))
    return {
        "median_delta_rt_min": median,
        "median_delta_rt_sec": median * 60.0,
        "mad_delta_rt_min": mad,
        "mad_delta_rt_sec": mad * 60.0,
        "mean_abs_delta_rt_min": float(np.mean(absolute)),
        "mean_abs_delta_rt_sec": float(np.mean(absolute) * 60.0),
        "p95_abs_delta_rt_min": float(np.percentile(absolute, 95)),
        "p95_abs_delta_rt_sec": float(np.percentile(absolute, 95) * 60.0),
    }


def classify_rt_distribution(df, mad_k=3.0):
    out = df.copy()
    vals = pd.to_numeric(out["delta_rt_sec"], errors="coerce")
    finite = vals.dropna().to_numpy(float)
    if len(finite) < 3:
        out["rt_center_sec"] = np.nan
        out["rt_mad_sec"] = np.nan
        out["rt_deviation_sec"] = np.nan
        out["rt_robust_z"] = np.nan
        out["rt_status"] = "Insufficient data"
        return out, {}
    center = float(np.median(finite))
    mad = float(np.median(np.abs(finite-center)))
    dev = np.abs(vals-center)
    if mad > 0:
        lower, upper = center-mad_k*mad, center+mad_k*mad
        rz = dev/mad
        # Inclusive MAD boundary with a tiny numerical tolerance.
        # Prevents values theoretically equal to ±k MAD from being rejected
        # because of binary floating-point representation.
        eps = max(1e-9, np.finfo(float).eps * max(1.0, abs(center), abs(lower), abs(upper)) * 32)
        status = np.where(
            vals.isna(),
            "Missing RT",
            np.where(dev <= (float(mad_k) * mad + eps), "Inlier", "Outlier")
        )
    else:
        lower = upper = center
        rz = pd.Series(np.nan, index=out.index)
        status = np.where(vals.isna(), "Missing RT",
                          np.where(np.isclose(vals, center), "Inlier", "Outlier"))
    out["rt_center_sec"], out["rt_mad_sec"] = center, mad
    out["rt_deviation_sec"], out["rt_robust_z"] = dev, rz
    out["rt_status"] = status
    return out, {"center_sec":center,"mad_sec":mad,"lower_sec":float(lower),"upper_sec":float(upper),
                 "mad_k":float(mad_k),"n_inlier":int(np.sum(status=="Inlier")),
                 "n_outlier":int(np.sum(status=="Outlier"))}


def rt_agreement_plot(df):
    fig = go.Figure()
    if df is None or df.empty:
        return fig
    x = pd.to_numeric(df["pos_rt_min"], errors="coerce")
    y = pd.to_numeric(df["neg_rt_min"], errors="coerce")
    ok = x.notna() & y.notna()
    x, y = x[ok], y[ok]
    if not len(x):
        return fig
    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    fig.add_trace(go.Scatter(
        x=x, y=y, mode="markers",
        customdata=np.column_stack([
            df.loc[ok, "pair_id"].astype(str),
            df.loc[ok, "delta_rt_sec"].astype(float)
        ]),
        hovertemplate="Pair %{customdata[0]}<br>POS RT %{x:.4f} min<br>NEG RT %{y:.4f} min<br>ΔRT %{customdata[1]:+.2f} s<extra></extra>",
        name="Polarity pairs"
    ))
    fig.add_trace(go.Scatter(
        x=[lo, hi], y=[lo, hi], mode="lines",
        name="1:1", hoverinfo="skip"
    ))
    fig.update_layout(
        height=430, xaxis_title="POS RT (min)", yaxis_title="NEG RT (min)",
        title="Retention-time agreement", margin=dict(l=60, r=20, t=55, b=55)
    )
    return fig


def rt_delta_plot(df, info=None):
    fig = go.Figure()
    if df is None or df.empty:
        return fig
    vals = pd.to_numeric(df["delta_rt_sec"], errors="coerce").dropna()
    fig.add_trace(go.Histogram(x=vals, nbinsx=min(40,max(10,int(np.sqrt(max(len(vals),1))*2)))))
    fig.add_vline(x=0, line_dash="dot", annotation_text="0 s")
    if info:
        c, mad = info["center_sec"], info["mad_sec"]
        fig.add_vline(x=c, line_dash="solid", annotation_text=f"Median {c:+.2f} s")
        if mad > 0:
            for k, dash in [(2,"dot"),(3,"dash")]:
                fig.add_vline(x=c-k*mad, line_dash=dash, annotation_text=f"-{k} MAD")
                fig.add_vline(x=c+k*mad, line_dash=dash, annotation_text=f"+{k} MAD")
    fig.update_layout(height=410,title="Signed retention-time shift",
                      xaxis_title="ΔRT = RT NEG − RT POS (s)",yaxis_title="Pair count",
                      showlegend=False,margin=dict(l=60,r=20,t=55,b=55))
    return fig


# ----------------------------- Graph export -----------------------------

def merged_graph(pos_g, neg_g, accepted_edges):
    G = nx.Graph()
    for nid, d in pos_g.nodes(data=True):
        G.add_node(f"POS::{nid}", **{**d, "original_node_id": str(nid), "polarity": "POS"})
    for nid, d in neg_g.nodes(data=True):
        G.add_node(f"NEG::{nid}", **{**d, "original_node_id": str(nid), "polarity": "NEG"})

    for u, v, d in pos_g.edges(data=True):
        G.add_edge(f"POS::{u}", f"POS::{v}", **{**d, "edge_type": d.get("edge_type", "MolecularNetworking"), "polarity": "POS"})
    for u, v, d in neg_g.edges(data=True):
        G.add_edge(f"NEG::{u}", f"NEG::{v}", **{**d, "edge_type": d.get("edge_type", "MolecularNetworking"), "polarity": "NEG"})

    for _, e in accepted_edges.iterrows():
        attrs = {}
        for k, v in e.to_dict().items():
            if k in ("source", "target", "accepted"):
                continue
            if pd.isna(v):
                continue
            if isinstance(v, (np.integer,)): v = int(v)
            if isinstance(v, (np.floating,)): v = float(v)
            if isinstance(v, (np.bool_,)): v = bool(v)
            attrs[k] = v
        G.add_edge(e.source, e.target, **attrs)
    return G


def graphml_bytes(G):
    # networkx GraphML cannot serialize None; all values above are sanitized.
    bio = io.BytesIO()
    nx.write_graphml(G, bio)
    return bio.getvalue()


# ----------------------------- Spectrum plotting -----------------------------

def find_spectrum(idx, node_id):
    if not idx:
        return None
    if str(node_id) in idx:
        return idx[str(node_id)]
    # fallback: FEATURE_ID may not equal GraphML node id; caller can still inspect metadata
    return None


def spectrum_plot(sp, title, top_n=None, min_rel=0.0, label_n=5):
    fig = go.Figure()
    if sp is None or not sp["mz"]:
        fig.add_annotation(text="Spectrum not found for this node ID", x=.5, y=.5, showarrow=False)
        fig.update_xaxes(visible=False); fig.update_yaxes(visible=False)
        fig.update_layout(height=390, title=title)
        return fig

    mz = np.asarray(sp["mz"], float)
    inten = np.asarray(sp["intensity"], float)
    rel = 100 * inten / max(float(np.max(inten)), 1e-12)

    # First apply the minimum relative-intensity threshold.
    keep = rel >= float(min_rel)
    mz, rel = mz[keep], rel[keep]

    # Then retain the N most intense peaks, while plotting them in m/z order.
    if top_n is not None and len(mz) > top_n:
        idx = np.argsort(rel)[-int(top_n):]
        mz, rel = mz[idx], rel[idx]
    order = np.argsort(mz)
    mz, rel = mz[order], rel[order]

    xs, ys = [], []
    for x, y in zip(mz, rel):
        xs += [x, x, None]
        ys += [0, y, None]
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", hoverinfo="skip"))

    # Label only the most intense displayed peaks to avoid clutter.
    if label_n and len(mz):
        nlab = min(int(label_n), len(mz))
        lab_idx = np.argsort(rel)[-nlab:]
        fig.add_trace(go.Scatter(
            x=mz[lab_idx], y=rel[lab_idx],
            mode="text",
            text=[f"{x:.4f}" for x in mz[lab_idx]],
            textposition="top center",
            hovertemplate="m/z %{x:.5f}<br>Relative intensity %{y:.1f}%<extra></extra>"
        ))

    fig.update_layout(
        height=390, title=title,
        xaxis_title="m/z", yaxis_title="Relative intensity (%)",
        yaxis=dict(range=[0, 110]),
        showlegend=False, margin=dict(l=55, r=20, t=55, b=50)
    )
    return fig


# ----------------------------- UI -----------------------------

# LAABio branding
_logo_path = Path(__file__).parent / "assets" / "LAABio_logo.png"
if _logo_path.exists():
    st.sidebar.image(str(_logo_path), use_container_width=True)

with st.sidebar:
    st.markdown("## Polarity Merge")
    st.caption("Cross-polarity molecular-network integration")
    st.markdown("**by Ricardo · IPPN/UFRJ**")
    st.divider()

    st.markdown("### Input networks")
    pos_file = st.file_uploader("Positive GraphML", type=["graphml"], key="pos_graph")
    neg_file = st.file_uploader("Negative GraphML", type=["graphml"], key="neg_graph")

    st.markdown("### MS/MS spectra · optional")
    pos_mgf = st.file_uploader("Positive MGF", type=["mgf"], key="pos_mgf")
    neg_mgf = st.file_uploader("Negative MGF", type=["mgf"], key="neg_mgf")

    st.markdown("### Matching")
    rt_tol = st.number_input("RT tolerance (min)", min_value=0.0, value=0.10, step=0.01, format="%.3f")
    mass_mode = st.selectbox("Mass tolerance", ["ppm", "Da"])
    if mass_mode == "ppm":
        ppm_tol = st.number_input("MS1 tolerance (ppm)", min_value=0.1, value=5.0, step=0.5)
        da_tol = 0.005
    else:
        da_tol = st.number_input("MS1 tolerance (Da)", min_value=0.0001, value=0.0050, step=0.0010, format="%.4f")
        ppm_tol = 5.0

    st.markdown("### Optional abundance evidence")
    use_corr = st.checkbox("Use cross-polarity intensity correlation", value=False)
    pos_quant = neg_quant = None
    if use_corr:
        st.caption("Requires feature × sample abundance tables. The supplied GraphML files do not contain per-sample abundance vectors.")
        pos_quant = st.file_uploader("POS abundance table", type=["csv","tsv","txt","xlsx","xls"])
        neg_quant = st.file_uploader("NEG abundance table", type=["csv","tsv","txt","xlsx","xls"])

st.title("Polarity Merge")
st.write("Integrate positive- and negative-ion molecular networks by adding cross-polarity edges while preserving every original node and molecular-network edge.")

if not pos_file or not neg_file:
    st.info("Upload the positive and negative GraphML networks in the sidebar, choose the parameters, and press **Run Polarity Merge**.")
    st.stop()

st.sidebar.divider()
st.sidebar.subheader("RT distribution filter")
use_rt_distribution = st.sidebar.checkbox("Filter by observed ΔRT distribution", value=True)
mad_k = st.sidebar.select_slider("Acceptance band (± k × MAD)",
    options=[1.5,2.0,2.5,3.0,3.5,4.0,5.0], value=3.0,
    disabled=not use_rt_distribution)
st.sidebar.caption("The initial RT window generates candidates. Final RT support is evaluated around the observed NEG−POS offset. The ±k×MAD boundary is inclusive.")
st.sidebar.divider()
run_merge = st.sidebar.button("Run Polarity Merge", type="primary", use_container_width=True)

# Do not process automatically after upload or parameter changes.
# The last successful run is kept in session_state until the user explicitly runs again.
if run_merge:
    try:
        Gp = read_graphml(pos_file)
        Gn = read_graphml(neg_file)
        P = node_table(Gp, "POS")
        N = node_table(Gn, "NEG")
        edges = make_candidates(P, N, ppm_tol, rt_tol, use_abs_da=(mass_mode=="Da"), da_tol=da_tol)
        if use_rt_distribution:
            edges, rt_filter_info = classify_rt_distribution(edges, mad_k)
        else:
            edges["rt_status"] = "Not filtered"
            rt_filter_info = {}
        st.session_state["merge_run"] = {
            "Gp": Gp,
            "Gn": Gn,
            "P": P,
            "N": N,
            "edges": edges,
            "params": {
                "rt_tol": rt_tol,
                "mass_mode": mass_mode,
                "ppm_tol": ppm_tol,
                "da_tol": da_tol,
                "pos_name": pos_file.name,
                "neg_name": neg_file.name,
                "use_rt_distribution": use_rt_distribution,
                "mad_k": mad_k,
            },
            "rt_filter_info": rt_filter_info,
        }
        # Reset downstream edited results whenever a new merge is executed.
        st.session_state.pop("edges_current", None)
    except Exception as e:
        st.error(f"Could not run Polarity Merge: {e}")
        st.stop()

if "merge_run" not in st.session_state:
    st.info("Inputs are ready. Press **Run Polarity Merge** in the sidebar to calculate cross-polarity matches.")
    st.stop()

run = st.session_state["merge_run"]
Gp, Gn, P, N, edges = run["Gp"], run["Gn"], run["P"], run["N"], run["edges"]
used = run["params"]
rt_filter_info = run.get("rt_filter_info", {})

# Warn if controls were changed after the last explicit run.
changed = (
    used["rt_tol"] != rt_tol or
    used["mass_mode"] != mass_mode or
    used["ppm_tol"] != ppm_tol or
    used["da_tol"] != da_tol or
    used.get("use_rt_distribution", True) != use_rt_distribution or
    used.get("mad_k", 3.0) != mad_k or
    used["pos_name"] != pos_file.name or
    used["neg_name"] != neg_file.name
)
if changed:
    st.warning("Inputs or matching parameters changed after the last run. Press **Run Polarity Merge** to recalculate the results.")

m1, m2, m3, m4 = st.columns(4)
m1.metric("POS nodes", f"{Gp.number_of_nodes():,}")
m2.metric("NEG nodes", f"{Gn.number_of_nodes():,}")
m3.metric("Original edges", f"{Gp.number_of_edges()+Gn.number_of_edges():,}")
m4.metric("Candidate polarity edges", f"{len(edges):,}")

st.caption(
    f'Last run: RT tolerance = {used["rt_tol"]:.3f} min · '
    + (f'MS1 tolerance = {used["ppm_tol"]:.2f} ppm' if used["mass_mode"] == "ppm"
       else f'MS1 tolerance = {used["da_tol"]:.4f} Da')
)

tabs = st.tabs(["Polarity matches", "Sample mapping & correlation", "Spectral Pair Viewer", "Export", "Method"])

with tabs[0]:
    st.subheader("Candidate [M+H]⁺ ↔ [M−H]⁻ links")
    st.caption("Candidates are generated from neutral-mass agreement plus retention-time proximity. Nodes are never collapsed.")
    if edges.empty:
        st.warning("No candidate pairs satisfy the current mass and RT tolerances.")
    else:
        edited = st.data_editor(
            edges,
            use_container_width=True,
            hide_index=True,
            disabled=[c for c in edges.columns if c != "accepted"],
            column_config={"accepted": st.column_config.CheckboxColumn("Keep edge")},
            key="edge_editor"
        )
        st.session_state["edges_current"] = edited
        st.download_button("Download candidate edges (CSV)", edited.to_csv(index=False).encode(),
                           "polarity_merge_edges.csv", "text/csv")

        st.divider()
        st.subheader("Retention-time agreement")
        user_kept = edited[edited["accepted"].astype(bool)].copy()
        if used.get("use_rt_distribution", True):
            rt_review, review_info = classify_rt_distribution(user_kept, used.get("mad_k",3.0))
        else:
            rt_review, review_info = user_kept.copy(), {}
            rt_review["rt_status"] = "Not filtered"

        qc = rt_qc_summary(user_kept)
        if qc:
            q1,q2,q3,q4 = st.columns(4)
            q1.metric("Median ΔRT", f'{qc["median_delta_rt_sec"]:+.2f} s')
            q2.metric("MAD ΔRT", f'{qc["mad_delta_rt_sec"]:.2f} s')
            q3.metric("Mean |ΔRT|", f'{qc["mean_abs_delta_rt_sec"]:.2f} s')
            q4.metric("95th percentile |ΔRT|", f'{qc["p95_abs_delta_rt_sec"]:.2f} s')
            if review_info:
                st.info(f'Observed offset **{review_info["center_sec"]:+.2f} s** · MAD **{review_info["mad_sec"]:.2f} s** · '
                        f'accepted interval **{review_info["lower_sec"]:+.2f} to {review_info["upper_sec"]:+.2f} s** · '
                        f'**{review_info["n_inlier"]} inliers / {review_info["n_outlier"]} outliers**')
            st.caption("ΔRT = RT NEG − RT POS. Filtering is centered on the observed chromatographic offset, not on zero.")
            r1,r2=st.columns(2)
            with r1: st.plotly_chart(rt_agreement_plot(user_kept),use_container_width=True)
            with r2: st.plotly_chart(rt_delta_plot(user_kept,review_info),use_container_width=True)
            cols=[c for c in ["pair_id","pos_node","neg_node","mass_error_ppm","pos_rt_min","neg_rt_min",
                              "delta_rt_sec","rt_deviation_sec","rt_robust_z","rt_status"] if c in rt_review.columns]
            st.dataframe(rt_review[cols],use_container_width=True,hide_index=True)
            st.caption("1:N and N:1 relationships are deliberately preserved; this step filters only chromatographic plausibility.")


with tabs[1]:
    st.subheader("Sample mapping & optional abundance correlation")
    st.write("This layer is optional. It tests whether the POS and NEG features follow the same abundance pattern across matched samples.")
    if not use_corr:
        st.info("Enable **Use cross-polarity intensity correlation** in the sidebar if you want this evidence layer.")
    elif not pos_quant or not neg_quant:
        st.warning("Upload both POS and NEG abundance tables.")
    else:
        try:
            pq, nq = read_table(pos_quant), read_table(neg_quant)
            c1, c2 = st.columns(2)
            with c1:
                pos_id_col = st.selectbox("POS feature-ID column", pq.columns.tolist(), key="pid")
            with c2:
                neg_id_col = st.selectbox("NEG feature-ID column", nq.columns.tolist(), key="nid")

            pcols = [c for c in pq.columns if c != pos_id_col]
            ncols = [c for c in nq.columns if c != neg_id_col]
            auto = []
            n_norm = {normalize_sample_name(c): c for c in ncols}
            for pc in pcols:
                auto.append({"POS_sample": pc, "NEG_sample": n_norm.get(normalize_sample_name(pc), None)})
            mapping = pd.DataFrame(auto)

            with st.expander("Review/edit POS ↔ NEG sample names", expanded=True):
                mapping = st.data_editor(
                    mapping, use_container_width=True, hide_index=True,
                    column_config={
                        "POS_sample": st.column_config.SelectboxColumn("POS sample", options=pcols),
                        "NEG_sample": st.column_config.SelectboxColumn("NEG sample", options=ncols),
                    },
                    num_rows="dynamic", key="sample_mapping"
                )
            base = st.session_state.get("edges_current", edges)
            corr_edges = correlate_pairs(base, pq, nq, pos_id_col, neg_id_col, mapping)
            st.session_state["edges_current"] = corr_edges
            st.dataframe(corr_edges[["pair_id","pos_node","neg_node","intensity_pearson_r","n_paired_samples"]],
                         use_container_width=True, hide_index=True)
        except Exception as e:
            st.error(f"Correlation step failed: {e}")

with tabs[2]:
    st.subheader("Spectral Pair Viewer")
    st.caption("Independent POS and NEG spectra are shown side by side. No cross-polarity cosine score is calculated.")
    current = st.session_state.get("edges_current", edges)
    if current.empty:
        st.info("No polarity pair is available to inspect.")
    elif not pos_mgf or not neg_mgf:
        st.info("Upload both MGF files in the sidebar to enable the viewer.")
    else:
        ps = parse_mgf(pos_mgf); ns = parse_mgf(neg_mgf)
        pi = spectrum_index(ps); ni = spectrum_index(ns)
        choices = current["pair_id"].astype(str).tolist()
        chosen = st.selectbox("Polarity pair", choices)
        row = current.loc[current["pair_id"].astype(str) == chosen].iloc[0]

        with st.expander("Peak display settings", expanded=True):
            c1, c2, c3 = st.columns(3)
            with c1:
                peak_mode = st.selectbox(
                    "Peaks to display",
                    ["All peaks", "Top 10", "Top 20", "Top 50", "Top 100"],
                    index=2,
                    help="Selects the most intense peaks independently in each spectrum."
                )
            with c2:
                min_rel = st.selectbox(
                    "Minimum relative intensity (%)",
                    [0.0, 1.0, 2.0, 5.0, 10.0],
                    index=1
                )
            with c3:
                label_n = st.selectbox(
                    "Label most intense peaks",
                    [0, 3, 5, 10, 20],
                    index=2
                )

        top_n = None if peak_mode == "All peaks" else int(peak_mode.split()[-1])

        st.write(f"**POS node:** {row.pos_node} · m/z {row.pos_mz:.6f} · RT {row.pos_rt_min:.3f} min  |  "
                 f"**NEG node:** {row.neg_node} · m/z {row.neg_mz:.6f} · RT {row.neg_rt_min:.3f} min")
        a, b = st.columns(2)
        with a:
            sp = find_spectrum(pi, row.pos_node)
            st.plotly_chart(
                spectrum_plot(sp, "Positive-ion MS/MS", top_n=top_n, min_rel=min_rel, label_n=label_n),
                use_container_width=True
            )
            if sp: st.json(sp["meta"], expanded=False)
        with b:
            sn = find_spectrum(ni, row.neg_node)
            st.plotly_chart(
                spectrum_plot(sn, "Negative-ion MS/MS", top_n=top_n, min_rel=min_rel, label_n=label_n),
                use_container_width=True
            )
            if sn: st.json(sn["meta"], expanded=False)
        st.caption("If a spectrum is not found, the GraphML node identifier and MGF FEATURE_ID are not the same identifier in that export. The app deliberately does not guess the mapping.")

with tabs[3]:
    st.subheader("Cytoscape export")
    current = st.session_state.get("edges_current", edges)
    if current.empty:
        st.warning("There are no Polarity Merge edges to export.")
    else:
        keep_all = current[current.get("accepted", True).astype(bool)].copy()
        if used.get("use_rt_distribution", True):
            audited, export_info = classify_rt_distribution(keep_all, used.get("mad_k",3.0))
            keep = audited[audited["rt_status"] == "Inlier"].copy()
            outliers = audited[audited["rt_status"] == "Outlier"].copy()
        else:
            audited = keep_all.copy()
            audited["rt_status"] = "Not filtered"
            keep, outliers = audited, audited.iloc[0:0].copy()

        st.write(f"The merged network will contain **{Gp.number_of_nodes()+Gn.number_of_nodes():,} nodes**, "
                 f"**{Gp.number_of_edges()+Gn.number_of_edges():,} original edges**, and "
                 f"**{len(keep):,} RT-supported PolarityMerge edges**.")
        if len(outliers):
            st.warning(f"**{len(outliers)} RT outlier candidate(s)** will not become PolarityMerge edges, but remain in the audit CSV.")
        st.success("All original POS and NEG molecular-network edges, including their existing cosine/spectral-similarity attributes, are preserved.")
        try:
            GM = merged_graph(Gp, Gn, keep)
            blob = graphml_bytes(GM)
            st.download_button("Download merged Cytoscape GraphML", blob,
                               "polarity_merge_network.graphml", "application/xml")
            st.download_button("Download exported RT-supported edges (CSV)",
                               keep.to_csv(index=False).encode(), "polarity_merge_edges_accepted.csv", "text/csv")
            st.download_button("Download all reviewed candidates + RT status (CSV)",
                               audited.to_csv(index=False).encode(), "polarity_merge_edges_rt_audit.csv", "text/csv")
        except Exception as e:
            st.error(f"Could not create GraphML export: {e}")


with tabs[4]:
    st.subheader("Method")
    st.markdown(r"""
The initial implementation targets the chemically explicit pair **[M+H]⁺ ↔ [M−H]⁻**.

For a positive-ion node:

\[
M_{POS}=m/z_{POS}-1.007276466621
\]

For a negative-ion node:

\[
M_{NEG}=m/z_{NEG}+1.007276466621
\]

A cross-polarity edge is proposed when the neutral-mass error is within the selected MS1 tolerance and:

\[
|RT_{NEG}-RT_{POS}| \leq RT_{tolerance}
\]

For every candidate, the app also retains the **signed** retention-time shift:

\[
\Delta RT = RT_{NEG}-RT_{POS}
\]

Thus, positive ΔRT means the NEG feature eluted after the POS feature.

The initial RT tolerance is a permissive candidate-generation window. The app then estimates the observed cross-polarity offset using the median and MAD of the signed ΔRT distribution. When enabled, a candidate is exported when:

\[
|\Delta RT-\mathrm{median}(\Delta RT)| \leq k\times MAD
\]

The plot displays the median and ±2/±3 MAD guides. The acceptance boundary is inclusive (`≤ k × MAD`) and uses a small numerical tolerance to avoid floating-point boundary artifacts. RT outliers remain in the audit CSV but are excluded from the merged GraphML. No one-to-one constraint is imposed: 1:N and N:1 relationships are intentionally preserved for later investigation.

The expected raw precursor difference for this ion pair is approximately **2.014553 Da**, but matching by reconstructed neutral mass makes the criterion explicit and extensible to other adduct pairs.

If abundance tables are supplied, Pearson correlation across manually reviewed POS↔NEG sample mappings is added as **supporting evidence**, not as an obligatory matching rule.

The exported GraphML prefixes node IDs with `POS::` or `NEG::`, preserving both original networks and preventing identifier collisions. All original within-polarity edges and their existing cosine/spectral-similarity attributes are retained unchanged; cross-polarity links are added as `edge_type = PolarityMerge`.
""")

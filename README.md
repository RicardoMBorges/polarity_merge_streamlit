# Polarity Merge

**Cross-polarity molecular-network integration for LC–MS/MS**  
**LAABio — Biodiversity Analysis and Evaluation Laboratory | IPPN–UFRJ**  
Developed by Ricardo M. Borges

Polarity Merge is a [Streamlit](https://streamlit.io/) application for linking **positive-ion (POS)** and **negative-ion (NEG)** molecular networks using neutral-mass agreement and chromatographic retention-time evidence. It **preserves the original networks** and adds explicit cross-polarity edges that can be explored in [Cytoscape](https://cytoscape.org/).

> **Current implementation: v0.7.** The supported ion-pair hypothesis is **`[M+H]+ ↔ [M−H]−`**. The software generates *candidate associations*, not definitive chemical identifications.

## Contents

1. [Why Polarity Merge?](#why-polarity-merge)
2. [Features](#features)
3. [Installation and deployment](#installation-and-deployment)
4. [Input files](#input-files)
5. [Quick-start tutorial](#quick-start-tutorial)
6. [Detailed matching method](#detailed-matching-method)
7. [Retention-time offset and MAD filtering](#retention-time-offset-and-mad-filtering)
8. [Abundance correlation](#abundance-correlation)
9. [MS/MS spectral viewer](#msms-spectral-viewer)
10. [Exports and Cytoscape](#exports-and-cytoscape)
11. [Interpreting ambiguous matches](#interpreting-ambiguous-matches)
12. [Worked interpretation](#worked-interpretation)
13. [Limitations and quality control](#limitations-and-quality-control)
14. [Troubleshooting](#troubleshooting)
15. [Reproducibility and citation](#reproducibility-and-citation)

## Why Polarity Merge?

ESI positive and negative ionization can reveal complementary portions of the same metabolome. Molecular networking is commonly performed separately for each polarity, producing two graphs whose node IDs and edge scores are not directly comparable.

Polarity Merge asks a narrower, chemically testable question:

> **Could a POS feature and a NEG feature correspond to the same neutral molecule?**

It answers by comparing neutral masses and retention times, while keeping POS and NEG fragmentation networks intact. It does **not** merge nodes, force a one-to-one assignment, or equate fragmentation similarity across ionization modes.

## Features

- Upload independent POS and NEG **GraphML** molecular networks.
- Generate cross-polarity candidates using a defined `[M+H]+ / [M−H]−` mass relation.
- Choose MS1 mass tolerance in **ppm** or **Da**.
- Set an initial, permissive retention-time window.
- Calculate signed `ΔRT = RT_NEG − RT_POS`, in minutes and seconds.
- Characterize the observed ΔRT distribution with **median** and **MAD**.
- Filter candidates using an adjustable **median ± k × MAD** interval (default `k = 3`).
- Visualize retention-time agreement and the signed ΔRT distribution, with median and ±2/±3 MAD guides.
- Review candidate pairs in an editable table.
- Optionally calculate cross-sample **Pearson intensity correlations**, using editable POS↔NEG sample mappings.
- Optionally inspect POS and NEG **MGF** spectra side by side, with Top-N and relative-intensity controls.
- Export a merged **GraphML**, accepted-edge **CSV**, and RT-audit **CSV**.
- Preserve original within-polarity network edges and their stored attributes.
- Display LAABio branding in the sidebar.
- Execute matching only when **Run Polarity Merge** is pressed.

## Installation and deployment

### Requirements

- Python 3.10+ recommended.
- Dependencies declared in `requirements.txt`: `streamlit`, `pandas`, `numpy`, `networkx`, `plotly`, and `openpyxl`.

### Run locally

```bash
git clone https://github.com/YOUR-ACCOUNT/YOUR-REPOSITORY.git
cd YOUR-REPOSITORY

python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# macOS / Linux:
# source .venv/bin/activate

pip install -r requirements.txt
streamlit run app.py
```

Replace the placeholder repository URL with your actual GitHub URL. Streamlit will display a local address, usually `http://localhost:8501`.

### Deploy on Streamlit Community Cloud

1. Create a GitHub repository and upload `app.py`, `requirements.txt`, `README.md`, and `assets/LAABio_logo.png`.
2. Open [share.streamlit.io](https://share.streamlit.io/).
3. Select **Create app** and connect the GitHub repository.
4. Set the main file path to **`app.py`** and deploy.
5. Open the resulting URL and upload your own POS/NEG files.

No database or credentials are required by the current application. Uploaded data are processed in the running Streamlit session; use deployment settings appropriate to your data-confidentiality requirements.

### Repository structure

```text
.
├── app.py
├── requirements.txt
├── README.md
├── .gitignore
└── assets/
    └── LAABio_logo.png
```

Do not commit large, sensitive, or unpublished sample data unless you intend to distribute them.

## Input files

### Required: two GraphML files

| File | Role |
|---|---|
| POS GraphML | Positive-ion molecular network |
| NEG GraphML | Negative-ion molecular network |

Each candidate node must provide numeric `mz` and a retention-time attribute readable by the app (`rt_min` or `rt`), **in minutes**. Nodes missing mass or RT cannot participate in candidate generation. Source GraphML node IDs are treated as identifiers, not as chemically meaningful numbers.

Existing edges may contain attributes such as `score`, `EdgeScore`, `matched_peaks`, `EdgeType`, `deltamz`, and annotations. These are carried forward as original within-polarity edges.

**Important:** Confirm the time units in your source files. MGF headers often use `RTINSECONDS`, while this application's GraphML matching uses **minutes**. Do not mix the two.

### Optional: two MGF files

| File | Role |
|---|---|
| POS MGF | Positive-ion MS/MS spectra |
| NEG MGF | Negative-ion MS/MS spectra |

The viewer uses MGF spectrum identifiers (including `FEATURE_ID` and related fields) to locate spectra associated with GraphML nodes. Ensure that node IDs and MGF feature IDs are compatible. A missing spectrum does not invalidate an MS1/RT candidate.

### Optional: abundance matrices

Upload a POS and a NEG feature-by-sample table (`.csv`, `.tsv`, `.txt`, `.xlsx`, or `.xls` where supported). Each table should have:

- One feature-ID column, whose values match the corresponding network's node IDs.
- One numeric intensity column per sample.
- Corresponding biological samples across polarities, even if their column names differ.

Example POS:

| FeatureID | Sample_A_POS | Sample_B_POS | Sample_C_POS |
|---|---:|---:|---:|
| 1721 | 150000 | 220000 | 180000 |
| 5427 | 80000 | 65000 | 99000 |

Example NEG:

| FeatureID | Sample_A_NEG | Sample_B_NEG | Sample_C_NEG |
|---|---:|---:|---:|
| 85 | 140000 | 210000 | 170000 |
| 186 | 76000 | 62000 | 95000 |

These values are **illustrative**, not experimental results.

**A GraphML network and an MGF file do not substitute for full feature-by-sample abundance tables.** A single MGF `FEATURE_MS1_HEIGHT` is not a sample-level abundance profile.

## Quick-start tutorial

### 1. Upload the networks

In the sidebar, upload **Positive GraphML** and **Negative GraphML**. You may also upload POS/NEG MGF files for spectral inspection.

### 2. Set candidate-generation tolerances

Start with:

| Control | Suggested starting value | Meaning |
|---|---:|---|
| RT tolerance | `0.100 min` | Initial absolute POS/NEG RT difference ≤ 6 s |
| MS1 tolerance | `5 ppm` | Maximum neutral-mass disagreement |
| RT distribution filter | On | Evaluate the observed offset |
| Acceptance band | `±3 MAD` | Final RT-supported interval |

These are starting points, **not universal analytical acceptance criteria**. The initial RT window must be wide enough to include the actual cross-polarity offset and enough candidates to estimate its distribution. If the observed shift is larger than six seconds, enlarge the initial window before running.

### 3. Press **Run Polarity Merge**

The app does **not** run automatically when files are uploaded or settings are changed. Click **Run Polarity Merge** to calculate candidates.

The overview shows POS nodes, NEG nodes, original network edges, and candidate cross-polarity edges.

### 4. Review **Polarity matches**

Inspect the candidate table, including:

- POS and NEG node identifiers.
- Precursor m/z values and reconstructed neutral masses.
- Signed and absolute mass error.
- POS RT, NEG RT, and signed/absolute ΔRT.
- `pair_id` and the editable `accepted` field.

Use the acceptance checkbox to exclude clearly unsuitable candidates. This manual review is separate from the RT distribution filter.

### 5. Explore the ΔRT distribution

In the same tab, inspect:

- Median signed ΔRT.
- Median absolute deviation (MAD).
- Mean absolute ΔRT and its 95th percentile.
- POS RT versus NEG RT agreement plot.
- Signed ΔRT histogram, with offset and ±2/±3 MAD guides.
- Inlier/Outlier classifications.

A **negative** ΔRT means NEG eluted **earlier** than POS. A **positive** ΔRT means NEG eluted **later**.

### 6. Optionally inspect abundance correlations

Enable **Use cross-polarity intensity correlation** and upload the POS and NEG quantification tables. Choose feature-ID columns, review the proposed POS↔NEG sample-name mapping, and inspect Pearson `r` and the number of paired samples.

The mapping table is editable; verify it carefully. Do not pair technical injections or biological samples incorrectly.

### 7. Optionally inspect MS/MS spectra

Open **Spectral Pair Viewer**, select a candidate pair, and inspect its POS and NEG spectra side by side. Use:

- **Peaks to display:** All / Top 10 / Top 20 / Top 50 / Top 100.
- **Minimum relative intensity:** 0, 1, 2, 5, or 10%.
- **Label most intense peaks:** 0, 3, 5, 10, or 20.

Each polarity is plotted independently. This is **not a mirror plot**, and the app does **not** compute POS-versus-NEG spectral cosine.

### 8. Export

Open **Export** and download:

1. `polarity_merge_network.graphml` — original POS and NEG networks plus RT-supported cross-polarity edges.
2. `polarity_merge_edges_accepted.csv` — cross-polarity pairs included in the exported network.
3. `polarity_merge_edges_rt_audit.csv` — reviewed candidates with RT classification, including excluded outliers.

**Note:** The RT-audit export reflects the currently reviewed/accepted candidate subset. To retain the full initial candidate list, also download **Download candidate edges (CSV)** from **Polarity matches** before manual exclusions.

## Detailed matching method

### Neutral-mass reconstruction

The initial implementation assumes that the POS feature is a protonated molecule and the NEG feature is a deprotonated molecule:

\[
[M+H]^+ \longleftrightarrow [M-H]^-
\]

Using the proton mass:

\[
m_p = 1.007276466621\ \mathrm{Da}
\]

the neutral masses are reconstructed as:

\[
M_{POS} = (m/z)_{POS} - m_p
\]

\[
M_{NEG} = (m/z)_{NEG} + m_p
\]

and:

\[
\Delta M = M_{NEG} - M_{POS}
\]

The signed error in ppm is:

\[
\mathrm{error}_{ppm}
= \frac{\Delta M}{|M_{POS}|}\times 10^6
\]

The candidate passes the MS1 criterion if its **absolute** error is within the selected ppm threshold, or if `|ΔM|` is within the selected Da threshold.

For this specific ion-pair hypothesis, the expected *raw* precursor difference is:

\[
(m/z)_{POS}-(m/z)_{NEG}
\approx 2m_p = 2.014552933242\ \mathrm{Da}
\]

The app evaluates **neutral-mass agreement**, not merely whether the raw difference resembles 2.01455 Da.

**Chemical caution:** This model does not automatically validate ion/adduct identities. Sodium, ammonium, formate, chloride, multimers, in-source fragments, isotopologues, or incorrect feature/adduct assignments may create misleading associations. Current matching is specific to `[M+H]+` versus `[M−H]−`.

### Initial chromatographic window

Candidates must also satisfy:

\[
|RT_{NEG}-RT_{POS}| \le RT_{tolerance}
\]

RT values are in minutes. This is the **candidate-generation** filter; it is not the final offset-aware criterion.

## Retention-time offset and MAD filtering

Positive and negative LC–MS acquisitions may exhibit a small systematic retention-time offset. Consequently, a pair with ΔRT near the **observed** offset can be more plausible than one with ΔRT near zero.

Define signed retention-time shift:

\[
\Delta RT_i = RT_{NEG,i}-RT_{POS,i}
\]

The robust center is:

\[
c = \mathrm{median}(\Delta RT_i)
\]

The median absolute deviation is:

\[
MAD = \mathrm{median}(|\Delta RT_i-c|)
\]

For a user-selected multiplier `k`, the accepted interval is:

\[
c-k\,MAD \le \Delta RT_i \le c+k\,MAD
\]

or equivalently:

\[
|\Delta RT_i-c| \le k\,MAD
\]

The boundary is **inclusive**; v0.7 includes a small numerical tolerance to avoid excluding points that fall exactly on the theoretical boundary because of floating-point rounding.

### Why median and MAD?

Mean and standard deviation can be distorted by a small number of mismatched features. Median and MAD provide robust descriptive estimates for a distribution dominated by plausible pairs.

The app uses **raw MAD**, not `1.4826 × MAD` (the commonly used normal-consistency-scaled MAD). Thus, `3 MAD` here must **not** be interpreted automatically as a classical `3σ` interval or as a calibrated false-discovery threshold.

### Practical implications

- A systematic shift is **retained**, not forced to zero.
- A candidate far from the observed offset is marked `Outlier`.
- An `Outlier` is omitted from the **new cross-polarity edges** in the GraphML.
- Original POS and NEG nodes and original within-polarity edges remain untouched.
- The excluded pair is retained in the RT-audit CSV.
- Candidate multiplicity is not automatically resolved.

### Important implementation details

- The initial candidate generation still uses the **absolute** RT window around zero. Therefore the window must be broad enough to capture the offset; the MAD step cannot recover candidates excluded earlier.
- The distribution is estimated from **candidate edges**, not from a set of verified unique chemical identities. Multiple candidates from the same node can influence the estimated offset.
- The RT distribution is recalculated from the currently **manually accepted** candidate subset in the review/export workflow. Manual exclusions may therefore change the estimated center and MAD.
- The app does **not** perform repeated automatic outlier-removal iterations.
- If fewer than three finite ΔRT values are available, the classifier returns `Insufficient data`; no RT-supported edges are exported with filtering enabled. In that case, review the data or disable the distribution filter deliberately.
- If `MAD = 0`, the classifier treats values effectively equal to the median as inliers. Inspect tied/discretized RT values carefully.
- The default ±3 MAD is an exploratory heuristic, not a validated statistical confidence interval.

### Recommended QC practice

Before interpreting the offset as an instrument or acquisition effect, inspect the RT distribution across the chromatographic gradient, evaluate repeated acquisitions and QC samples, and check whether candidate redundancy is distorting the distribution. For high-confidence use, estimate the offset from independently verified POS↔NEG pairs or QC standards and then apply it to the full network.

## Abundance correlation

If paired POS and NEG feature-by-sample matrices are available, the app can calculate Pearson correlation for each candidate:

\[
r_{POS,NEG} =
\mathrm{corr}(\mathbf{x}_{POS},\mathbf{x}_{NEG})
\]

where the vectors contain intensities for the **same mapped samples** across polarities.

This is **supporting evidence only**. The current workflow does not require a minimum `r` to create a cross-polarity edge. Correlation may be affected by batch effects, missing values, signal saturation, ion suppression, and differences in ionization efficiency. Verify sample alignment and preprocessing before interpretation.

A high `r` does not prove identical chemical identity; a low `r` does not necessarily disprove it.

## MS/MS spectral viewer

The viewer is designed for **qualitative cross-polarity inspection**.

- POS and NEG MS/MS spectra appear **side by side**.
- Intensities are displayed relative to each spectrum's own base peak.
- Top-N and intensity filters are applied independently to each spectrum.
- Peak labels aid manual comparison.
- No cross-polarity cosine or mirror-spectrum assumption is imposed.

Positive- and negative-ion fragmentation pathways can differ substantially. Existing within-polarity molecular-network cosine or other spectral similarity attributes remain part of their original edges; they are **not** recomputed as cross-polarity similarity.

## Exports and Cytoscape

### Merged network design

To prevent collisions between POS and NEG GraphML node IDs, the exported graph prefixes node IDs:

```text
POS::1721
NEG::85
```

Original IDs are retained in the `original_node_id` node attribute, and polarity is stored as `polarity = POS` or `NEG`.

Within-polarity edges are retained with their original attributes and assigned `edge_type = MolecularNetworking` where that lowercase attribute was absent. Added cross-polarity edges use `edge_type = PolarityMerge`.

Representative cross-polarity attributes include:

| Attribute | Meaning |
|---|---|
| `pair_id` | Candidate pair identifier |
| `source`, `target` | POS/NEG prefixed node IDs (in CSV) |
| `pos_node`, `neg_node` | Original feature IDs |
| `pos_mz`, `neg_mz` | Precursor m/z values |
| `neutral_mass_pos`, `neutral_mass_neg` | Reconstructed neutral masses |
| `delta_neutral_mass_da` | Signed neutral-mass difference |
| `mass_error_ppm` | Signed ppm error |
| `delta_rt_sec` | Signed RT difference |
| `abs_delta_rt_sec` | Absolute RT difference |
| `rt_center_sec`, `rt_mad_sec` | Distribution parameters |
| `rt_deviation_sec`, `rt_robust_z` | Offset-centered RT deviation and deviation/raw MAD |
| `rt_status` | `Inlier`, `Outlier`, or another applicable status |
| `edge_type` | `PolarityMerge` |

Optional `intensity_pearson_r` and `n_paired_samples` are included when computed.

### Import into Cytoscape

1. Launch Cytoscape.
2. Use **File → Import → Network from File**.
3. Select `polarity_merge_network.graphml`.
4. In the Style panel, map:
   - **Node fill color** → `polarity` (POS vs NEG).
   - **Node label** → `original_node_id`.
   - **Edge color / line type** → `edge_type`.
5. Use different visual treatments for `MolecularNetworking` and `PolarityMerge`.
6. Inspect components and ambiguous 1:N / N:1 links without collapsing nodes.

A separately supplied **`Polarity_Merge_Cytoscape_Style.xml`** can also be imported through Cytoscape's style-import controls. This XML is **not automatically bundled** into the v0.7 app export; add it to the repository separately if you wish to distribute it.

**Graph preservation caveat:** The current exporter builds a NetworkX `Graph` (not a `MultiGraph`). It is appropriate for ordinary simple molecular networks; if an input contains parallel edges between the same two nodes, those parallel edges may not be preserved separately. Validate edge counts for multigraph inputs.

## Interpreting ambiguous matches

A POS feature may match several NEG features, or vice versa. The application deliberately preserves these **1:N** and **N:1** cases.

Potential explanations include:

- Peak splitting during chromatographic processing.
- Different peak integration/deconvolution between polarities.
- Closely eluting isomers or partially resolved features.
- Redundant features or imperfect adduct assignment.
- Accidental neutral-mass and RT agreement.

**Do not automatically choose the smallest ppm error as the only valid edge.** Examine chromatographic peak shapes, raw extracted-ion chromatograms, MS/MS evidence, and sample-level abundance behavior first.

## Worked interpretation

**Illustrative example**, not a benchmark: assume the candidate ΔRT distribution has median `−4.5 s` and raw MAD `2.7 s`. With `k = 3`:

\[
-4.5 - 3(2.7) \le \Delta RT \le -4.5 + 3(2.7)
\]

so the RT-supported interval is:

\[
\boxed{-12.6\ \mathrm{s}\ \le \Delta RT\ \le +3.6\ \mathrm{s}}
\]

A pair with `ΔRT = −5.0 s` is close to the observed offset and passes. A pair with `ΔRT = +7.2 s` falls outside and is flagged. A pair exactly at `+3.6 s` passes the **inclusive** boundary.

**Note:** With the app's default initial RT tolerance of 0.100 min (6 s), a `+7.2 s` pair would **not enter candidate generation** in the first place. To explore such candidates, increase the initial RT tolerance (for example to 0.15 min) before running. This illustrates why the two RT filters must be distinguished.

## Limitations and quality control

1. **Candidate matching, not identification.** Matching m/z and RT is insufficient to establish structural identity.
2. **Ion-pair assumption.** Only `[M+H]+ ↔ [M−H]−` is currently modeled.
3. **No automatic RT alignment.** The offset is estimated for filtering; RT values in the input graphs are not transformed or warped.
4. **Potential circularity.** The offset is estimated from the candidates being filtered; an independent reference set would be stronger.
5. **Multiplicity can bias the distribution.** Several candidate edges may originate from one chromatographic feature.
6. **MAD is descriptive, not a probability model.** Do not interpret the selected cutoff as a formal p-value or FDR.
7. **No mandatory correlation cutoff.** Abundance correlation is optional supporting evidence.
8. **MGF mapping depends on IDs.** Inconsistent GraphML/MGF feature IDs may prevent spectrum retrieval.
9. **Original edges are not cross-polarity cosine evidence.** Existing within-mode similarity attributes retain their original meaning.
10. **Potential multigraph limitation.** Parallel edges may be collapsed by the simple-graph exporter.
11. **Session-state behavior.** Changing settings does not automatically rerun matching; press **Run Polarity Merge** again.
12. **Export after manual review.** Save the initial candidate table if you need a full record of candidates rejected manually.

## Troubleshooting

| Symptom | Check / action |
|---|---|
| No candidates found | Confirm RT units, `mz` fields, ion/adduct hypothesis, and mass tolerance; widen the initial RT window if needed. |
| RT offset looks truncated | The initial ±RT window may be too narrow to capture the full shift distribution. |
| Many 1:N pairs | Review chromatographic peak splitting, integration, coelution, and duplicate features. |
| No spectrum in viewer | Verify GraphML node ID ↔ MGF `FEATURE_ID` correspondence. |
| Correlation unavailable | Provide full POS and NEG abundance matrices and verify sample mappings. |
| `Insufficient data` | Fewer than three finite candidate ΔRT values are available for the filter. |
| RT inliers unexpectedly change | Manual candidate acceptance may alter the center/MAD used at export. |
| Cytoscape colors not applied | Map the lowercase `edge_type` and node `polarity` columns; import the optional style XML separately. |
| Updated parameters have no effect | Press **Run Polarity Merge** to recompute. |
| GraphML export fails | Check unsupported/non-scalar attributes and graph type; inspect the Streamlit error message. |

## Reproducibility and citation

For a reproducible analysis, record:

- App version or Git commit.
- Source GraphML and optional MGF/abundance filenames and hashes.
- Ion-pair hypothesis.
- Mass tolerance (ppm or Da).
- Initial RT window.
- RT filter enabled/disabled, multiplier `k`, median ΔRT, and MAD.
- Number of candidates, inliers, outliers, and manually excluded pairs.
- POS↔NEG sample mapping and preprocessing, if abundance correlation is used.
- Final GraphML and both CSV exports.

### Suggested methods wording

> Positive- and negative-ion molecular networks were integrated using Polarity Merge (LAABio, IPPN–UFRJ). Candidate cross-polarity links were generated by comparing neutral masses reconstructed under the `[M+H]+` and `[M−H]−` ion-pair hypothesis, subject to predefined MS1 mass-error and initial retention-time tolerances. The signed retention-time offset (`RT_NEG − RT_POS`) was characterized by its median and median absolute deviation (MAD). When enabled, candidate links outside the median ± k × MAD interval were excluded from the cross-polarity export. Original within-polarity nodes, edges, and their similarity attributes were retained, and multiple candidate correspondences were not forcibly resolved. Candidate assignments were treated as putative associations rather than confirmed metabolite identifications.

**Software citation:** Until a versioned release, DOI, or publication exists, cite the GitHub repository URL, author, and exact release/commit used. Do not invent a DOI.

### License

No software license is declared here. Add a `LICENSE` file to the repository after choosing appropriate terms; without an explicit license, standard copyright restrictions generally apply.

---

**LAABio — Biodiversity Analysis and Evaluation Laboratory**  
**Instituto de Pesquisas de Produtos Naturais (IPPN), Universidade Federal do Rio de Janeiro (UFRJ)**  
**Polarity Merge | v0.7**

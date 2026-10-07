# Polarity Merge

Streamlit app for connecting positive- and negative-ion molecular networks without collapsing nodes.

## Core workflow

1. Upload POS and NEG Cytoscape GraphML files.
2. Choose the mass/RT tolerances and explicitly press **Run Polarity Merge**; uploading files or changing controls does not automatically recompute the merge.
3. Reconstruct neutral masses for the initial ion pair `[M+H]+ / [M-H]-`.
4. Match candidates by MS1 mass tolerance and retention-time tolerance.
5. Review/accept candidate cross-polarity edges.
6. Optionally upload POS/NEG abundance tables, review sample-name mapping, and calculate cross-polarity Pearson correlations.
7. Optionally upload POS/NEG MGF files and inspect the two MS/MS spectra side by side, with Top-N peak filtering, a minimum relative-intensity threshold, and selective m/z labels.
8. Export a Cytoscape-ready GraphML containing both original networks plus the new `PolarityMerge` edges.

## Important design choices

- Original POS and NEG nodes are preserved; they are **never collapsed**.
- Original molecular-network edges are preserved, including existing cosine/spectral-similarity edge attributes.
- POS and NEG node IDs receive `POS::` and `NEG::` prefixes only in the merged export to avoid ID collisions.
- No POS-vs-NEG MS/MS cosine similarity is calculated.
- Intensity correlation is optional supporting evidence.
- The app does not silently infer an MGF-to-GraphML feature-ID mapping when the identifiers differ.

## Local run

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Streamlit Community Cloud

Push `app.py` and `requirements.txt` to a GitHub repository. In Streamlit Community Cloud, create a new app from that repository and select `app.py` as the entry point.

## Input assumptions

The GraphML nodes should contain precursor mass in `mz` and retention time in `rt_min` or `rt`. The current implementation assumes retention time is in minutes.

For optional abundance correlation, upload a feature × sample table for each polarity. One column must contain the feature/node identifier; the remaining selected numeric columns are sample abundances.

## Current scope

Version 0.1 implements `[M+H]+ ↔ [M-H]-`. The matching engine is intentionally structured so that additional adduct-pair rules can be added later.


## Retention-time QC

For each POS↔NEG candidate, the app exports both signed and absolute retention-time differences:

- `delta_rt_min` / `delta_rt_sec`: `RT_NEG - RT_POS`
- `abs_delta_rt_min` / `abs_delta_rt_sec`: absolute difference

The interface reports median signed ΔRT, MAD, mean absolute ΔRT, 95th percentile absolute ΔRT, an RT POS vs RT NEG scatter plot with a 1:1 reference, and a signed ΔRT distribution.


## Robust observed ΔRT filtering
The initial RT tolerance generates candidates. Final chromatographic support is evaluated around the observed signed offset `ΔRT = RT_NEG - RT_POS`, using `median ± k×MAD` (default k=3). The app shows median, ±2 MAD and ±3 MAD guides. RT outliers are excluded from GraphML cross-polarity edges but retained in an audit CSV. 1:N/N:1 relationships are not collapsed.

The `±k×MAD` criterion is inclusive. A small numerical epsilon is applied so candidates exactly on the theoretical boundary are not incorrectly rejected by floating-point representation.

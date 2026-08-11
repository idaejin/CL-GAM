# pennLC nested ATA — Case A / B / C

**Status:** PRELIMINARY (2026-08-11)  
**Data:** [`SpatialEpi::pennLC`](https://rdrr.io/cran/SpatialEpi/man/pennLC.html) — PA lung cancer 2002, 67 counties + smoking  
**Scripts:** `R/07_load_pennLC.R`, `scripts/17_penn_case_ABC.R`, `scripts/18_plot_penn_case_ABC.R`

## Design

| | |
|--|--|
| Fine | 67 counties (centroids, stratified expected via `expected(..., n.strata=16)`) |
| Coarse | 20 k-means regions of centroids (`seed=1`) — not admin districts |
| **A** | `s(smoking)` at county |
| **B** | linear `north_agg` = region-mean latitude expanded |
| **C** | A + B |

Why not “all NY”: no CRAN bundle with statewide NY tracts+counties comparable to NY8. pennLC gives **67 fine units** and **20 coarse** \(y\), enough for stable `pois_SAP` and a non-flat county-oracle scatter.

## Outputs

- `output/penn_case_ABC_fit.rds`, `penn_case_ABC_summary.csv`
- `output/penn_case_ABC_maps.png`, `penn_case_ABC_diagnostics.png`
- Competitor metrics: `scripts/20_compare_pennLC.R` → `penn_compare_metrics.csv` / `.png`

## Competitor ranking (full suite, PRELIMINARY)

Script: `scripts/20_compare_pennLC.R` → `penn_compare_metrics.csv` / `.png`

| Layer | What we see |
|-------|-------------|
| Coarse deviance / loglik | Baseline saturates (~0); Case C ≲ A ≲ B ≈ ATA ≪ Malone |
| AIC (`pois_SAP` only) | ATA **22.2** &lt; Case A 23.4 &lt; B 24.2 &lt; C 25.2 |
| Mass max/L2 | Baseline ~0; CL-GAMM/ATA ~23–24 / ~48–53; Malone **64 / 96** |
| Fine MSE/MAE vs oracle | All CL-GAMM ≈ ATA ≈ 0.045; Malone worse; baseline 0.044 |
| cover95 vs oracle | 0.40–0.49 (`pois_SAP`); Malone 0.25 — oracle is noisy (not truth) |
| LOO coarse log-score | Baseline best; ATA/Case A next; Malone worst |
| Truth MSE / cover | **NA** (no known η) |

Structural: within-region oracle variance ~74%; cor ceiling ≈ 0.51.

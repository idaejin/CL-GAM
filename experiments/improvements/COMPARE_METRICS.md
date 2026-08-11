# Model comparison metrics

**Status:** superseded as narrative by the **DECISION** note  
→ **`notes_metrics_ATA.md`** (canonical) · Brain [[ATA-evaluation-metrics]]

**Scripts:** `19_compare_metrics.R` (Madrid+penn short) · `20_compare_pennLC.R` (full suite) · helper `R/08_metrics.R`

## Primary metric (cross-method)

\[
\mathrm{MSE}_\eta = \frac{1}{m}\sum_{j=1}^{m}\bigl(\hat\eta_j - \log((y_j+0.5)/e_j)\bigr)^2
\]

at the **fine** support, where \(y_j,e_j\) are fine counts/expected **not used as response** in ATA fits (oracle / hold-in descriptive). Lower is better.

## Secondary

| Metric | Role |
|--------|------|
| `cor_eta` | association with oracle log-rate |
| `rmse_counts` | \(\sqrt{\mathrm{mean}((e e^{\hat\eta}-y_{\mathrm{fine}})^2)}\) |
| `mass_err` | \(\max\|C\hat\mu - y_{\mathrm{coarse}}\|\) (calibration of aggregates) |
| `aic` | **only** among `pois_SAP` on the same dataset |

## Competitors included

**Madrid:** ATA spatial, CL-GAMM Case A, Malone Method 1, ST-PCLM ATA (if RDS present), mun-SMR baseline.  
**pennLC:** ATA spatial, Case A/B/C, Malone, region-SMR baseline.

`disaggregation` (Madrid raster) deferred: needs CT aggregation of cell predictions for the same `mse_eta`.

## Outputs

- `output/compare_metrics.csv` / `.rds` / `.png`

## First run (2026-08-11) — interpretation

**Factual (tables):** On both Madrid and pennLC, the **region-SMR expanded** baseline has the lowest `mse_eta`. Among likelihood-based models, Madrid ranks ATA spatial &lt; Malone ≈ Case A ≪ ST-PCLM ATA on `mse_eta`; Case A still wins **AIC** among `pois_SAP` (390.9 vs ATA 396.9). pennLC: Case B / ATA ≈ Case C / Case A ≪ Malone on `mse_eta`; spatial-only has best AIC.

**Inferred:** Low fine MSE ≠ best coarse AIC — Case A improves the composite-link fit while the spatial field + covariates need not track noisy fine log-SMR. The SMR-expansion baseline is a strong descriptive competitor when within-unit heterogeneity is modest; it is **not** a smooth covariate-aware risk surface.

**Caveat:** `mse_eta` uses fine counts as an oracle diagnostic; those counts are not the ATA response. Do not treat baseline “wins” as evidence against CL-GAMM’s scientific goal (covariate + spatial smooth under mass constraint).

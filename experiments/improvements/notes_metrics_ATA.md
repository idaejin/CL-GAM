# Evaluation metrics for CL-GAMM (ATA)

**Status:** **DECISION** (2026-08-11) — hierarchy for simulations, real ATA, and competitor comparison.  
**Provenance:** user specification + pennLC/Madrid scorecards; coherent with SMiMR `PAPER - Simulations - version 3` (`mse(eta.est, eta.true)`).

Canonical helper: `experiments/R/08_metrics.R`  
Scripts: `scripts/19_compare_metrics.R` (Madrid+penn), `scripts/20_compare_pennLC.R` (full suite)  
Brain: [[ATA-evaluation-metrics]] (in `00_RESEARCH_BRAIN/methods/`)

---

## Hierarchy

### 1. Known fine truth \(\eta\) (simulations) — **primary for the paper / `clgam`**

| Metric | Role |
|--------|------|
| **MSE / MAE of \(\hat\eta\) vs \(\eta_{\mathrm{true}}\)** (fine support) | **Primary verdict** |
| MSE / MAE of \(\exp(\hat\eta)\) vs \(\exp(\eta_{\mathrm{true}})\) | Optional risk scale |
| Coverage of ±1.96 SE bands for \(\eta\) vs **truth** | Uncertainty calibration |
| Bias by strata (urban/rural, covariate extremes) | Optional diagnostics |

Paper v3: `mse(eta.est, eta.true)` over 100 Poisson replicates; boxplots vs Malone / multistep.

**Do not** treat `cor(\hat\eta,\eta)` as primary (amplitude/shape can fail while cor looks OK).

For additive \(g(z)\): do **not** trust raw `cor(\hat g, g)` if \(g\) is correlated with \(z\) (e.g. \(\sin(2\pi z)\)). Use MSE(\(\hat g\)) and/or correlation of residuals after regressing out linear \(z\).

Code fields (when truth known): `mse_eta_truth`, `mae_eta_truth`, `cover95_truth`.

### 2. Only aggregated data (real ATA) — **primary on Madrid / open ATA toys**

Comparable across methods that deliver fine \(\hat\mu\) and can form \(C\hat\mu\):

| Metric | Role |
|--------|------|
| Coarse Poisson **log-score / deviance** of \(y\) vs \(C\hat\mu\) | Fit under composite link (fair across engines) |
| **AIC / BIC / ED** | **Only within** the same likelihood family (e.g. `pois_SAP` Case A/B/C vs ATA spatial) |
| **Mass calibration** \(\max_i\|(C\hat\mu)_i - y_i\|\) (also \(\|C\hat\mu-y\|_2\)) | Whether aggregates recover observed coarse counts |
| Leave-one-coarse-unit-out predictive log-score | Preferable if feasible (pennLC: IDW holdout η approx.) |

**Mass calibration (definition).** Only coarse \(y\) is observed. Fine means \(\hat\mu_j=e_j\exp(\hat\eta_j)\). Region totals are \((C\hat\mu)_i=\sum_{j\in i}\hat\mu_j\). Residual of mass in region \(i\): \((C\hat\mu)_i-y_i\). Report

\[
\max_i\bigl|(C\hat\mu)_i-y_i\bigr|
\qquad\text{and/or}\qquad
\|C\hat\mu-y\|_2.
\]

≈0 ⇒ predicted totals match observed coarse counts (baseline SMR expansion saturates by construction). Large ⇒ map does not conserve mass (typical of Malone-style heuristics without \(C\mu\approx y\) in the likelihood). Does **not** measure within-region risk quality.

**Do not** compare AIC Malone (GAM on disaggregated proxies) vs `pois_SAP` (CLMM).

### 3. Fine counts as diagnosis only (Madrid MEDEA / pennLC)

Secondary; structural ceiling (noisy oracle; change-of-support):

| Metric | Role |
|--------|------|
| MSE/MAE vs oracle \(\log((y_j+0.5)/e_j)\) | Rough fine check |
| `cor_eta` | **Supplement only** (pennLC ceiling ≈0.51; ~74% oracle var within region) |
| RMSE / MAE of fine counts | Secondary |
| `cover95_oracle` using `sd.eta` | Diagnostic only — oracle ≠ truth; expect miscalibration |

---

## Implications for simulations (`clgam::simulate_ata`, paper v3)

1. Generate smooth \(\eta_{\mathrm{true}}\) → fine (or latent) means → \(y=C y_f\) (or \(y\sim\mathrm{Poisson}(C\gamma)\)).
2. Monte Carlo: **fix \(\eta_{\mathrm{true}}\)**, replicate only counts; primary table = MSE(\(\hat\eta\)).
3. Report SE coverage vs truth when `elements=TRUE` / `sd.eta` available.
4. Correlations only as supplements; recovery of nonlinear \(g(z)\) is a **separate** software/identifiability experiment, not the SMiMR sim claim.
5. Align competitor re-runs with the same truth and the same fine grid / \(C\).

## Implications for fitted models (competitors)

| Method class | What to report |
|--------------|----------------|
| `pois_SAP` / CL-GAMM A–C / ATA spatial | deviance, AIC within family, mass, LOO if cheap; fine oracle secondary |
| Malone / dissever heuristics | deviance & mass under delivered \(\hat\mu\) (scoring rule); **not** their GAM AIC vs SAP |
| Baseline region-SMR expanded | mass ≈0, deviance ≈0 (saturated coarse); useful descriptive foil, not a smooth cov surface |
| Grid `disaggregation` | aggregate cell preds to fine areal support before MSE/mass |

## pennLC scorecard *(PRELIMINARY, 2026-08-11)*

See `output/penn_compare_metrics.csv`. Summary: baseline saturates coarse; among `pois_SAP`, ATA best AIC; Case A best fine RMSE counts; Malone worst on deviance, mass, LOO, fine MSE. Truth metrics NA.

## Code pointers

- Helper: `experiments/R/08_metrics.R` — `clgam_score_ata()`, `clgam_oracle_lograte()`, `clgam_baseline_region_smr()`
- Full pennLC: `scripts/20_compare_pennLC.R`
- Madrid+penn short: `scripts/19_compare_metrics.R`
- Paper sims: `Diego Ayma/SMiMR/.../PAPER - Simulations - version 3.R`
- Package: `clgam::simulate_ata()`; upgrade demos away from cor-only (`clgam_recovery_caseA.png`)

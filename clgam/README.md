# `clgam` R package

**Version 0.1.25.** Composite link **generalized additive (mixed) models** for areal counts:
disaggregate coarse Poisson (or quasi-Poisson) observations to a nested fine support via a composition
matrix \(C\), with anisotropic spatial P-splines and optional fine-scale smooth
covariates. Estimation uses PIRLS + separation of overlapping precision matrices (SOP).

This package is a clean public rewrite of the experimental `spclmm` fork used in
the CL-GAM / Ayma line of work. **No MEDEA or other proprietary health data are
included.** Examples and tests use [`simulate_ata()`](R/simulate_ata.R).
Call `simulate_ata_scenarios()` for every named DGP and the recommended formula.

## Install

```r
remotes::install_github("idaejin/CL-GAM", subdir = "clgam")
# from a local clone of this package directory:
install.packages(".", repos = NULL, type = "source")
```

## Quick start

```r
library(clgam)

simulate_ata_scenarios()[, c("scenario", "formula", "orth.smooth")]

dat <- simulate_ata(scenario = "A", n_coarse = 16, n_fine_per = 8, seed = 1)
# dat$sf_coarse / dat$sf_fine are nested Voronoi polygons (requires sf)

# Formula interface (y is coarse; s() covariates are fine).
# Bare z is linear; s(z) is a P-spline. See ?clgam::s (not mgcv::s).
# P-spline settings: ndx (intervals), bdeg (degree, default 3),
# pord (penalty order, default 2). k is an alias for ndx.
fit <- clgam(
  y ~ s(x1, x2, ndx = c(10, 10), bdeg = 3, pord = 2) +
    s(z_f, ndx = 12),
  data = dat,
  orth.smooth = TRUE
)
# Optional: Laplace (requires TMB) and/or coarse iid overdispersion
# clgam(y ~ s(x1, x2) + s(z_f), data = dat, method = "Laplace")
# clgam(y ~ s(x1, x2) + s(z_f), data = dat, re = "coarse")

# Same model with C and exposure named:
# clgam(y ~ s(x1, x2) + s(z), C = C, exposure = ef, orth.smooth = TRUE)

# Positional API (still supported)
fit_pos <- clgam(
  y = dat$y,
  coords = cbind(dat$x1, dat$x2),
  C = dat$C,
  exposure = dat$efine,
  smooth = dat$nlcovfine,
  knots = c(10, 10),
  knots_nl = 12,
  orth.smooth = TRUE
)
# Overdispersion: same fitted surface, SEs scaled by sqrt(phi)
fit_qp <- clgam(
  y = dat$y,
  coords = cbind(dat$x1, dat$x2),
  C = dat$C,
  exposure = dat$efine,
  smooth = dat$nlcovfine,
  family = quasipoisson()
)
fit_qp$phi

plot(fit, which = 1, sf_fine = dat$sf_fine, sf_coarse = dat$sf_coarse)
plot(fit, which = 4, g_true = dat$g_true)
```

Canonical low-level fitters: `pois_SOP` / `pois_incat_SOP`
(legacy aliases `pois_SAP` / `pois_incat_SAP` remain available).
Two-group contrasts: `clgam_contrast()` (wraps `pois_incat_SOP`).
(Optional Laplace/TMB coverage benchmark lives under
`experiments/benchmark_tmb/` — not part of this package.)

S3 methods: `print`, `summary`, `coef`, `fitted`, `residuals`, `predict`,
`plot`, `AIC`, `BIC`, `logLik`, `nobs`, `deviance`.

```r
fitted(fit, type = "eta")
residuals(fit)
predict(fit, type = "eta", se.fit = TRUE)
plot(fit, which = 1:3)   # eta map, obs vs fit, residuals
```

Toggle compiled SOP kernels with `options(clgam.use_rcpp = FALSE)` (default: compiled TRUE).

## Evaluation (simulations and models)

Evaluation hierarchy (consistent with the CL-GAMM/ATA scoring notes in
`experiments/`):

1) Simulations (truth known) — *primary*

When \(\eta_{\mathrm{true}}\) is known (`simulate_ata()`), report
**MSE/MAE of \(\hat\eta\)** on the fine support.
If `elements=TRUE` is available, also report nominal 95% coverage of the
pointwise \(\hat\eta \pm 1.96\,SE\) band (vs truth), using the
**unconditional** SE that accounts for variance-component uncertainty
(Wood–Pya–Säfken; `experiments/R/10_se_unconditional.R`). The default
`fit$sd.eta` is Bayesian conditional on \(\hat\tau^2\).

Do not use correlation with a noisy fine oracle (e.g. SMR) as the primary
verdict: `cor()` can look acceptable while amplitude/shape are wrong.

2) Real ATA / fitted models (only coarse \(y\) observed) — *primary*

Models are compared by their fit to the **aggregated** Poisson likelihood under
the same likelihood family. Use:

- Coarse Poisson deviance / log-score for \(y\) vs \(C\hat\mu\).
- **Mass calibration**: \(\max_i |(C\hat\mu)_i - y_i|\) (and/or \(\|C\hat\mu-y\|_2\)).
- AIC only within the same composite-link engine (e.g. `pois_SOP` variants).

3) Fine diagnostics (Madrid / pennLC)

If fine oracle information is available only for diagnosis, use MSE/MAE of
fine log-rate vs oracle as a secondary check. Treat it as descriptive, not a
formal ranking criterion.

Full hierarchy: `../experiments/improvements/notes_metrics_ATA.md`
(Research Brain: `ATA-evaluation-metrics`).

## Relation to papers

- Spatial CLMM (ATP/ATA): Ayma, Durbán, Lee, Eilers — *Spatial Statistics* (2016)
- Unpublished CL-GAMM (multi-resolution covariates): working draft in the CL-GAM project

Method name in papers: **CL-GAM**. R package name: **`clgam`**.

For the “ATP via nested ATA” narrative used in the experiments, see
`experiments/improvements/ST_PCLM_ATA.md` (renamed to “ATP (area-to-point)” in
plots/tables, while the underlying estimator is still `pois_SOP`).

## License

GPL (>= 2)

# `clgam` R package

**Version 0.1.25.** Composite-link generalized additive (mixed) models for areal
counts: disaggregate coarse Poisson (or quasi-Poisson) observations to a nested
fine support via a composition matrix `C`, with anisotropic spatial P-splines
and optional multi-resolution smooth covariates (Cases A--C). Estimation is
PIRLS + separation of overlapping precision matrices (SOP). Optional Laplace
(TMB) uses the same mixed-model design.

This package is a public rewrite of the experimental `spclmm` fork used in the
CL-GAM / Ayma line of work. **No MEDEA or other proprietary health data are
included.** Examples and tests use [`simulate_ata()`](R/simulate_ata.R).
Call `simulate_ata_scenarios()` for every named DGP and the recommended formula.

## Install

```r
remotes::install_github("idaejin/CL-GAM", subdir = "clgam")
# from this package directory:
install.packages(".", repos = NULL, type = "source")
```

`method = "Laplace"` needs [TMB](https://github.com/kaskr/adcomp) (`Suggests`).
Voronoi geometry in `simulate_ata()` needs `sf`.

## Quick start

```r
library(clgam)

simulate_ata_scenarios()[, c("scenario", "formula", "orth.smooth")]

dat <- simulate_ata(scenario = "A", n_coarse = 16, n_fine_per = 8, seed = 1)
# dat$sf_coarse / dat$sf_fine are nested Voronoi polygons (requires sf)

# Formula: y is coarse; s() covariates are fine.
# Bare z is linear; s(z) is a P-spline. See ?clgam::s (not mgcv::s).
# ndx = internal B-spline intervals; bdeg = degree (default 3);
# pord = penalty order (default 2). k is an alias for ndx.
fit <- clgam(
  y ~ s(x1, x2, ndx = c(10, 10), bdeg = 3, pord = 2) +
    s(z_f, ndx = 12),
  data = dat,
  orth.smooth = TRUE
)
# Optional: Laplace (TMB) and/or coarse iid overdispersion
# clgam(y ~ s(x1, x2) + s(z_f), data = dat, method = "Laplace")
# clgam(y ~ s(x1, x2) + s(z_f), data = dat, re = "coarse")

plot(fit, which = 1, sf_fine = dat$sf_fine, sf_coarse = dat$sf_coarse)
plot(fit, which = 4, g_true = dat$g_true)
```

Cases B/C: `s(z_a, level = "coarse")`. Identified A/C simulations use
`identifying = "perp"` (default) so the spatial truth is `f_perp`, matching
`orth.smooth = TRUE`. Spatial coordinates `x1`, `x2` are used as supplied
(a small range-relative pad for the knot domain); rescale to `[0, 1]` yourself
if the two axes have very different units.

Low-level fitters: `pois_SOP` / `pois_incat_SOP` (aliases `pois_SAP` /
`pois_incat_SAP`). Two-group models: `clgam_contrast()`; pointwise difference
inference: `clgam_contrast_infer(fit)`.

S3 methods: `print`, `summary`, `coef`, `fitted`, `residuals`, `predict`,
`plot`, `AIC`, `BIC`, `logLik`, `nobs`, `deviance`.

```r
fitted(fit, type = "eta")
predict(fit, type = "eta", se.fit = TRUE)
plot(fit, which = 1:3)
```

Toggle compiled SOP kernels with `options(clgam.use_rcpp = FALSE)`
(default: compiled `TRUE`).

## Evaluation

When the fine truth is known (`simulate_ata()`), the primary score is MSE/MAE
of `eta_hat` vs `eta_true`. Coverage of pointwise 95% bands should use
unconditional SEs (Wood--Pya--Safken), not `fit$sd.eta` (Bayesian, conditional
on `tau^2`). Do not rank methods by correlation with a noisy fine oracle.

When only coarse `y` is observed, compare models on the aggregated Poisson
likelihood: coarse deviance / log-score, mass calibration
`max |C mu_hat - y|`, and AIC only within the same composite-link engine.

## Papers

- Spatial CLMM (ATP/ATA): Ayma, Durban, Lee, Eilers, *Spatial Statistics* (2016)
- Unpublished CL-GAMM (multi-resolution covariates): working draft

Method name in papers: **CL-GAM**. R package name: **`clgam`**.

## License

GPL-2

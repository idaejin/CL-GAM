# clgam

Composite link **generalized additive mixed models** (CL-GAMM) for areal counts:
disaggregate coarse Poisson observations to a nested fine support via a composition
matrix \(C\), with anisotropic spatial P-splines and optional fine-scale smooth
covariates. Estimation uses PIRLS + separation of penalties (SOP/SAP).

This package is a clean public rewrite of the experimental `spclmm` fork used in
the CL-GAM / Ayma line of work. **No MEDEA or other proprietary health data are
included.** Examples and tests use [`simulate_ata()`](R/simulate_ata.R).

## Install

```r
# remotes::install_github("idaejin/clgam")  # when published
# from a local clone:
install.packages(".", repos = NULL, type = "source")
```

## Quick start

```r
library(clgam)

dat <- simulate_ata(
  n_coarse = 16, n_fine_per = 8, seed = 1,
  include_covariate = TRUE
)
# dat$sf_coarse / dat$sf_fine are nested Voronoi polygons (requires sf)

fit <- clgam(
  y = dat$y,
  coords = cbind(dat$x1, dat$x2),  # fine centroids
  C = dat$C,
  exposure = dat$efine,
  smooth = dat$nlcovfine,
  knots = c(10, 10),
  knots_nl = 12
)

plot(fit, which = 1, sf_fine = dat$sf_fine, sf_coarse = dat$sf_coarse)
plot(fit, which = 4, g_true = dat$g_true)
```

Legacy paper names `pois_SAP` / `pois_incat_SAP` remain available.
Two-group contrasts: `clgam_contrast()` (alias of `pois_incat_SAP`).

S3 methods: `print`, `summary`, `coef`, `fitted`, `residuals`, `predict`,
`plot`, `AIC`, `BIC`, `logLik`, `nobs`, `deviance`.

```r
fitted(fit, type = "eta")
residuals(fit)
predict(fit, type = "eta", se.fit = TRUE)
plot(fit, which = 1:3)   # eta map, obs vs fit, residuals
```

Toggle compiled SOP kernels with `options(clgam.use_rcpp = FALSE)`.

## Evaluation (simulations and models)

When \(\eta_{\mathrm{true}}\) is known (`simulate_ata()`), report **MSE/MAE of \(\hat\eta\)**
first (optional ±1.96 SE coverage vs truth). On real ATA with only coarse \(y\),
use coarse Poisson deviance/log-score, mass calibration \(\max\|C\hat\mu-y\|\),
and AIC only among fits that share the same likelihood. Do not treat
correlation with a noisy fine log-SMR as the primary score, and do not compare
AIC across heuristic Malone-style GAMs vs SOP CLMM.

Full hierarchy: `../experiments/improvements/notes_metrics_ATA.md`
(Research Brain: ATA-evaluation-metrics).

## Relation to papers

- Spatial CLMM (ATP/ATA): Ayma, Durbán, Lee, Eilers — *Spatial Statistics* (2016)
- Unpublished CL-GAMM (multi-resolution covariates): working draft in the CL-GAM project

Method name in papers: **CL-GAMM**. R package name: **`clgam`**.

## License

GPL (>= 2)

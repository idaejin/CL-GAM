# clgam 0.1.23

* Formula interface: `clgam(y ~ s(x1, x2) + s(z), C = C, exposure = ef,
  orth.smooth = TRUE)`. `s()` is parsed by the package (not `mgcv`).
  `data=` accepts a `simulate_ata()` list (`C` / `efine` filled in when
  omitted). Case B/C: `s(z_a, level = "coarse")` expands a coarse
  covariate through partition `C`.
* `simulate_ata_scenarios()` catalogues every DGP preset (default, A, B,
  C, confounding, matern1, matern2, jump) with recommended formula and
  `orth.smooth`. Documented in `?simulate_ata`.
* `R CMD check --as-cran` is clean of ERROR/WARNING (Rd math in
  `\eqn{}`, ASCII sources, documented S3 arguments).

# clgam 0.1.22

* `simulate_ata()` misspecification DGPs whose spatial truth is not an
  anisotropic P-spline: Matérn Gaussian fields (`scenario="matern1"` /
  `"matern2"`, \(\nu=1,2\)) and a piecewise-constant jump between two
  coarse-unit groups (`scenario="jump"`). Also `spatial_truth=`,
  `matern_nu=`, `matern_range=`, `jump_amp=`.

# clgam 0.1.21

* `simulate_ata()` can reproduce the manuscript Monte Carlo with one call:
  `scenario="A"|"B"|"C"|"confounding"`, covariate-truth presets
  (`nl_fun="sine"|"tanh"|...` or a custom function), spatial confounding
  `rho`, an `intercept` / `exposure_scale` for the Poisson intensity, and
  `family="poisson"` (default) or `"none"` (mean counts, no Poisson noise).
  Existing defaults (no `scenario`) are unchanged.

# clgam 0.1.20

* `summary.clgam()` is reviewer-facing: Pearson `phi` (and SE scaling
  under quasi-Poisson), overlap `kappa` by fine-scale covariate, SOP
  `tau^2` / ED labelled by smooth, AIC/BIC, and PIRLS convergence
  (`niter`, last relative change in `eta`, `diverged`). Fits store
  `tol`, `maxit`, `thr`, `converged`, and `orth.info$kappa_by`.
  Contrast summaries print the joint AIC/BIC (not the per-group vector).

# clgam 0.1.19

* Export `kappa_diagnostic()`: functional-space overlap
  \(\kappa=\|P_{A_f} f\|^2/\|f\|^2\) (identifiability diagnostic, manuscript
  eq. 21). Fits store `orth.info$kappa` for the fitted spatial field. For the
  paper's design overlap, pass `f` (unrestricted spatial field or known
  \(f_{\mathrm{raw}}\)).

# clgam 0.1.18

* Kronecker–Schur W3: `options(clgam.sop.backend = "kron_hybrid")`
  inverts the inner SOP Schur complement `S` by an exact B3-vs-rest
  block factorization (same estimator as dense `solve(S)`). Default
  remains `"dense"` (Rcpp when compiled). Contrast models
  (`pois_incat_SOP`) do not use the hybrid path.

# clgam 0.1.17

* Kronecker–Schur W2: isolated B3 solver `.sop_solve_b3()` (`dense` /
  Jacobi-`pcg` / diagnostic `diag`) and `.sop_solve_P2D_diag()` for the
  original-tensor `sparse_P2D + W` path. Production `pois_SOP` Schur
  solve is unchanged.

# clgam 0.1.16

* Internal Kronecker–Schur W1 helpers (`.sop_kron_meta()`,
  `.sop_ginv_spatial()`, `.sop_N_blocks()`): spatial B1/B2/B3 index sets
  matching `pois_SOP`, with tests that `N[idx_B3, idx_B3]` equals the
  dense Gram of the interaction block. The inner Schur solve is unchanged.

# clgam 0.1.15

* When `elements=TRUE`, cache the PIRLS-weighted Bayesian blocks (`M1`, `M2`)
  at the final outer PIRLS step and reuse them for AIC/BIC/SEs instead of
  recomputing `clmm_mat` + `inv_bblock2` after convergence (`pois_SOP`,
  `pois_incat_SOP`).

# clgam 0.1.14

* Fix Schur SOP dispatch: `.sop_solve_schur()` now detects
  `sop_solve_schur_cpp` with `exists(..., mode = "function")` (same as
  `clmm_mat()`). The previous `envir = environment(), inherits = FALSE`
  check always failed, so the compiled Schur solver was never used despite
  `options(clgam.use_rcpp = TRUE)`.

# clgam 0.1.13

* Add `family = poisson()` / `family = quasipoisson()` to `clgam()`,
  `clgam_contrast()`, `pois_SOP()`, and `pois_incat_SOP()`.
  Quasi-Poisson keeps Poisson PIRLS+SOP point estimates and variance
  components; a Pearson dispersion `fit$phi` is estimated at convergence
  and all standard errors are multiplied by `sqrt(phi)`. Fit objects now
  store GLM `family` (`"poisson"` / `"quasipoisson"`) and model `type`
  (`"spatial"` / `"contrast"`); legacy `family = "spatial"|"contrast"`
  in `.as_clgam()` remains accepted.
* `summary()` reports the Pearson dispersion and notes SE scaling under
  quasi-Poisson.

# clgam 0.1.12

* Correctness and robustness from code review: AIC/BIC no longer
  double-count shared `ed` on contrast fits; non-finite variance-component
  updates abort with `warning()` and `fit$diverged`; `pois_incat_SOP()`
  validates matching fine dimensions; linear-covariate SEs
  (`sdleffects`) use the same PIRLS-weighted Bayesian covariance as
  nonlinear SEs when `elements=TRUE`.
* Performance: sparse `Diagonal` for `Ginv`; `mm_basis()` uses
  `eigen(..., symmetric=TRUE)`; SOP Schur solver caches G-free blocks
  across inner iterations and uses a general solve for `S` (not
  `likely_sympd`; `S` is not symmetric in general).
* New regression tests in `tests/testthat/test-fixes.R`.

# clgam 0.1.11

* Identifiability projection (`orth.smooth`) now uses the space of
  **fine-scale** covariate effects only:
  raw fine B-spline bases \(B(z)\) plus optional linear fine covariates
  (`lcovfine`). Coarse-scale smooths are excluded via `nl.level` /
  `clgam(..., smooth_level=)`.
* Case B: pass `smooth_level="coarse"`; Case C:
  `smooth_level=c("fine","coarse")`. `simulate_ata()` returns matching
  `nl_level`.
* Fit objects store `orth.info` with numerical orthogonality diagnostics.
* SOP penalty blocks and `.sop_solve_schur()` unchanged.

# clgam 0.1.10

* Remove optional `pois_TMB()` from the package. Laplace composite-link
  benchmark lives under `experiments/benchmark_tmb/` only; public engine
  remains PIRLS + SOP.

# clgam 0.1.9

* Rename the PIRLS+SOP engine to `pois_SOP` / `pois_incat_SOP`
  (separation of overlapping precision matrices). Legacy names
  `pois_SAP` / `pois_incat_SAP` remain as aliases.
* Rename Rcpp Schur kernel to `sop_solve_schur_cpp`.
* Document SOP estimation in DESCRIPTION and README.
* Optional `pois_TMB()` Laplace check (Suggests: TMB) — moved out in 0.1.10.
* Generate Rd help pages via roxygen2.

# clgam 0.1.8

* Previous release (SOP kernels under historical SAP names).

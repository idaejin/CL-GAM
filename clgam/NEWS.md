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

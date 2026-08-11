# Package `clgam` — public CL-GAMM software (no MEDEA)

**DECISION (2026-08-11):** new package **`clgam`** (not rename of archival `spclmm`); method name in papers remains CL-GAMM.

## Location
`01_PROJECTS/CL-GAM/clgam/` — source tree ready for GitHub (`idaejin/clgam`).

## Rules
- **No MEDEA / proprietary health data** in the repo.
- Examples and tests use `simulate_ata()`.
- Madrid replication stays under `experiments/` (local only).

## Seeded from
`experiments/pkg/spclmm` 0.2.0 core: `pois_SAP`, `pois_incat_SAP`, Rcpp SOP kernels, sparse `C` helpers, `pois_PCLM_latent`. Dropped choromap/isomap/REML variants for a lean first public surface.

## Public API (cleaned 2026-08-11)
- `simulate_ata()`, `clgam()`, `clgam_contrast()`
- Legacy: `pois_SAP`, `pois_incat_SAP`
- S3: `print` / `summary` / `coef` / `fitted` / `residuals` / `predict` / `plot` / `AIC` / `BIC` / `logLik` / `nobs` / `deviance`
- Fit banner only when `trace=TRUE`
- Explicit NAMESPACE exports
- `predict(..., newdata=)` not yet (fitted grid only)

## Covariate smooth \(g(z)\) *(Factual / DECISION 2026-08-11, updated)*

- **`clgam(..., smooth = Z)`** Case A: univariate P-spline (`nl.basis="pspline"`).
  Null space via `mm_basis(..., decom=2)` **without intercept** + penalized \(\mathbf{Z}_k\).
- **`orth.smooth=TRUE`** (pspline default): spatial bases projected ⊥ \(\mathrm{span}\{B(z)\}\)
  so the spatial field cannot absorb the additive smooth (fixes linear-only \(\hat g\)).
- **`pois_SAP`** default remains **`nl.basis="legacy"`**, `orth.smooth=FALSE` (Madrid/SMiMR).
- Variance components floored at `1e-50`.
- Recovery demos: need enough coarse units for joint \(f(s)+g(z)\) recovery under ATA.
  Example Case A: `n_coarse=40`, tuned amps → `clgam_recovery_caseA.png`.
  Example Case B: `covariate_level="coarse"` → `clgam_recovery_caseB.png`
  (`μ=exp(log(Cγ)+h(z_a))` via piecewise-constant expansion of \(z_a\)).
  Example Case C: `covariate_level="both"` → `clgam_recovery_caseC.png`
  (fine \(g(z_f)\) + coarse \(h(z_a)\); scripts `21_` / `22_`).

**Generative story (correct):**
1. Smooth \(\eta(s)\) on fine centroids (piecewise constant on subpolygons)
2. \(\gamma = e_f\exp(\eta)\), \(e_f \propto\) area
3. \(y_f\sim\mathrm{Poisson}(\gamma)\)
4. Raw map \(y=C y_f\)

**Geometry:** nested Voronoi partitions (exact `n_fine_per` per coarse); `C` 0-1; coverage checked (empty symdiff = OK).

**Fixes:** no longer drop Voronoi cells; no buffer fallback; fixed `n_fine_per`; empty-gap NA bug; Rcpp Schur falls back to R on SVD failure.

## NEXT
- `git init` + push to GitHub when approved
- roxygen man pages; CITATION
- optional vignette with only simulated data
- wire `experiments/R` to prefer installed `clgam` once stable

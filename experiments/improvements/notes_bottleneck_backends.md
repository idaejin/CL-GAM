# Bottleneck backends: SparseMatrix vs LMMsolver vs Schur-SOP

## Sparse backends in `spclmm` 0.1.8

- `as_comp_C(C, backend=)` / `pois_SAP(..., sparse.backend=)`:
  `"Matrix"` (default), `"spam"`, `"dense"`, or `"auto"` /
  `options(spclmm.sparse.backend=...)`.
- Madrid-style **partition** `C` (one coarse cell per fine cell) uses
  `rowsum` aggregation (backend-agnostic; usually fastest).
- `sparse_P2D()` / `sparse_chol_solve()`: Kronecker anisotropic precision
  for a future sparse-CLMM path (Matrix or spam).

Script: `scripts/13_sparse_backend_bench.R`.

After PIRLS weights, Ayma SOP repeatedly inverts a **dense** nonsymmetric system
`(V+D)` of size `(p+q)≈529` (ndx=20). That matrix is dense because
`mm_basis` SVD-reparameterizes B-splines into dense `Z`.

Wrapping `V+D` in `Matrix::sparseMatrix` / `spam` does **not** help: density ≈ 1.

## Three ways forward

### A. Schur SOP (in `spclmm` 0.1.7) — same estimator

Cache `solve(XtX)` (tiny) and invert only the Schur complement of size `q`;
recover ED weights from block rows of `(V+D)^{-1}` without a full `(p+q)` inverse.
**Same algebra / same AIC targets.** Script: re-run `08_speed_check.R`.

### B. Sparse Kronecker precision (SparseMatrix) — reformulation

Keep coefficients on the tensor B-spline basis with
`P = λ1 (I⊗D'D) + λ2 (D'D⊗I)` ([Boer 2023](https://doi.org/10.1177/1471082X231178591)).
Then `Matrix::Cholesky` is fast (see `scripts/09_sparse_precision_demo.R`).
**Requires rewriting CLMM** (composite-link weights + sparse normal equations).
This is the route that actually matches what [LMMsolver](https://biometris.github.io/LMMsolver/) does internally.

### C. LMMsolver — fast for `C = I`, not for ATA/CLMM

[`LMMsolve`](https://biometris.github.io/LMMsolver/) + `spl2D()` uses sparse REML
P-splines ([site](https://biometris.github.io/LMMsolver/)).

**Empirical (Madrid mun, ndx/nseg=20):** with **standardized coordinates**,
`family=poisson()`, `offset=log(e)`, `lmm$yhat` matches `pois_SAP` η to
corr ≈ 1 (RMSE ~1e-4) in **~1.3 s vs ~24 s** (~20×). Raw UTM coords → Cholesky singularity.

**No composition `C`:** cannot express Case A mun→ct / sex contrast without a custom
composite-link layer. Path B (sparse Kronecker + `C`) is how to merge the ideas.

Script: `scripts/10_lmmsolver_mun_baseline.R`.

## Recommended programme

| Priority | Action | Claim status |
|----------|--------|--------------|
| Now | Schur in fork (A) — same AIC; modest vs ED cost | ESTABLISHED |
| Now | LMMsolver mun `C=I` baseline (~1 s vs ~24 s SOP; corr η high if coords scaled) | PRELIMINARY |
| NEXT | Sparse CLMM (B): PIRLS + Kronecker `P` + composition `C` | HYPOTHESIS |

Do not cite LMMsolver fits as CL-GAMM Case A / ATA results.

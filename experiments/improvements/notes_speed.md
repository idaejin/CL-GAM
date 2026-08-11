# Speed improvements (spclmm experiment fork)

## Provenance
Factual timings on Madrid MEDEA (ndx=20,20, `elements=TRUE`), Apple Silicon, `OMP_NUM_THREADS=1`.

| Fit | Baseline (pre-tune) | R Schur / sparse C | **0.2.0 Rcpp** | Match paper? |
|-----|---------------------|--------------------|----------------|--------------|
| model0 mun | ~35–37 s | ~23–24 s | **~2.9 s** | AIC / τ² exact |
| model1 mun→ct | ~95 s | ~66–69 s | **~6.8 s** | AIC / τ² exact |
| sex `pois_incat_SAP` | ~158 s | ~118–120 s | **~12 s** | τ² + dif range exact |

Rcpp vs pure-R on same install (`scripts/15_rcpp_bench.R` + sex check): **~8× model0, ~10× model1, ~10× sex**; AIC / τ² / dif match paper.

Toggle: `options(spclmm.use_rcpp = FALSE)` for pure-R fallback.

## What changed (same estimates)

### Through 0.1.x (R)
1. Sparse `C` / `rowsum` partition; working `z` without nested `t(t(C)*…)`.
2. `sd.eta` / `sd.dif`: `rowSums((A%*%S)*A)`.
3. Schur SOP + early-stop on relative Δτ.
4. Matrix/spam backends: same AIC, **no meaningful further speedup** (aggregation already cheap).

### 0.2.0 (RcppArmadillo)
Hot path in `src/clmm_kernels.cpp`:
- `sap_solve_schur_cpp` — Schur complement + ED diagonal
- `comp_mul_groups_cpp` — partition aggregation
- `clmm_crossprod_cpp` — working Gram blocks
- `btWb_cpp` — Camarda latent B′WB

**Why Rcpp over sparse-basis reformulation (Boer / LMMsolver-style):** Madrid fine grid is scattered; tensor B-spline bases stay effectively dense in the SOP Gram `B′WB`. Sparse *precision* `P` helps storage (`sparse_P2D`) but does not sparsify the PIRLS/SOP cross-products that dominate. Compiled dense BLAS on the Schur path is the direct win without changing the estimator.

Sparse-basis / LMMsolver remains useful for **identity-C** municipal fits (see `10_*`); not a drop-in for ATA composition `C`.

## Remaining
- Optional: sparse-basis SOP rewrite only if a structured fine grid + local support is introduced (different design).

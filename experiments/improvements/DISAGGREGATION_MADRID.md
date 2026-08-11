# Madrid → disaggregation bridge

**Status:** PRELIMINARY scaffold in `experiments/` (2026-08-11).

## What it does

Adapts [`disaggregation`](https://doi.org/10.18637/jss.v106.i11) (TMB Poisson disaggregation + spatial field) to MEDEA Madrid:

| Input | Source |
|---|---|
| Response polygons | Municipality `sf`, `y = ym`, joined by `GEOCODIGO` (tidy order ≠ shapefile order) |
| Covariates | CT `unemployed`, `ageing` rasterized (`touches=TRUE`) |
| Aggregation raster | Expected deaths per cell = `ec / n_cells_in_CT` |

## Scripts

- `R/04_disaggregation.R` — helpers
- `scripts/06_disaggregation_madrid.R` — fit + RDS/PDF

```r
Sys.setenv(CLGAM_FAST = "1")  # 200 m grid, 80 iters, no iid
# full: 100 m, 300 iters, iid=TRUE
Rscript scripts/06_disaggregation_madrid.R
```

Package install (once):

```r
install.packages("disaggregation", lib = "experiments/R_libs")
# also needs terra, sf, TMB, fmesher (pulled as deps)
```

## Caveats *(Factual / Inferred)*

- Covariates are **areal CT painted on a grid**, not true environmental rasters.
- `disaggregation` uses **linear** covariate effects + GP; not P-spline \(g_k\).
- Tiny CTs with no intersecting cell claim the nearest empty cell (`n_empty_painted`).
- Prefer **simulation MSE** for method comparison; Madrid run is a pipeline/illustration.

## Do not

- Treat FAST smoke AIC/maps as manuscript evidence.
- Claim identical estimand to CL-GAMM Case A without stating the raster bridge.

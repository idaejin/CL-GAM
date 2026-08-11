# NY leukemia — CL-GAMM Case A / B / C (reproducible)

**Status:** PRELIMINARY illustration (not MEDEA).  
**Script:** `scripts/15_ny_case_ABC.R` · loader `R/06_load_ny_leukemia.R`  
**Data:** `SpatialEpi::NYleukemia` + `DClusterm::NY8` (CRAN)

## Default: true counties + poor bases

| | |
|--|--|
| \(\mathbf{C}\) | **8 counties × 281 tracts** (FIPS) |
| Spatial | `ndx = (3,3)` (Madrid-style `20` is singular here) |
| **A** | linear `PCTOWNHOME` at tract |
| **B** | linear county-mean `PEXPOSURE` expanded to tracts |
| **C** | A + B |

Nonlinear fine smooths / two fine linears (ageing+ownhome) often hit `pinv(): svd failed` with only 8 county observations.

```bash
cd experiments
Rscript scripts/15_ny_case_ABC.R
```

## Optional: k-means mid-level (richer smooths)

```bash
CLGAM_NY_COARSE=kmeans Rscript scripts/15_ny_case_ABC.R
```

Uses 40 k-means aggregates so `s(ageing)+s(ownhome)` is numerically feasible. Documented as synthetic mid-level, not administrative counties.

## Nested Case B equivalence *(Factual for partition C)*

Expanding a coarse covariate as constant within each coarse unit and putting it in `lcovfine` matches manuscript Case B when \(C\) is a partition.

## Outputs

- `output/ny_case_ABC_fit.rds` (`mode`, `table`, `C`, …)
- `output/ny_case_ABC_summary.csv`

# CL-GAM experiments

Replicate (and where useful, modernise) the Madrid MEDEA analyses from the CL-GAM / SMiMR draft, without duplicating the multi-GB `Diego Ayma/` trees.

## Layout

| Path | Role |
|------|------|
| `R/00_paths.R` | Absolute paths to tidy data, cartography, local `spclmm` |
| `R/01_load_spclmm.R` | Install-from-source (once) + `library(spclmm)` |
| `R/02_load_madrid.R` | Load mun/ct counts, covariates, composition matrices, maps (`sf` → `sp` for `choromap`) |
| `R/03_malone.R` | Corrected Malone / dissever ATA helper (`clgam_malone_fit`) |
| `R/06_load_ny_leukemia.R` | Public NY leukemia nested ATA (`SpatialEpi` + `DClusterm::NY8`) |
| `scripts/15_ny_case_ABC.R` | Case A / B / C `pois_SAP` on NY (see `improvements/NY_CASE_ABC.md`) |
| `R/07_load_pennLC.R` | `SpatialEpi::pennLC` nested ATA (67 counties, k-means regions) |
| `scripts/17_penn_case_ABC.R` | pennLC Case A/B/C `pois_SAP` |
| `scripts/18_plot_penn_case_ABC.R` | pennLC maps + diagnostics |
| `scripts/00_smoke_test.R` | Data + package sanity check (seconds) |
| `scripts/01_replicate_ATA_spatial.R` | Municipal GAM vs CLMM downscaling (paper Codes2 models 0–1) |
| `scripts/02_replicate_CaseA_covariates.R` | Case A: ageing + unemployment at CT (paper model 2.5) |
| `scripts/03_replicate_sex_contrast.R` | Optional: sex contrast scaffold (Codes1) |
| `scripts/04_mgcv_explore_covariates.R` | Oracle `mgcv` screening at CT (no `C`) |
| `scripts/05_malone_corrected.R` | Malone Method 1 (paper defaults) + order-safe mass-balance |
| `scripts/06_disaggregation_madrid.R` | `disaggregation` competitor on Madrid (raster bridge) |
| `scripts/07_malone_intermediates.R` | Method 1 vs revive heuristic intermediates (no `C` in lik) |
| `scripts/08_speed_check.R` | Paper AIC/τ² + timing for sped-up `spclmm` |
| `scripts/09_sparse_precision_demo.R` | Kronecker sparse `P` vs dense (SparseMatrix path) |
| `scripts/10_lmmsolver_mun_baseline.R` | [LMMsolver](https://biometris.github.io/LMMsolver/) vs SAP for **C=I** |
| `scripts/11_st_pclm_ata_nested.R` | ST-PCLM spatial ATA nested + raster \(C_{mun,grid}=C_{mun,ct}C_{ct,grid}\) |
| `scripts/12_plot_madrid_comparison.R` | Maps + partials PDF/PNG (ATA, Case A, Malone, oracle) |
| `scripts/13_sparse_backend_bench.R` | Composition backends: dense / Matrix / spam |
| `scripts/14_camarda_latent_pclm.R` | Camarda–Durbán latent PCLM vs `pois_SAP` ([arXiv:2412.04956](https://arxiv.org/abs/2412.04956)) |
| `improvements/notes_camarda_2412.04956.md` | What was ported / not ported from that paper |
| `improvements/notes_bottleneck_backends.md` | Schur / SparseMatrix / spam / LMMsolver map |
| `improvements/MALONE_ADAPTATIONS.md` | Malone 2012 vs Method 1 adaptations vs bugfixes |
| `improvements/DISAGGREGATION_MADRID.md` | Notes on CT→raster bridge + caveats |
| `improvements/ST_PCLM_ATA.md` | Lee et al. 2022 → nested ATA / raster bridge |
| `improvements/` | Notes on package fixes and API clean-ups |
| `output/` | Figures / RDS from runs |

**Data source (read-only):**  
`../Diego Ayma/SMiMR/Community of Madrid data analysis/`

**Package source (fork, installable on modern R):**  
`pkg/spclmm/` (v0.1.9 — `pois_PCLM_latent`, Matrix/spam, Schur SOP; ATP in `R-legacy/`)

**Archival original:**  
`../Diego Ayma/SMiMR/spclmm/` (v0.1.0 — requires retired `rgeos`)

## Quick start (R)

```r
setwd(".../01_PROJECTS/CL-GAM/experiments")  # or open this folder as project
source("scripts/00_smoke_test.R")
```

Full paper-like fits use `ndx = c(20, 20)` and can take **several minutes** each. For a fast dry-run:

```r
Sys.setenv(CLGAM_FAST = "1")   # smaller ndx inside scripts
source("scripts/01_replicate_ATA_spatial.R")
```

## What is *not* copied here

- Raw/sim trees under `Diego/` and duplicate zips — use SMiMR hub only.
- Simulation Monte Carlo (Codes4 / `eta.sim.*`) — wire later when needed.
- `maptools` / Windows `setwd` from original PAPER scripts — replaced by helpers.

## Provenance

- Scripts rewrite paper logic against the **same tidy files** used in SMiMR.
- Estimator remains `pois_SAP` / `pois_incat_SAP` (SOP/PIRLS), consistent with `NOTES_ESTIMATION.md` and manuscript V1.

# ST-PCLM → nested ATA (+ raster)

**Status:** DECISION scaffold in `experiments/` (2026-08-11).  
**Ref:** Lee, Durbán, Ayma, van de Kassteele (2022) PLoS ONE e0263711 (ST-PCLM).

## Adaptation *(report)*

| Lee et al. 2022 (Q-fever) | Madrid nested ATA |
|---|---|
| \(\mathbf{C}_s\): mun → **grid** (ATP) | \(\mathbf{C}_s\): mun → **CT** (ATA 0/1) |
| \(\mathbf{C}_t\): month → week | Not used in CVD tidy (no time) |
| \(\mathbf{C}_{st}=\mathbf{C}_t\otimes\mathbf{C}_s\) | Spatial section only (= CL-GAMM engine) |
| Fine coords: cell centroids | CT centroids (ATA) or cell centroids (raster path) |

**Factual:** nested ATA spatial fit **is** the spatial block of ST-PCLM with a different \(\mathbf{C}_s\); estimator `pois_SAP` matches the GLMM/SOP line used in the PCLM papers.

## Raster nested bridge *(new helper, not a new likelihood)*

Using the disaggregation CT→raster paint:

\[
\mathbf{C}_{\mathrm{mun,grid}}
=
\mathbf{C}_{\mathrm{mun,ct}}\,
\mathbf{C}_{\mathrm{ct,grid}}
\]

- `C_mun_ct`: existing `spC_mun-ct`
- `C_ct_grid`: cell ∈ CT (from `ct_index` raster)
- Fits:
  - **ATA:** `pois_SAP(..., C = C_mun_ct)` — Case A CL-GAMM / ST-PCLM spatial
  - **ATP-via-nested (optional):** `pois_SAP(..., C = C_mun_grid)` — PLoS-style fine grid but built as nested product; set `CLGAM_FIT_GRID=1`

Post-aggregate grid \(\boldsymbol{\mu}\) to CT with `C_ct_grid %*% mu_grid`.

## Code

| File | Role |
|---|---|
| `R/05_st_pclm_ata.R` | `clgam_C_st`, `clgam_stpclm_ata_raster_C`, `clgam_stpclm_ata_fit` |
| `scripts/11_st_pclm_ata_nested.R` | Madrid run |
| `R/04_disaggregation.R` | returns `ct_index` for \(\mathbf{C}_{ct,grid}\) |

```r
Sys.setenv(CLGAM_FAST = "1")
Rscript scripts/11_st_pclm_ata_nested.R
# optional heavier grid fit:
Sys.setenv(CLGAM_FAST = "1", CLGAM_FIT_GRID = "1")
Rscript scripts/11_st_pclm_ata_nested.R
```

## NEXT (not this scaffold)

- Real \(\mathbf{C}_t\) when MEDEA / other data have time.
- Full `sclm::spt_sclm` path for ST (Q-fever-style) with nested \(\mathbf{C}_s\).
- Do not claim Madrid CVD run is spatio-temporal ST-PCLM.

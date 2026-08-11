#!/usr/bin/env Rscript
# ST-PCLM adapted to nested ATA (Lee et al. 2022 composition with C_s = mun→CT),
# plus nested raster bridge C_mun_grid = C_mun_ct %*% C_ct_grid.
#
# Madrid CVD tidy has no time → fit spatial section only (pois_SAP).
# Optional: FIT_GRID=1 also fits mun→grid PCLM (ATP-via-nested; can be heavy).
#
# From experiments/:
#   Sys.setenv(CLGAM_FAST = "1")
#   Rscript scripts/11_st_pclm_ata_nested.R
#
# See improvements/ST_PCLM_ATA.md

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))
source(file.path(root, "R/04_disaggregation.R"))
source(file.path(root, "R/05_st_pclm_ata.R"))

dat <- clgam_load_madrid(load_maps = TRUE)
ndx <- CLGAM_NDX_SPATIAL
ndxnl <- CLGAM_NDX_NL
res_m <- if (CLGAM_FAST) 500 else 200
fit_grid <- identical(Sys.getenv("CLGAM_FIT_GRID", unset = "0"), "1")

cov <- as.data.frame(dat$covariates_cen)
nl <- cbind(ageing = cov$ageing, unemployed = cov$unemployed)

message("=== ST-PCLM spatial ATA nested (mun→CT) ===")
fit_ata <- clgam_stpclm_ata_fit(
  y_mun = dat$ym,
  e_ct = dat$ec,
  C_mun_ct = dat$C_m,
  lon_ct = dat$xxc[, 1],
  lat_ct = dat$xxc[, 2],
  nlcovfine = nl,
  ndx = ndx,
  ndxnl = c(ndxnl, ndxnl),
  trace = TRUE
)

message(sprintf(
  "ATA  AIC=%.3f  BIC=%.3f  ed=%.3f  mass|C mu - y|=%.4g",
  fit_ata$aic, fit_ata$bic, fit_ata$ed,
  max(abs(as.numeric(as.matrix(dat$C_m) %*% (dat$ec * exp(fit_ata$eta))) - dat$ym))
))

message("=== Build nested raster compositions (res_m=", res_m, ") ===")
comp <- clgam_stpclm_ata_raster_C(dat, res_m = res_m)
message(sprintf(
  "cells=%d  C_mun_grid %d x %d",
  ncol(comp$C_mun_grid),
  nrow(comp$C_mun_grid), ncol(comp$C_mun_grid)
))
# sanity: each mun's cells count vs CT counts
message(sprintf(
  "rowSums C_mun_ct in [%s]; rowSums C_mun_grid in [%s]",
  paste(range(rowSums(as.matrix(dat$C_m))), collapse = ","),
  paste(range(rowSums(comp$C_mun_grid)), collapse = ",")
))

# Example C_t (identity in time) → C_st for documentation / future ST
C_t_id <- diag(1)
C_st_demo <- clgam_C_st(as.matrix(dat$C_m), C_t_id)
stopifnot(identical(dim(C_st_demo), dim(as.matrix(dat$C_m))))

out_grid <- NULL
if (fit_grid) {
  message("=== ST-PCLM spatial ATP via nested grid (C_mun_grid) ===")
  # covariates constant on cells within CT
  nl_cell <- nl[comp$ct_id, , drop = FALSE]
  ndx_g <- if (CLGAM_FAST) c(8L, 8L) else ndx
  fit_g <- pois_SAP(
    y = dat$ym,
    x1 = comp$xy_grid[, 1],
    x2 = comp$xy_grid[, 2],
    efine = comp$e_cell,
    nlcovfine = nl_cell,
    C = comp$C_mun_grid,
    ndx = ndx_g,
    ndxnl = c(ndxnl, ndxnl),
    elements = TRUE,
    trace = TRUE
  )
  # aggregate grid eta rates to CT for comparison: mu_cell = e_cell*exp(eta)
  mu_cell <- comp$e_cell * exp(fit_g$eta)
  mu_ct <- clgam_grid_to_ct(mu_cell, comp$C_ct_grid)
  eta_ct <- log(pmax(mu_ct, .Machine$double.xmin) / pmax(dat$ec, .Machine$double.xmin))
  out_grid <- list(
    aic = fit_g$aic,
    bic = fit_g$bic,
    ed = fit_g$ed,
    eta_grid = fit_g$eta,
    eta_ct_agg = eta_ct,
    mse_log_vs_ata = mean((eta_ct - fit_ata$eta)^2)
  )
  message(sprintf(
    "GRID AIC=%.3f  MSE(eta_ct_agg vs ATA eta)=%.5f",
    fit_g$aic, out_grid$mse_log_vs_ata
  ))
}

saveRDS(
  list(
    note = comp$note,
    st_pclm_ref = "Lee et al. PLoS ONE 2022 e0263711",
    adaptation = c(
      "C_s = nested mun→CT (ATA) instead of mun→grid indicator only",
      "Optional nested raster: C_mun_grid = C_mun_ct %*% C_ct_grid",
      "No C_t in Madrid CVD run (spatial section of ST-PCLM only)",
      "Estimator pois_SAP (SOP) — same family as CL-GAMM / ST-PCLM GLMM"
    ),
    fast = CLGAM_FAST,
    res_m = res_m,
    fit_ata = fit_ata[intersect(
      names(fit_ata),
      c("aic", "bic", "ed", "edf", "eta", "tau2", "composition", "st_pclm_note")
    )],
    composition = list(
      dim_C_mun_ct = dim(comp$C_mun_ct),
      dim_C_ct_grid = dim(comp$C_ct_grid),
      dim_C_mun_grid = dim(comp$C_mun_grid),
      n_cell = ncol(comp$C_mun_grid),
      n_empty_painted = comp$rasters$n_empty_painted
    ),
    fit_grid = out_grid
  ),
  file.path(CLGAM_OUTPUT, "st_pclm_ata_nested.rds")
)

message("Wrote ", file.path(CLGAM_OUTPUT, "st_pclm_ata_nested.rds"))
message("Set CLGAM_FIT_GRID=1 to also fit mun→grid PCLM (heavier).")

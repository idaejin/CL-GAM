#!/usr/bin/env Rscript
# LMMsolver baseline for municipal spatial smooth (C = I only).
# NOT a drop-in for composite-link CLMM (Case A mun→ct): LMMsolve has no C.
# Docs: https://biometris.github.io/LMMsolver/
#
# Coordinates are standardized (LMMsolver Cholesky is sensitive to UTM scale).

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))

.libPaths(c(file.path(CLGAM_EXP, "R_libs"), .libPaths()))
suppressPackageStartupMessages(library(LMMsolver))

dat <- clgam_load_madrid(FALSE)
lon_s <- as.numeric(scale(dat$xxm[, 1]))
lat_s <- as.numeric(scale(dat$xxm[, 2]))
df <- data.frame(
  y = dat$ym,
  lon = lon_s,
  lat = lat_s,
  loge = log(pmax(dat$em, 1e-8))
)

message("=== pois_SAP model0 (SOP CLMM, C = I, raw coords) ===")
t0 <- system.time({
  m0 <- pois_SAP(
    y = dat$ym, x1 = dat$xxm[, 1], x2 = dat$xxm[, 2],
    efine = dat$em, C = diag(length(dat$ym)),
    ndx = c(20, 20), elements = TRUE
  )
})
message(sprintf("pois_SAP: AIC=%.3f time=%.2fs", m0$aic, t0["elapsed"]))

message("=== LMMsolver Poisson + spl2D (REML, scaled coords, C = I) ===")
t1 <- system.time({
  lmm <- LMMsolve(
    fixed = y ~ 1,
    spline = ~spl2D(x1 = lon, x2 = lat, nseg = c(20, 20)),
    family = poisson(),
    offset = "loge",
    data = df,
    trace = FALSE
  )
})
pr_eta <- as.numeric(lmm$yhat) # linear predictor = log relative risk (offset separate)

message(sprintf("LMMsolver time=%.2fs", t1["elapsed"]))
message(sprintf(
  "corr(eta_SAP, eta_LMM)=%.6f  RMSE=%.3e",
  cor(m0$eta, pr_eta),
  sqrt(mean((m0$eta - pr_eta)^2))
))
message(sprintf("speedup vs pois_SAP (C=I): %.1fx", t0["elapsed"] / t1["elapsed"]))

saveRDS(
  list(
    pois_SAP_sec = unname(t0["elapsed"]),
    LMMsolver_sec = unname(t1["elapsed"]),
    pois_SAP_aic = m0$aic,
    corr = cor(m0$eta, pr_eta),
    rmse = sqrt(mean((m0$eta - pr_eta)^2)),
    note = "Use lmm$yhat (not predict$ypred). C=I only; scale coords for LMMsolver."
  ),
  file.path(CLGAM_OUTPUT, "lmmsolver_mun_baseline.rds")
)
message("Wrote ", file.path(CLGAM_OUTPUT, "lmmsolver_mun_baseline.rds"))

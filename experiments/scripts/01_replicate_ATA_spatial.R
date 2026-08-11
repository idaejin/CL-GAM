#!/usr/bin/env Rscript
# Replicate PAPER Codes2 models 0–1: municipal GAM vs ATA CLMM (no covariates).
# Full run: minutes. Fast: Sys.setenv(CLGAM_FAST = "1")
#
# From experiments/:  Rscript scripts/01_replicate_ATA_spatial.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))

dat <- clgam_load_madrid(load_maps = TRUE)
ndx <- CLGAM_NDX_SPATIAL

message("=== model0: Poisson GLMM at municipality (C = I) ===")
model0 <- pois_SAP(
  y = dat$ym,
  x1 = dat$xxm[, 1],
  x2 = dat$xxm[, 2],
  efine = dat$em,
  C = diag(length(dat$ym)),
  ndx = ndx,
  elements = TRUE
)
message(sprintf("AIC=%.3f  BIC=%.3f  ed=%.3f", model0$aic, model0$bic, model0$ed))

message("=== model1: Poisson CLMM mun→ct (true e_fine) ===")
model1 <- pois_SAP(
  y = dat$ym,
  x1 = dat$xxc[, 1],
  x2 = dat$xxc[, 2],
  efine = dat$ec,
  C = dat$C_m,
  ndx = ndx,
  elements = TRUE
)
message(sprintf("AIC=%.3f  BIC=%.3f  ed=%.3f", model1$aic, model1$bic, model1$ed))

# Paper reported (ndx=20,20): model0 AIC ~374.6; model1 AIC ~396.9
brks <- quantile(dat$logSMR_mun, seq(0, 1, 1 / 10))
brks[1] <- brks[1] - 0.001
brks[length(brks)] <- brks[length(brks)] + 0.001

pdf(file.path(CLGAM_OUTPUT, "ATA_raw_logSMR_mun.pdf"))
choromap(
  sf = dat$map_mun$sf, values = dat$logSMR_mun, breaks = brks,
  border = "grey30", lwd = 0.1, title = "Raw log SMR (mun)"
)
dev.off()

pdf(file.path(CLGAM_OUTPUT, "ATA_model0_eta_mun.pdf"))
choromap(
  sf = dat$map_mun$sf, values = model0$eta, breaks = brks,
  border = "grey30", lwd = 0.1, title = "model0 eta (mun)"
)
dev.off()

pdf(file.path(CLGAM_OUTPUT, "ATA_model1_eta_ct.pdf"))
choromap(
  sf = dat$map_ct$sf, values = model1$eta, breaks = brks,
  border = "grey30", lwd = 0.1, title = "model1 eta (ct)"
)
dev.off()

saveRDS(
  list(
    ndx = ndx, fast = CLGAM_FAST,
    model0 = model0[c("aic", "bic", "ed", "edf", "dev", "eta", "tau2")],
    model1 = model1[c("aic", "bic", "ed", "edf", "dev", "eta", "tau2")],
    paper_ref = list(model0_aic = 374.5841, model1_aic = 396.9185)
  ),
  file.path(CLGAM_OUTPUT, "ATA_spatial_fit.rds")
)

message("Wrote figures + RDS under ", CLGAM_OUTPUT)

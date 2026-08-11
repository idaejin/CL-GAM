#!/usr/bin/env Rscript
# Replicate PAPER Codes2 model 2.5 (Case A): spatial + ageing + unemployment at CT.
# Columns: expl[,5]=ageing, expl[,2]=unemployed (centered).
#
# From experiments/:  Rscript scripts/02_replicate_CaseA_covariates.R
# Fast: Sys.setenv(CLGAM_FAST = "1")

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))

dat <- clgam_load_madrid(load_maps = FALSE)
ndx <- CLGAM_NDX_SPATIAL
ndxnl <- CLGAM_NDX_NL

ageing <- dat$covariates_cen[, "ageing", drop = FALSE]
unemp <- dat$covariates_cen[, "unemployed", drop = FALSE]

message("=== model2.5: nlcovfine = ageing + unemployment ===")
model25 <- pois_SAP(
  y = dat$ym,
  x1 = dat$xxc[, 1],
  x2 = dat$xxc[, 2],
  efine = dat$ec,
  nlcovfine = cbind(ageing, unemp),
  C = dat$C_m,
  ndx = ndx,
  ndxnl = c(ndxnl, ndxnl),
  elements = TRUE,
  trace = TRUE
)

message(sprintf(
  "AIC=%.3f  BIC=%.3f  ed=%.3f  edf=%s",
  model25$aic, model25$bic, model25$ed,
  paste(round(model25$edf, 3), collapse = ", ")
))

saveRDS(
  list(
    ndx = ndx, ndxnl = ndxnl, fast = CLGAM_FAST,
    fit = model25[intersect(
      names(model25),
      c("aic", "bic", "ed", "edf", "dev", "eta", "tau2", "beta", "theta")
    )]
  ),
  file.path(CLGAM_OUTPUT, "CaseA_model25_fit.rds")
)

message("Wrote ", file.path(CLGAM_OUTPUT, "CaseA_model25_fit.rds"))

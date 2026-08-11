#!/usr/bin/env Rscript
# Malone Method 1 competitor — defaults match CL-GAMM manuscript / Codes3 science.
# Bugfix: order-safe mass-balance. Adaptations documented in
# improvements/MALONE_ADAPTATIONS.md and in the fit$adaptations vector.
#
# Optional (report if set):
#   CLGAM_MALONE_SPATIAL = none|additive|te   (default none = Method 1)
#   CLGAM_MALONE_FAMILY  = poisson|quasipoisson (default poisson)
#
# From experiments/:  Rscript scripts/05_malone_corrected.R
# Fast: Sys.setenv(CLGAM_FAST = "1")

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/02_load_madrid.R"))
source(file.path(root, "R/03_malone.R"))

dat <- clgam_load_madrid(load_maps = FALSE)
cov_use <- as.data.frame(dat$covariates_cen)[, c("ageing", "unemployed"), drop = FALSE]

spatial <- Sys.getenv("CLGAM_MALONE_SPATIAL", unset = "none")
family <- Sys.getenv("CLGAM_MALONE_FAMILY", unset = "poisson")
k_sp <- if (CLGAM_FAST) 40L else 80L
k_cov <- if (CLGAM_FAST) 6L else 10L
max_iter <- if (CLGAM_FAST) 50L else 300L

message("=== Malone Method 1: spatial=", spatial, " family=", family, " ===")
fit <- clgam_malone_fit(
  y_mun = dat$ym,
  e_fine = dat$ec,
  C = dat$C_m,
  cov_fine = cov_use,
  lon = dat$xxc[, 1],
  lat = dat$xxc[, 2],
  cov_smooth = c("ageing", "unemployed"),
  spatial = spatial,
  family = family,
  round_counts = identical(family, "poisson"),
  max_iter = max_iter,
  tol = 1e-3,
  k_spatial = k_sp,
  k_cov = k_cov,
  trace = TRUE
)

mse_log <- mean((log((dat$yc + 1) / dat$ec) - fit$eta)^2)
mse_y <- mean((dat$yc - fit$fitted)^2)

message(sprintf(
  "converged=%s  n_iter=%d  crit_map=%.5g  mass_fitted=%.4g  mass_balanced=%.4g",
  fit$converged, fit$n_iter, fit$crit_map, fit$mass_error_fitted, fit$mass_error_balanced
))
message(sprintf("MSE logSMR vs CT obs: %.6f", mse_log))
message(sprintf("MSE counts vs CT obs:  %.4f", mse_y))
message("Adaptations vs Malone 2012:")
for (a in fit$adaptations) message("  - ", a)
message("Bugfix: ", fit$bugfixes)

tr <- predict(fit$model, type = "terms")
cn <- colnames(tr)
ix_age <- grep("ageing", cn)
ix_unemp <- grep("unemployed", cn)

partial <- data.frame(
  ageing_pct = dat$covariates[, "ageing"] * 100,
  unemployed_pct = dat$covariates[, "unemployed"] * 100,
  s_ageing = if (length(ix_age)) as.numeric(scale(tr[, ix_age[1]], scale = FALSE)) else NA_real_,
  s_unemployed = if (length(ix_unemp)) as.numeric(scale(tr[, ix_unemp[1]], scale = FALSE)) else NA_real_
)

out <- list(
  method = "Malone Method 1 (CL-GAMM draft)",
  reference = fit$reference,
  adaptations = fit$adaptations,
  bugfixes = fit$bugfixes,
  fast = CLGAM_FAST,
  spatial = fit$spatial,
  family = fit$family,
  round_counts = fit$round_counts,
  converged = fit$converged,
  n_iter = fit$n_iter,
  crit_map = fit$crit_map,
  crit_balance = fit$crit_balance,
  crit_fitgap = fit$crit_fitgap,
  mass_error_fitted = fit$mass_error_fitted,
  mass_error_balanced = fit$mass_error_balanced,
  mse_logSMR = mse_log,
  mse_counts = mse_y,
  formula = paste(deparse(fit$formula), collapse = " "),
  eta = fit$eta,
  fitted = fit$fitted,
  edf = sum(fit$model$edf),
  aic = tryCatch(AIC(fit$model), error = function(e) NA_real_),
  partial = partial,
  history = fit$history
)

saveRDS(out, file.path(CLGAM_OUTPUT, "malone_corrected_fit.rds"))
write.csv(
  data.frame(
    spatial = fit$spatial,
    family = fit$family,
    converged = fit$converged,
    n_iter = fit$n_iter,
    crit_map = fit$crit_map,
    mass_error_fitted = fit$mass_error_fitted,
    mass_error_balanced = fit$mass_error_balanced,
    mse_logSMR = mse_log,
    mse_counts = mse_y,
    aic = out$aic,
    edf = out$edf
  ),
  file.path(CLGAM_OUTPUT, "malone_corrected_summary.csv"),
  row.names = FALSE
)
if (!is.null(fit$history)) {
  write.csv(fit$history, file.path(CLGAM_OUTPUT, "malone_corrected_history.csv"), row.names = FALSE)
}

message("Wrote ", file.path(CLGAM_OUTPUT, "malone_corrected_fit.rds"))
message("See improvements/MALONE_ADAPTATIONS.md")

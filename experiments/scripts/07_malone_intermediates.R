#!/usr/bin/env Rscript
# Compare Method 1 Malone vs revive *intermediate* improvements (still heuristic;
# C not in the likelihood). All REVIVE flags are labeled in fit$adaptations.
#
# From experiments/:  Rscript scripts/07_malone_intermediates.R
# Fast: Sys.setenv(CLGAM_FAST = "1")

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/02_load_madrid.R"))
source(file.path(root, "R/03_malone.R"))

dat <- clgam_load_madrid(load_maps = FALSE)
cov_use <- as.data.frame(dat$covariates_cen)[, c("ageing", "unemployed"), drop = FALSE]

k_cov <- if (CLGAM_FAST) 6L else 10L
max_iter <- if (CLGAM_FAST) 40L else 200L
tol <- 1e-3

variants <- list(
  list(
    id = "M1_method1",
    label = "Method 1 (poisson+round, equal init, deliver fitted)",
    family = "poisson", round_counts = TRUE, init = "equal",
    deliver = "fitted", stop_on = "map"
  ),
  list(
    id = "I1_deliver_balanced",
    label = "REVIVE: Method1 + deliver mass-balanced map",
    family = "poisson", round_counts = TRUE, init = "equal",
    deliver = "balanced", stop_on = "map"
  ),
  list(
    id = "I2_quasi_noround",
    label = "REVIVE: quasipoisson, no round, deliver fitted",
    family = "quasipoisson", round_counts = FALSE, init = "equal",
    deliver = "fitted", stop_on = "map"
  ),
  list(
    id = "I3_quasi_balanced",
    label = "REVIVE: quasipoisson + deliver balanced",
    family = "quasipoisson", round_counts = FALSE, init = "equal",
    deliver = "balanced", stop_on = "balanced_map"
  ),
  list(
    id = "I4_exposure_init",
    label = "REVIVE: exposure init (mun SMR x e), poisson+round",
    family = "poisson", round_counts = TRUE, init = "exposure",
    deliver = "fitted", stop_on = "map"
  ),
  list(
    id = "I5_exposure_quasi_bal",
    label = "REVIVE: exposure init + quasi + deliver balanced",
    family = "quasipoisson", round_counts = FALSE, init = "exposure",
    deliver = "balanced", stop_on = "balanced_map"
  )
)

rows <- list()
fits <- list()

for (v in variants) {
  message("\n======== ", v$id, " ========")
  fit <- clgam_malone_fit(
    y_mun = dat$ym,
    e_fine = dat$ec,
    C = dat$C_m,
    cov_fine = cov_use,
    lon = dat$xxc[, 1],
    lat = dat$xxc[, 2],
    cov_smooth = c("ageing", "unemployed"),
    spatial = "none",
    family = v$family,
    round_counts = v$round_counts,
    init = v$init,
    deliver = v$deliver,
    stop_on = v$stop_on,
    max_iter = max_iter,
    tol = tol,
    k_cov = k_cov,
    trace = TRUE
  )
  mse_log_gam <- mean((log((dat$yc + 1) / dat$ec) - fit$eta_gam)^2)
  mse_log_del <- mean((log((dat$yc + 1) / dat$ec) - fit$eta)^2)
  mse_y_del <- mean((dat$yc - fit$fitted)^2)
  mse_y_gam <- mean((dat$yc - fit$fitted_gam)^2)
  rows[[v$id]] <- data.frame(
    id = v$id,
    label = v$label,
    family = fit$family,
    round_counts = fit$round_counts,
    init = fit$init,
    deliver = fit$deliver,
    stop_on = fit$stop_on,
    converged = fit$converged,
    n_iter = fit$n_iter,
    crit_map = fit$crit_map,
    crit_bal_map = fit$crit_bal_map,
    mass_delivered = fit$mass_error_delivered,
    mass_fitted_gam = fit$mass_error_fitted,
    mse_logSMR_gam = mse_log_gam,
    mse_logSMR_delivered = mse_log_del,
    mse_counts_gam = mse_y_gam,
    mse_counts_delivered = mse_y_del,
    stringsAsFactors = FALSE
  )
  fits[[v$id]] <- fit[c(
    "eta", "eta_gam", "fitted", "fitted_gam", "n_iter", "converged", "adaptations",
    "mass_error_delivered", "deliver", "init", "family"
  )]
  message(sprintf(
    "→ MSE log(gam)=%.5f  MSE log(del)=%.5f  MSE y(gam)=%.3f  MSE y(del)=%.3f  mass_del=%.3g  iters=%d",
    mse_log_gam, mse_log_del, mse_y_gam, mse_y_del, fit$mass_error_delivered, fit$n_iter
  ))
}

cmp <- do.call(rbind, rows)
cmp <- cmp[order(cmp$mse_logSMR_gam), ]
print(cmp[, c(
  "id", "mse_logSMR_gam", "mse_logSMR_delivered",
  "mse_counts_gam", "mse_counts_delivered", "mass_delivered", "n_iter"
)])

saveRDS(
  list(
    note = paste(
      "Intermediate Malone improvements (heuristic only; C not in likelihood).",
      "REVIVE variants labeled in adaptations. Method 1 = M1_method1."
    ),
    fast = CLGAM_FAST,
    comparison = cmp,
    fits = fits
  ),
  file.path(CLGAM_OUTPUT, "malone_intermediates.rds")
)
write.csv(cmp, file.path(CLGAM_OUTPUT, "malone_intermediates.csv"), row.names = FALSE)

message("\nWrote ", file.path(CLGAM_OUTPUT, "malone_intermediates.csv"))
message("Best by MSE logSMR (GAM eta): ", cmp$id[1])
message("Best by MSE counts (delivered): ", cmp$id[order(cmp$mse_counts_delivered)[1]])

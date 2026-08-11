#!/usr/bin/env Rscript
# CL-GAMM Case A/B/C on SpatialEpi::pennLC (67 counties as fine units).
#
# Coarse: 20 k-means regions of county centroids (more y than NY8's 8 counties).
#   A: s(smoking) at county
#   B: linear north_agg (region-mean latitude expanded)
#   C: A + B
#
# From experiments/:  Rscript scripts/17_penn_case_ABC.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/07_load_pennLC.R"))

dat <- clgam_load_pennLC(n_coarse = 20L, seed = 1L)
ndx <- if (CLGAM_FAST) c(6L, 6L) else c(10L, 10L)
ndxnl <- if (CLGAM_FAST) 5L else 8L

message(
  "pennLC nested ATA: ", dat$n_coarse, " coarse -> ", dat$n_fine,
  " counties | scheme=", dat$coarse_scheme,
  " | sum(y)=", sum(dat$y),
  " | max|C e - e_g|=", dat$e_mass_check
)
message(dat$note)

fit_summary <- function(fit, label) {
  list(
    case = label, aic = fit$aic, bic = fit$bic, ed = fit$ed,
    edf = fit$edf, niter = fit$niter, var.comp = fit$var.comp,
    eta = fit$eta, elapsed = fit$elapsed.time,
    nleffects = fit$nleffects, sdnleffects = fit$sdnleffects,
    leffects = fit$leffects, sdleffects = fit$sdleffects
  )
}

message("=== Case A: s(smoking) ===")
fit_a <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_fine,
  nlcovfine = dat$covariates_fine, C = dat$C,
  ndx = ndx, ndxnl = ndxnl, elements = TRUE, trace = TRUE
)

message("=== Case B: linear north_agg ===")
fit_b <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_fine,
  lcovfine = dat$covariates_agg_expanded, C = dat$C,
  ndx = ndx, elements = TRUE, trace = TRUE
)

message("=== Case C: s(smoking) + north_agg ===")
fit_c <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_fine,
  lcovfine = dat$covariates_agg_expanded,
  nlcovfine = dat$covariates_fine, C = dat$C,
  ndx = ndx, ndxnl = ndxnl, elements = TRUE, trace = TRUE
)

log_smr <- log((dat$y_fine + 0.5) / dat$e_fine)
mse <- function(eta) mean((eta - log_smr)^2)
mass <- function(eta) {
  max(abs(as.numeric(dat$C %*% (dat$e_fine * exp(eta))) - dat$y))
}

tab <- data.frame(
  case = c("A", "B", "C"),
  aic = c(fit_a$aic, fit_b$aic, fit_c$aic),
  bic = c(fit_a$bic, fit_b$bic, fit_c$bic),
  ed = c(fit_a$ed, fit_b$ed, fit_c$ed),
  mse_vs_county = c(mse(fit_a$eta), mse(fit_b$eta), mse(fit_c$eta)),
  cor_vs_county = c(
    cor(fit_a$eta, log_smr), cor(fit_b$eta, log_smr), cor(fit_c$eta, log_smr)
  ),
  mass_err = c(mass(fit_a$eta), mass(fit_b$eta), mass(fit_c$eta)),
  elapsed_s = c(fit_a$elapsed.time, fit_b$elapsed.time, fit_c$elapsed.time),
  stringsAsFactors = FALSE
)

message("=== Summary ===")
print(tab)

out <- list(
  mode = "pennLC_kmeans20",
  data = list(
    n_coarse = dat$n_coarse, n_fine = dat$n_fine,
    coarse_scheme = dat$coarse_scheme, y = dat$y, C = dat$C,
    county = dat$county, group = dat$group,
    y_fine = dat$y_fine, e_fine = dat$e_fine,
    smoking = dat$smoking, x1 = dat$x1, x2 = dat$x2,
    source = dat$source, note = dat$note
  ),
  ndx = ndx, ndxnl = ndxnl, fast = CLGAM_FAST,
  case_a = fit_summary(fit_a, "A"),
  case_b = fit_summary(fit_b, "B"),
  case_c = fit_summary(fit_c, "C"),
  table = tab
)

saveRDS(out, file.path(CLGAM_OUTPUT, "penn_case_ABC_fit.rds"))
write.csv(tab, file.path(CLGAM_OUTPUT, "penn_case_ABC_summary.csv"), row.names = FALSE)
message("Wrote ", file.path(CLGAM_OUTPUT, "penn_case_ABC_fit.rds"))
message("Wrote ", file.path(CLGAM_OUTPUT, "penn_case_ABC_summary.csv"))

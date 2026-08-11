#!/usr/bin/env Rscript
# CL-GAMM Case A / B / C on NY leukemia (nested ATA).
#
# Data: SpatialEpi::NYleukemia + DClusterm::NY8
#
# Default (true counties, poor bases — numerically stable with n=8):
#   C: 8 counties × 281 tracts
#   ndx = (3,3); linear covariates only
#   A: ownhome (tract)
#   B: tce_agg = county-mean PEXPOSURE expanded
#   C: ownhome + tce_agg
#
# Optional richer mid-level (k-means 40) for nonlinear fine smooths:
#   Sys.setenv(CLGAM_NY_COARSE = "kmeans")
#
# From experiments/:  Rscript scripts/15_ny_case_ABC.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/06_load_ny_leukemia.R"))

use_kmeans <- identical(Sys.getenv("CLGAM_NY_COARSE", unset = "county"), "kmeans")

if (use_kmeans) {
  dat <- clgam_load_ny_leukemia(n_coarse = 40L, seed = 1L)
  ndx <- if (CLGAM_FAST) c(6L, 6L) else c(10L, 10L)
  ndxnl <- if (CLGAM_FAST) 4L else 8L
  mode <- "kmeans40_nl"
  message("Mode: k-means mid-level with nonlinear fine smooths")
} else {
  dat <- clgam_load_ny_leukemia(n_coarse = "county")
  ndx <- c(3L, 3L)
  ndxnl <- NA_integer_
  mode <- "county8_linear_poor_bases"
  message("Mode: true 8 counties + poor spatial bases + linear A/B/C")
}

message(
  "Coarse scheme: ", dat$coarse_scheme,
  "  ndx=(", paste(ndx, collapse = ","), ")",
  if (use_kmeans) paste0(" ndxnl=", ndxnl) else " (linear covariates)"
)
message(
  "NY nested ATA: ", dat$n_coarse, " coarse → ", dat$n_tract,
  " tracts  |  sum(y)=", sum(dat$y)
)
message("Sources: ", paste(dat$source, collapse = " + "))

fit_summary <- function(fit, label) {
  list(
    case = label,
    aic = fit$aic,
    bic = fit$bic,
    ed = fit$ed,
    edf = fit$edf,
    niter = fit$niter,
    var.comp = fit$var.comp,
    eta = fit$eta,
    elapsed = fit$elapsed.time
  )
}

own <- dat$covariates_fine[, "ownhome", drop = FALSE]
age <- dat$covariates_fine[, "ageing", drop = FALSE]
tce <- dat$covariates_agg_expanded

if (use_kmeans) {
  message("=== Case A: fine s(ageing)+s(ownhome) ===")
  fit_a <- pois_SAP(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_tract,
    nlcovfine = dat$covariates_fine, C = dat$C,
    ndx = ndx, ndxnl = c(ndxnl, ndxnl), elements = TRUE, trace = TRUE
  )
  message("=== Case B: aggregated linear tce_agg ===")
  fit_b <- pois_SAP(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_tract,
    lcovfine = tce, C = dat$C, ndx = ndx, elements = TRUE, trace = TRUE
  )
  message("=== Case C: fine smooths + tce_agg ===")
  fit_c <- pois_SAP(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_tract,
    lcovfine = tce, nlcovfine = dat$covariates_fine, C = dat$C,
    ndx = ndx, ndxnl = c(ndxnl, ndxnl), elements = TRUE, trace = TRUE
  )
  cov_note <- "A: s(ageing)+s(ownhome); B: linear tce_agg; C: A+B"
} else {
  # Poor bases + linear only (nonlinear / 2 fine linears often singular with n=8)
  message("=== Case A: linear ownhome (tract) ===")
  fit_a <- pois_SAP(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_tract,
    lcovfine = own, C = dat$C, ndx = ndx, elements = TRUE, trace = TRUE
  )
  message("=== Case B: linear tce_county (expanded) ===")
  fit_b <- pois_SAP(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_tract,
    lcovfine = tce, C = dat$C, ndx = ndx, elements = TRUE, trace = TRUE
  )
  message("=== Case C: ownhome + tce_county ===")
  fit_c <- pois_SAP(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_tract,
    lcovfine = cbind(own, tce), C = dat$C, ndx = ndx, elements = TRUE, trace = TRUE
  )
  cov_note <- paste(
    "A: linear PCTOWNHOME; B: linear county-mean PEXPOSURE expanded;",
    "C: A+B. ageing alone / 2 fine linears / P-splines often fail with n=8."
  )
}

log_smr_tr <- log((dat$cases_tract + 0.5) / dat$e_tract)
mse <- function(eta) mean((eta - log_smr_tr)^2)

tab <- data.frame(
  case = c("A", "B", "C"),
  aic = c(fit_a$aic, fit_b$aic, fit_c$aic),
  bic = c(fit_a$bic, fit_b$bic, fit_c$bic),
  ed = c(fit_a$ed, fit_b$ed, fit_c$ed),
  mse_vs_tract = c(mse(fit_a$eta), mse(fit_b$eta), mse(fit_c$eta)),
  mass_err = c(
    max(abs(as.numeric(dat$C %*% (dat$e_tract * exp(fit_a$eta))) - dat$y)),
    max(abs(as.numeric(dat$C %*% (dat$e_tract * exp(fit_b$eta))) - dat$y)),
    max(abs(as.numeric(dat$C %*% (dat$e_tract * exp(fit_c$eta))) - dat$y))
  ),
  elapsed_s = c(fit_a$elapsed.time, fit_b$elapsed.time, fit_c$elapsed.time),
  stringsAsFactors = FALSE
)

message("=== Summary ===")
print(tab)
message("Covariates: ", cov_note)

out <- list(
  mode = mode,
  cov_note = cov_note,
  data = list(
    n_coarse = dat$n_coarse,
    n_tract = dat$n_tract,
    coarse_scheme = dat$coarse_scheme,
    y = dat$y,
    C = dat$C,
    group_levels = dat$group_levels,
    covariates_agg = dat$covariates_agg,
    note = dat$note,
    source = dat$source
  ),
  ndx = ndx,
  ndxnl = ndxnl,
  fast = CLGAM_FAST,
  case_a = fit_summary(fit_a, "A"),
  case_b = fit_summary(fit_b, "B"),
  case_c = fit_summary(fit_c, "C"),
  table = tab
)

saveRDS(out, file.path(CLGAM_OUTPUT, "ny_case_ABC_fit.rds"))
write.csv(tab, file.path(CLGAM_OUTPUT, "ny_case_ABC_summary.csv"), row.names = FALSE)
message("Wrote ", file.path(CLGAM_OUTPUT, "ny_case_ABC_fit.rds"))
message("Wrote ", file.path(CLGAM_OUTPUT, "ny_case_ABC_summary.csv"))

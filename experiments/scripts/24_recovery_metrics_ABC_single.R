#!/usr/bin/env Rscript
# Single-replica recovery metrics for Case A/B/C.
#
# Uses the exact geometry/seed/fit settings that generate:
#   experiments/output/clgam_recovery_caseA.png
#   experiments/output/clgam_recovery_caseB.png
#   experiments/output/clgam_recovery_caseC.png
#
# Outputs:
#   experiments/output/recovery_ABC_metrics_single.csv
#   experiments/output/recovery_ABC_metrics_single.rds

args_root <- if (dir.exists("clgam")) "." else if (dir.exists("../clgam")) ".." else {
  stop("Run from CL-GAM root or experiments/")
}
setwd(args_root)

.libPaths(c(normalizePath("experiments/R_libs"), .libPaths()))

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all("clgam", quiet = TRUE)
  } else {
    stop("Falta 'devtools' para cargar clgam; usa devtools::load_all().")
  }
  library(stats)
})

loglik_pois <- function(y, mu, eps = .Machine$double.xmin) {
  y <- as.numeric(y)
  mu <- pmax(as.numeric(mu), eps)
  sum(y * log(mu) - mu - lgamma(y + 1))
}

deviance_pois <- function(y, mu, eps = .Machine$double.xmin) {
  y <- as.numeric(y)
  mu <- pmax(as.numeric(mu), eps)
  sat <- ifelse(y > 0, y * log(y) - y, -y) # consistent with 0*log(0)=0
  fit <- y * log(mu) - mu
  2 * sum(sat - fit)
}

mse <- function(x, truth) mean((x - truth)^2)
rmse <- function(x, truth) sqrt(mse(x, truth))
mae <- function(x, truth) mean(abs(x - truth))

metric_row <- function(method, mse_eta, rmse_eta, mae_eta, mse_mu, rmse_mu, mae_mu,
                        mass_err_max, mass_err_l2, deviance_coarse, loglik_coarse,
                        mse_f = NA_real_, mse_nl = NA_real_, mse_g = NA_real_, mse_h = NA_real_) {
  data.frame(
    method = method,
    mse_eta = mse_eta,
    rmse_eta = rmse_eta,
    mae_eta = mae_eta,
    mse_mu = mse_mu,
    rmse_mu = rmse_mu,
    mae_mu = mae_mu,
    mass_err_max = mass_err_max,
    mass_err_l2 = mass_err_l2,
    deviance_coarse = deviance_coarse,
    loglik_coarse = loglik_coarse,
    mse_f = mse_f,
    mse_nl = mse_nl,
    mse_g = mse_g,
    mse_h = mse_h,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Case A/B (shared geometry/scales in scripts/21_recovery_caseAB.R)
# ---------------------------------------------------------------------------
nl_amp <- 1.2
nl_fun <- function(z) sin(2 * pi * z)
n_coarse <- 40L
n_fine_per <- 10L
spatial_amp_ab <- 0.8

simulate_fit_caseAB <- function(case = c("A", "B")) {
  case <- match.arg(case)
  if (identical(case, "A")) {
    dat <- simulate_ata(
      n_coarse = n_coarse, n_fine_per = n_fine_per, seed = 99L,
      include_covariate = TRUE, covariate_level = "fine",
      spatial_amp = spatial_amp_ab, nl_amp = nl_amp, nl_fun = nl_fun
    )
    fit <- clgam(
      dat$y, dat$x1, dat$x2, dat$C,
      exposure = dat$efine, smooth = dat$nlcovfine,
      knots = c(12L, 12L), knots_nl = 12L,
      elements = TRUE, trace = FALSE
    )
  } else {
    dat <- simulate_ata(
      n_coarse = n_coarse, n_fine_per = n_fine_per, seed = 6L,
      include_covariate = TRUE, covariate_level = "coarse",
      spatial_amp = spatial_amp_ab, nl_amp = nl_amp, nl_fun = nl_fun
    )
    fit <- clgam(
      dat$y, dat$x1, dat$x2, dat$C,
      exposure = dat$efine, smooth = dat$nlcovfine,
      knots = c(10L, 10L), knots_nl = 6L,
      elements = TRUE, trace = FALSE
    )
  }

  # Same centering convention as scripts/21_recovery_caseAB.R prep()
  nl_hat <- fit$nleffects[, 1] - mean(fit$nleffects[, 1])
  f_hat <- fit$eta - fit$nleffects[, 1]
  f_hat <- f_hat - mean(f_hat)

  eta_hat <- fit$eta - mean(fit$eta)
  eta_true <- dat$eta_true - mean(dat$eta_true)

  # Spatial truth:
  f_true <- dat$eta_spatial_true - mean(dat$eta_spatial_true)

  # NL truth depending on case:
  if (identical(case, "B")) {
    nl_true <- dat$h_true - mean(dat$h_true)
    # h_true lives on coarse-support values, but nl_hat is evaluated on fine grid.
    # We use the available direct definition in the script: compare nl_hat pointwise on the
    # fine grid proxy defined by dat$h_true (stored aligned to fine evaluations).
    mse_nl_val <- mse(nl_hat, nl_true)
  } else {
    nl_true <- dat$g_true - mean(dat$g_true)
    mse_nl_val <- mse(nl_hat, nl_true)
  }

  # Fine means for mu diagnostics.
  mu_hat <- fit$mu
  mu_true <- dat$efine * exp(dat$eta_true) # true fine means from oracle eta

  # Mass calibration on coarse counts.
  mu_c <- as.numeric(dat$C %*% mu_hat)
  mass_res <- mu_c - dat$y
  mass_err_max <- max(abs(mass_res))
  mass_err_l2 <- sqrt(sum(mass_res^2))
  dev_coarse <- deviance_pois(dat$y, mu_c)
  ll_coarse <- loglik_pois(dat$y, mu_c)

  data.frame(
    dat = dat, fit = fit,
    eta_hat = eta_hat, eta_true = eta_true,
    f_hat = f_hat, f_true = f_true,
    nl_hat = nl_hat, nl_true = nl_true,
    mu_hat = mu_hat, mu_true = mu_true,
    mu_c = mu_c,
    mass_err_max = mass_err_max, mass_err_l2 = mass_err_l2,
    deviance_coarse = dev_coarse, loglik_coarse = ll_coarse
  )
}

res_A <- simulate_fit_caseAB("A")
res_B <- simulate_fit_caseAB("B")

row_A <- metric_row(
  method = "Case A (fine covariate)",
  mse_eta = mse(res_A$eta_hat, res_A$eta_true),
  rmse_eta = rmse(res_A$eta_hat, res_A$eta_true),
  mae_eta = mae(res_A$eta_hat, res_A$eta_true),
  mse_mu = mse(res_A$mu_hat, res_A$mu_true),
  rmse_mu = rmse(res_A$mu_hat, res_A$mu_true),
  mae_mu = mae(res_A$mu_hat, res_A$mu_true),
  mass_err_max = res_A$mass_err_max,
  mass_err_l2 = res_A$mass_err_l2,
  deviance_coarse = res_A$deviance_coarse,
  loglik_coarse = res_A$loglik_coarse,
  mse_f = mse(res_A$f_hat, res_A$f_true),
  mse_nl = mse(res_A$nl_hat, res_A$nl_true)
)

row_B <- metric_row(
  method = "Case B (coarse covariate)",
  mse_eta = mse(res_B$eta_hat, res_B$eta_true),
  rmse_eta = rmse(res_B$eta_hat, res_B$eta_true),
  mae_eta = mae(res_B$eta_hat, res_B$eta_true),
  mse_mu = mse(res_B$mu_hat, res_B$mu_true),
  rmse_mu = rmse(res_B$mu_hat, res_B$mu_true),
  mae_mu = mae(res_B$mu_hat, res_B$mu_true),
  mass_err_max = res_B$mass_err_max,
  mass_err_l2 = res_B$mass_err_l2,
  deviance_coarse = res_B$deviance_coarse,
  loglik_coarse = res_B$loglik_coarse,
  mse_f = mse(res_B$f_hat, res_B$f_true),
  mse_nl = mse(res_B$nl_hat, res_B$nl_true)
)

# ---------------------------------------------------------------------------
# Case C (shared geometry/scales in scripts/23_recovery_ABC.R)
# ---------------------------------------------------------------------------
nl_funs <- list(
  function(z) -tanh(2.5 * (z - 0.5)),
  function(z) tanh(2.5 * (z - 0.5))
)
SEED <- 6L
spatial_amp_c <- 0.75
nl_amp_c <- c(1.2, 1.2)

dat_C <- simulate_ata(
  n_coarse = n_coarse, n_fine_per = n_fine_per, seed = SEED,
  include_covariate = TRUE, covariate_level = "both",
  spatial_amp = spatial_amp_c, nl_amp = nl_amp_c,
  nl_fun = nl_funs
)

fit_C <- clgam(
  dat_C$y, dat_C$x1, dat_C$x2, dat_C$C,
  exposure = dat_C$efine, smooth = dat_C$nlcovfine,
  knots = c(10L, 10L), knots_nl = c(8L, 6L),
  elements = TRUE, trace = FALSE
)

g_hat <- fit_C$nleffects[, 1] - mean(fit_C$nleffects[, 1])
h_hat <- fit_C$nleffects[, 2] - mean(fit_C$nleffects[, 2])
f_hat <- fit_C$eta - fit_C$nleffects[, 1] - fit_C$nleffects[, 2]
f_hat <- f_hat - mean(f_hat)
eta_hat <- fit_C$eta - mean(fit_C$eta)
eta_true <- dat_C$eta_true - mean(dat_C$eta_true)

g_true <- dat_C$g_true - mean(dat_C$g_true)
h_true <- dat_C$h_true - mean(dat_C$h_true)
f_true <- dat_C$eta_spatial_true - mean(dat_C$eta_spatial_true)

mu_hat <- fit_C$mu
mu_true <- dat_C$efine * exp(dat_C$eta_true)

mu_c <- as.numeric(dat_C$C %*% mu_hat)
mass_res <- mu_c - dat_C$y
mass_err_max <- max(abs(mass_res))
mass_err_l2 <- sqrt(sum(mass_res^2))
dev_coarse <- deviance_pois(dat_C$y, mu_c)
ll_coarse <- loglik_pois(dat_C$y, mu_c)

row_C <- metric_row(
  method = "Case C (two smooths)",
  mse_eta = mse(eta_hat, eta_true),
  rmse_eta = rmse(eta_hat, eta_true),
  mae_eta = mae(eta_hat, eta_true),
  mse_mu = mse(mu_hat, mu_true),
  rmse_mu = rmse(mu_hat, mu_true),
  mae_mu = mae(mu_hat, mu_true),
  mass_err_max = mass_err_max,
  mass_err_l2 = mass_err_l2,
  deviance_coarse = dev_coarse,
  loglik_coarse = ll_coarse,
  mse_f = mse(f_hat, f_true),
  mse_g = mse(g_hat, g_true),
  mse_h = mse(h_hat, h_true)
)

tab <- rbind(row_A, row_B, row_C)

out_dir <- "experiments/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(tab, file.path(out_dir, "recovery_ABC_metrics_single.csv"), row.names = FALSE)
saveRDS(
  list(
    table = tab,
    note = "Single-replica metrics. For paper tables, repeat over many seeds and report Monte Carlo mean/SD."
  ),
  file.path(out_dir, "recovery_ABC_metrics_single.rds")
)

cat("Wrote:\n",
    file.path(out_dir, "recovery_ABC_metrics_single.csv"), "\n",
    file.path(out_dir, "recovery_ABC_metrics_single.rds"), "\n", sep = "")


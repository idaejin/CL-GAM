# Shared evaluation metrics for nested ATA / CL-GAMM competitors.
#
# Metric roles (see improvements/COMPARE_METRICS.md):
#   Coarse (primary on real ATA): loglik/deviance, AIC/BIC (same family), mass error
#   Fine diagnostic (if y_fine known): mse/mae vs oracle log-rate, rmse counts
#   Uncertainty: nominal 95% coverage vs oracle when sd.eta available
#   Predictive: leave-one-coarse-unit-out Poisson log-score (optional)

#' Fine-scale oracle log-rate
clgam_oracle_lograte <- function(y_fine, e_fine, eps = 0.5) {
  log((as.numeric(y_fine) + eps) / pmax(as.numeric(e_fine), .Machine$double.xmin))
}

#' Fitted fine means from eta
clgam_mu_from_eta <- function(eta, e_fine) {
  as.numeric(e_fine) * exp(as.numeric(eta))
}

#' Coarse means mu_c = C (e * exp(eta))
clgam_mu_coarse <- function(eta, e_fine, C) {
  as.numeric(as.matrix(C) %*% clgam_mu_from_eta(eta, e_fine))
}

#' Poisson log-likelihood for coarse counts under composite link
clgam_poisson_loglik <- function(y, mu, eps = .Machine$double.xmin) {
  y <- as.numeric(y)
  mu <- pmax(as.numeric(mu), eps)
  sum(y * log(mu) - mu - lgamma(y + 1))
}

#' Poisson deviance (saturated vs fitted)
clgam_poisson_deviance <- function(y, mu, eps = .Machine$double.xmin) {
  y <- as.numeric(y)
  mu <- pmax(as.numeric(mu), eps)
  sat <- ifelse(y > 0, y * log(y) - y, 0)
  fit <- y * log(mu) - mu
  2 * sum(sat - fit)
}

#' Nested ATA full scorecard for one method
#'
#' @param sd_eta optional SE for eta (coverage vs oracle when supplied)
#' @param aic,bic,ed optional; AIC/BIC only comparable within same likelihood family
#' @param family_label character tag for AIC comparability group
clgam_score_ata <- function(
  eta,
  y_fine,
  e_fine,
  y_coarse = NULL,
  C = NULL,
  aic = NA_real_,
  bic = NA_real_,
  ed = NA_real_,
  sd_eta = NULL,
  family_label = NA_character_,
  method = "model"
) {
  eta <- as.numeric(eta)
  y_fine <- as.numeric(y_fine)
  e_fine <- as.numeric(e_fine)
  stopifnot(length(eta) == length(y_fine), length(e_fine) == length(y_fine))
  oracle <- clgam_oracle_lograte(y_fine, e_fine)
  mu <- clgam_mu_from_eta(eta, e_fine)

  mass_max <- mass_l2 <- loglik_c <- deviance_c <- NA_real_
  if (!is.null(C) && !is.null(y_coarse)) {
    mu_c <- clgam_mu_coarse(eta, e_fine, C)
    y_c <- as.numeric(y_coarse)
    resid <- mu_c - y_c
    mass_max <- max(abs(resid))
    mass_l2 <- sqrt(sum(resid^2))
    loglik_c <- clgam_poisson_loglik(y_c, mu_c)
    deviance_c <- clgam_poisson_deviance(y_c, mu_c)
  }

  cover <- mean_se <- NA_real_
  if (!is.null(sd_eta)) {
    sd_eta <- as.numeric(sd_eta)
    stopifnot(length(sd_eta) == length(eta))
    cover <- mean(abs(eta - oracle) <= 1.96 * sd_eta, na.rm = TRUE)
    mean_se <- mean(sd_eta, na.rm = TRUE)
  }

  data.frame(
    method = method,
    family = family_label,
    # --- coarse (primary on real ATA) ---
    loglik_coarse = loglik_c,
    deviance_coarse = deviance_c,
    aic = aic,
    bic = bic,
    ed = ed,
    mass_err_max = mass_max,
    mass_err_l2 = mass_l2,
    mass_err = mass_max, # alias for older scripts
    # --- fine diagnostic vs oracle ---
    mse_eta = mean((eta - oracle)^2),
    mae_eta = mean(abs(eta - oracle)),
    cor_eta = suppressWarnings(cor(eta, oracle)),
    rmse_counts = sqrt(mean((mu - y_fine)^2)),
    mae_counts = mean(abs(mu - y_fine)),
    # --- uncertainty vs noisy oracle ---
    cover95_oracle = cover,
    mean_se_eta = mean_se,
    # --- truth metrics (sims only; NA here) ---
    mse_eta_truth = NA_real_,
    cover95_truth = NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Baseline: constant log-rate within each coarse unit = log(y_i / e_i)
clgam_baseline_region_smr <- function(y_coarse, e_fine, C) {
  C <- as.matrix(C)
  e_c <- as.numeric(C %*% e_fine)
  eta_c <- log(pmax(as.numeric(y_coarse), 0.5) / pmax(e_c, .Machine$double.xmin))
  as.numeric(t(C) %*% eta_c)
}

#' Leave-one-coarse-unit-out Poisson log-score
#'
#' @param fit_fun function(y_train, C_train) -> list(eta = fine eta of length m)
#'   Must predict eta on the full fine grid using only training coarse rows.
#' @return data.frame with loo_logscore (sum), loo_logscore_mean, loo_mae_y, loo_rmse_y
clgam_loo_coarse_score <- function(y_coarse, e_fine, C, fit_fun, trace = FALSE) {
  C <- as.matrix(C)
  y_coarse <- as.numeric(y_coarse)
  e_fine <- as.numeric(e_fine)
  n <- length(y_coarse)
  stopifnot(nrow(C) == n, ncol(C) == length(e_fine))
  ll <- mae <- se <- numeric(n)
  for (i in seq_len(n)) {
    if (trace) message("  LOO fold ", i, "/", n)
    keep <- setdiff(seq_len(n), i)
    # Drop fine cells that belong only to held-out unit (exclusive partition)
    fine_train <- which(as.numeric(C[i, ]) == 0 | colSums(C[keep, , drop = FALSE]) > 0)
    # For exclusive 0/1 partition: train on counties not in region i
    fine_in_i <- which(C[i, ] > 0)
    fine_train <- setdiff(seq_len(ncol(C)), fine_in_i)
    if (length(fine_train) < 2L || length(keep) < 2L) {
      ll[i] <- NA_real_
      mae[i] <- NA_real_
      se[i] <- NA_real_
      next
    }
    C_tr <- C[keep, fine_train, drop = FALSE]
    y_tr <- y_coarse[keep]
    e_tr <- e_fine[fine_train]
    fit <- fit_fun(y_tr, e_tr, C_tr, fine_train)
    eta_full <- rep(NA_real_, ncol(C))
    eta_full[fine_train] <- as.numeric(fit$eta)
    # Predict held-out region: use spatial/cov model extrapolated to fine_in_i
    if (!is.null(fit$eta_holdout)) {
      eta_full[fine_in_i] <- as.numeric(fit$eta_holdout)
    } else if (!is.null(fit$predict_eta)) {
      eta_full[fine_in_i] <- as.numeric(fit$predict_eta(fine_in_i))
    } else {
      # Fallback: region-mean of neighbouring training etas (weak)
      eta_full[fine_in_i] <- mean(fit$eta, na.rm = TRUE)
    }
    mu_i <- sum(e_fine[fine_in_i] * exp(eta_full[fine_in_i]))
    ll[i] <- clgam_poisson_loglik(y_coarse[i], mu_i)
    mae[i] <- abs(mu_i - y_coarse[i])
    se[i] <- (mu_i - y_coarse[i])^2
  }
  data.frame(
    loo_logscore = sum(ll, na.rm = TRUE),
    loo_logscore_mean = mean(ll, na.rm = TRUE),
    loo_mae_y = mean(mae, na.rm = TRUE),
    loo_rmse_y = sqrt(mean(se, na.rm = TRUE)),
    loo_n = sum(is.finite(ll)),
    stringsAsFactors = FALSE
  )
}

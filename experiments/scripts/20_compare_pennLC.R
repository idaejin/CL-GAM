#!/usr/bin/env Rscript
# Full metric suite on pennLC (20 k-means regions -> 67 counties).
#
# Coarse: loglik, deviance, AIC/BIC/ed (pois_SAP), mass max/L2, LOO log-score
# Fine diagnostic: mse/mae/cor vs oracle, rmse/mae counts
# Uncertainty: 95% coverage of eta vs noisy oracle (sd.eta)
# Truth metrics: NA (no known eta)
#
# From experiments/:  Rscript scripts/20_compare_pennLC.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/03_malone.R"))
source(file.path(root, "R/07_load_pennLC.R"))
source(file.path(root, "R/08_metrics.R"))

suppressPackageStartupMessages(library(ggplot2))

dat <- clgam_load_pennLC(n_coarse = 20L, seed = 1L)
penn <- readRDS(file.path(CLGAM_OUTPUT, "penn_case_ABC_fit.rds"))
ndx <- penn$ndx
ndxnl <- penn$ndxnl
oracle <- clgam_oracle_lograte(dat$y_fine, dat$e_fine)
smk <- as.numeric(dat$covariates_fine[, "smoking"])
north <- as.numeric(dat$covariates_agg_expanded[, 1])

idw_eta <- function(x1, x2, eta_src, x1_src, x2_src, k = 5L) {
  k <- min(as.integer(k), length(eta_src))
  vapply(seq_along(x1), function(j) {
    d <- (x1_src - x1[j])^2 + (x2_src - x2[j])^2
    o <- order(d)[seq_len(k)]
    w <- 1 / pmax(d[o], 1e-12)
    sum(w * eta_src[o]) / sum(w)
  }, numeric(1))
}

message("=== Refitting pennLC models with sd.eta ===")
fit_sp <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_fine,
  C = dat$C, ndx = ndx, elements = TRUE, trace = FALSE
)
fit_a <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_fine,
  nlcovfine = dat$covariates_fine, C = dat$C,
  ndx = ndx, ndxnl = ndxnl, elements = TRUE, trace = FALSE
)
fit_b <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_fine,
  lcovfine = dat$covariates_agg_expanded, C = dat$C,
  ndx = ndx, elements = TRUE, trace = FALSE
)
fit_c <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_fine,
  lcovfine = dat$covariates_agg_expanded,
  nlcovfine = dat$covariates_fine, C = dat$C,
  ndx = ndx, ndxnl = ndxnl, elements = TRUE, trace = FALSE
)
fit_mal <- clgam_malone_fit(
  y_mun = dat$y,
  e_fine = dat$e_fine,
  C = dat$C,
  cov_fine = data.frame(smoking = smk),
  lon = dat$x1,
  lat = dat$x2,
  cov_smooth = "smoking",
  spatial = "none",
  family = "poisson",
  round_counts = TRUE,
  init = "equal",
  deliver = "fitted",
  stop_on = "map",
  max_iter = if (CLGAM_FAST) 40L else 200L,
  tol = 1e-3,
  k_cov = as.integer(ndxnl),
  trace = FALSE
)
eta_base <- clgam_baseline_region_smr(dat$y, dat$e_fine, dat$C)

# Malone pointwise SE from final GAM (terms+lp)
sd_mal <- tryCatch({
  se <- as.numeric(predict(fit_mal$model, type = "link", se.fit = TRUE)$se.fit)
  se
}, error = function(e) NULL)

models <- list(
  list(
    method = "Baseline region-SMR", family = "baseline",
    eta = eta_base, sd = NULL, aic = NA, bic = NA, ed = NA
  ),
  list(
    method = "ATA spatial", family = "pois_SAP",
    eta = fit_sp$eta, sd = fit_sp$sd.eta,
    aic = fit_sp$aic, bic = fit_sp$bic, ed = fit_sp$ed
  ),
  list(
    method = "CL-GAMM Case A", family = "pois_SAP",
    eta = fit_a$eta, sd = fit_a$sd.eta,
    aic = fit_a$aic, bic = fit_a$bic, ed = fit_a$ed
  ),
  list(
    method = "CL-GAMM Case B", family = "pois_SAP",
    eta = fit_b$eta, sd = fit_b$sd.eta,
    aic = fit_b$aic, bic = fit_b$bic, ed = fit_b$ed
  ),
  list(
    method = "CL-GAMM Case C", family = "pois_SAP",
    eta = fit_c$eta, sd = fit_c$sd.eta,
    aic = fit_c$aic, bic = fit_c$bic, ed = fit_c$ed
  ),
  list(
    method = "Malone Method 1", family = "malone_gam",
    eta = fit_mal$eta, sd = sd_mal,
    aic = tryCatch(AIC(fit_mal$model), error = function(e) NA_real_),
    bic = tryCatch(BIC(fit_mal$model), error = function(e) NA_real_),
    ed = tryCatch(sum(fit_mal$model$edf), error = function(e) NA_real_)
  )
)

tab <- do.call(rbind, lapply(models, function(m) {
  clgam_score_ata(
    m$eta, dat$y_fine, dat$e_fine, dat$y, dat$C,
    aic = m$aic, bic = m$bic, ed = m$ed, sd_eta = m$sd,
    family_label = m$family, method = m$method
  )
}))

# --------------------------------------------------------------------------
# LOO coarse predictive scores
# --------------------------------------------------------------------------
message("=== LOO coarse (20 folds; holdout eta via IDW of training surface) ===")
C <- as.matrix(dat$C)
n <- nrow(C)

loo_one <- function(label, fit_train) {
  ll <- mae <- se2 <- numeric(n)
  for (i in seq_len(n)) {
    fine_i <- which(C[i, ] > 0)
    fine_tr <- setdiff(seq_len(ncol(C)), fine_i)
    keep <- setdiff(seq_len(n), i)
    ft <- fit_train(keep, fine_tr)
    eta_i <- idw_eta(
      dat$x1[fine_i], dat$x2[fine_i], ft$eta,
      dat$x1[fine_tr], dat$x2[fine_tr], k = 5L
    )
    mu_i <- sum(dat$e_fine[fine_i] * exp(eta_i))
    ll[i] <- clgam_poisson_loglik(dat$y[i], mu_i)
    mae[i] <- abs(mu_i - dat$y[i])
    se2[i] <- (mu_i - dat$y[i])^2
  }
  data.frame(
    method = label,
    loo_logscore = sum(ll),
    loo_logscore_mean = mean(ll),
    loo_mae_y = mean(mae),
    loo_rmse_y = sqrt(mean(se2)),
    stringsAsFactors = FALSE
  )
}

loo_rows <- list(
  loo_one("Baseline region-SMR", function(keep, fine_tr) {
    list(eta = clgam_baseline_region_smr(dat$y[keep], dat$e_fine[fine_tr], C[keep, fine_tr, drop = FALSE]))
  }),
  loo_one("ATA spatial", function(keep, fine_tr) {
    f <- pois_SAP(
      y = dat$y[keep], x1 = dat$x1[fine_tr], x2 = dat$x2[fine_tr],
      efine = dat$e_fine[fine_tr], C = C[keep, fine_tr, drop = FALSE],
      ndx = ndx, elements = FALSE, trace = FALSE
    )
    list(eta = f$eta)
  }),
  loo_one("CL-GAMM Case A", function(keep, fine_tr) {
    f <- pois_SAP(
      y = dat$y[keep], x1 = dat$x1[fine_tr], x2 = dat$x2[fine_tr],
      efine = dat$e_fine[fine_tr],
      nlcovfine = matrix(smk[fine_tr], ncol = 1, dimnames = list(NULL, "smoking")),
      C = C[keep, fine_tr, drop = FALSE],
      ndx = ndx, ndxnl = ndxnl, elements = FALSE, trace = FALSE
    )
    list(eta = f$eta)
  }),
  loo_one("CL-GAMM Case B", function(keep, fine_tr) {
    f <- pois_SAP(
      y = dat$y[keep], x1 = dat$x1[fine_tr], x2 = dat$x2[fine_tr],
      efine = dat$e_fine[fine_tr],
      lcovfine = matrix(north[fine_tr], ncol = 1, dimnames = list(NULL, "north_agg")),
      C = C[keep, fine_tr, drop = FALSE],
      ndx = ndx, elements = FALSE, trace = FALSE
    )
    list(eta = f$eta)
  }),
  loo_one("CL-GAMM Case C", function(keep, fine_tr) {
    f <- pois_SAP(
      y = dat$y[keep], x1 = dat$x1[fine_tr], x2 = dat$x2[fine_tr],
      efine = dat$e_fine[fine_tr],
      lcovfine = matrix(north[fine_tr], ncol = 1, dimnames = list(NULL, "north_agg")),
      nlcovfine = matrix(smk[fine_tr], ncol = 1, dimnames = list(NULL, "smoking")),
      C = C[keep, fine_tr, drop = FALSE],
      ndx = ndx, ndxnl = ndxnl, elements = FALSE, trace = FALSE
    )
    list(eta = f$eta)
  }),
  loo_one("Malone Method 1", function(keep, fine_tr) {
    f <- clgam_malone_fit(
      y_mun = dat$y[keep],
      e_fine = dat$e_fine[fine_tr],
      C = C[keep, fine_tr, drop = FALSE],
      cov_fine = data.frame(smoking = smk[fine_tr]),
      lon = dat$x1[fine_tr],
      lat = dat$x2[fine_tr],
      cov_smooth = "smoking",
      spatial = "none",
      family = "poisson",
      round_counts = TRUE,
      init = "equal",
      deliver = "fitted",
      stop_on = "map",
      max_iter = 40L,
      tol = 1e-3,
      k_cov = as.integer(ndxnl),
      trace = FALSE
    )
    list(eta = f$eta)
  })
)

tab_loo <- do.call(rbind, loo_rows)
tab <- merge(tab, tab_loo, by = "method", all.x = TRUE, sort = FALSE)
tab <- tab[order(tab$deviance_coarse), ]
rownames(tab) <- NULL

# Structural ceiling note
reg <- apply(C, 2, function(z) which(z > 0)[1])
eta_ceil <- as.numeric(tapply(oracle, reg, mean)[as.character(reg)])
ceiling_cor <- cor(oracle, eta_ceil)
ss_tot <- sum((oracle - mean(oracle))^2)
ss_between <- sum((eta_ceil - mean(oracle))^2)
within_share <- 1 - ss_between / ss_tot

write.csv(tab, file.path(CLGAM_OUTPUT, "penn_compare_metrics.csv"), row.names = FALSE)
saveRDS(
  list(
    table = tab,
    design = list(n_coarse = 20L, n_fine = 67L, seed = 1L),
    notes = c(
      "mse_eta_truth / cover95_truth = NA (no known true eta on pennLC).",
      "cover95_oracle uses noisy county log-SMR; expect undercoverage if SE calibrated for true eta.",
      "AIC/BIC comparable only within family == pois_SAP.",
      "loglik_coarse / deviance_coarse are Poisson composite-link scores (fair across methods).",
      "LOO holdout eta via IDW(k=5) of training-county eta (approximation; not native predict).",
      sprintf("Structural cor ceiling (region-mean oracle): %.3f; within-region var share: %.3f",
              ceiling_cor, within_share)
    )
  ),
  file.path(CLGAM_OUTPUT, "penn_compare_metrics.rds")
)

# Console report
cols_primary <- c(
  "method", "loglik_coarse", "deviance_coarse", "aic", "bic", "ed",
  "mass_err_max", "mass_err_l2",
  "mse_eta", "mae_eta", "cor_eta", "rmse_counts", "mae_counts",
  "cover95_oracle", "mean_se_eta",
  "loo_logscore", "loo_logscore_mean", "loo_mae_y", "loo_rmse_y"
)
message("\n=== pennLC full metrics (sorted by deviance_coarse) ===")
print(tab[, cols_primary], digits = 4)

# Plot: three panels
tab$method_f <- factor(tab$method, levels = rev(tab$method))
p1 <- ggplot(tab, aes(method_f, deviance_coarse)) +
  geom_col(fill = "#4C78A8", width = 0.7) + coord_flip() +
  labs(title = "Coarse Poisson deviance", x = NULL, y = "deviance") +
  theme_bw(base_size = 11)
p2 <- ggplot(tab, aes(method_f, mse_eta)) +
  geom_col(fill = "#54A24B", width = 0.7) + coord_flip() +
  labs(title = "Fine MSE vs oracle", x = NULL, y = "mse_eta") +
  theme_bw(base_size = 11)
p3 <- ggplot(tab, aes(method_f, loo_logscore)) +
  geom_col(fill = "#F58518", width = 0.7) + coord_flip() +
  labs(title = "LOO coarse log-score (higher better)", x = NULL, y = "sum loglik") +
  theme_bw(base_size = 11)
p4 <- ggplot(subset(tab, family == "pois_SAP"), aes(reorder(method, -aic), aic)) +
  geom_col(fill = "#E45756", width = 0.7) + coord_flip() +
  labs(title = "AIC (pois_SAP only)", x = NULL, y = "AIC") +
  theme_bw(base_size = 11)

png(file.path(CLGAM_OUTPUT, "penn_compare_metrics.png"), width = 1600, height = 1100, res = 140)
if (requireNamespace("patchwork", quietly = TRUE)) {
  print(
    (p1 | p2) / (p3 | p4) +
      patchwork::plot_annotation(
        title = "pennLC full metric suite",
        subtitle = "Coarse deviance/LOO primary for real ATA; fine MSE diagnostic; AIC within pois_SAP only"
      )
  )
} else {
  print(p1)
}
dev.off()

message("Wrote ", file.path(CLGAM_OUTPUT, "penn_compare_metrics.csv"))
message("Wrote ", file.path(CLGAM_OUTPUT, "penn_compare_metrics.png"))

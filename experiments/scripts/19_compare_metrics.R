#!/usr/bin/env Rscript
# Compare CL-GAMM vs competitors with a common metric suite.
#
# Datasets:
#   - Madrid MEDEA (existing RDS): ATA spatial, Case A, Malone, ST-PCLM ATA
#   - pennLC open ATA: Case A/B/C + Malone + region-SMR baseline + spatial-only
#
# Primary metric: MSE of fine eta vs oracle log((y+0.5)/e)
# Secondary: cor_eta, rmse_counts, mass_err; AIC only within pois_SAP family
#
# From experiments/:  Rscript scripts/19_compare_metrics.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))
source(file.path(root, "R/03_malone.R"))
source(file.path(root, "R/07_load_pennLC.R"))
source(file.path(root, "R/08_metrics.R"))

suppressPackageStartupMessages(library(ggplot2))

# --------------------------------------------------------------------------
# Madrid
# --------------------------------------------------------------------------
message("=== Madrid MEDEA metrics ===")
dat_m <- clgam_load_madrid(load_maps = FALSE)
oracle_m <- clgam_oracle_lograte(dat_m$yc, dat_m$ec)

ata <- readRDS(file.path(CLGAM_OUTPUT, "ATA_spatial_fit.rds"))
casea <- readRDS(file.path(CLGAM_OUTPUT, "CaseA_model25_fit.rds"))
mal <- readRDS(file.path(CLGAM_OUTPUT, "malone_corrected_fit.rds"))
stp <- tryCatch(
  readRDS(file.path(CLGAM_OUTPUT, "st_pclm_ata_nested.rds")),
  error = function(e) NULL
)

madrid_rows <- list(
  clgam_score_ata(
    ata$model1$eta, dat_m$yc, dat_m$ec, dat_m$ym, dat_m$C_m,
    aic = ata$model1$aic, method = "ATA spatial (pois_SAP)"
  ),
  clgam_score_ata(
    casea$fit$eta, dat_m$yc, dat_m$ec, dat_m$ym, dat_m$C_m,
    aic = casea$fit$aic, method = "CL-GAMM Case A (pois_SAP)"
  ),
  clgam_score_ata(
    mal$eta, dat_m$yc, dat_m$ec, dat_m$ym, dat_m$C_m,
    aic = mal$aic, method = "Malone Method 1"
  )
)
if (!is.null(stp) && !is.null(stp$fit_ata$eta)) {
  madrid_rows[[length(madrid_rows) + 1L]] <- clgam_score_ata(
    stp$fit_ata$eta, dat_m$yc, dat_m$ec, dat_m$ym, dat_m$C_m,
    aic = stp$fit_ata$aic, method = "ST-PCLM ATA nested (pois_SAP)"
  )
}

# Region-SMR baseline at mun->CT
eta_base_m <- clgam_baseline_region_smr(dat_m$ym, dat_m$ec, dat_m$C_m)
madrid_rows[[length(madrid_rows) + 1L]] <- clgam_score_ata(
  eta_base_m, dat_m$yc, dat_m$ec, dat_m$ym, dat_m$C_m,
  method = "Baseline: mun SMR expanded"
)

tab_madrid <- do.call(rbind, madrid_rows)
tab_madrid$dataset <- "Madrid"
tab_madrid <- tab_madrid[order(tab_madrid$mse_eta), ]
rownames(tab_madrid) <- NULL
message("Madrid (ranked by mse_eta):")
print(tab_madrid[, c("method", "mse_eta", "cor_eta", "rmse_counts", "mass_err", "aic")])

# --------------------------------------------------------------------------
# pennLC
# --------------------------------------------------------------------------
message("=== pennLC metrics + competitors ===")
dat_p <- clgam_load_pennLC(n_coarse = 20L, seed = 1L)
penn <- readRDS(file.path(CLGAM_OUTPUT, "penn_case_ABC_fit.rds"))
ndx <- penn$ndx
ndxnl <- penn$ndxnl

# Spatial-only
message("Fitting pennLC spatial-only...")
fit_sp <- pois_SAP(
  y = dat_p$y, x1 = dat_p$x1, x2 = dat_p$x2, efine = dat_p$e_fine,
  C = dat_p$C, ndx = ndx, elements = TRUE, trace = FALSE
)

# Malone Method 1 with smoking
message("Fitting pennLC Malone Method 1...")
cov_p <- data.frame(smoking = as.numeric(dat_p$covariates_fine[, "smoking"]))
fit_mal_p <- clgam_malone_fit(
  y_mun = dat_p$y,
  e_fine = dat_p$e_fine,
  C = dat_p$C,
  cov_fine = cov_p,
  lon = dat_p$x1,
  lat = dat_p$x2,
  cov_smooth = "smoking",
  spatial = "none",
  family = "poisson",
  round_counts = TRUE,
  init = "equal",
  deliver = "fitted",
  stop_on = "map",
  max_iter = if (CLGAM_FAST) 40L else 200L,
  tol = 1e-3,
  k_cov = if (is.na(ndxnl)) 8L else as.integer(ndxnl),
  trace = FALSE
)

eta_base_p <- clgam_baseline_region_smr(dat_p$y, dat_p$e_fine, dat_p$C)

penn_rows <- list(
  clgam_score_ata(
    fit_sp$eta, dat_p$y_fine, dat_p$e_fine, dat_p$y, dat_p$C,
    aic = fit_sp$aic, method = "ATA spatial (pois_SAP)"
  ),
  clgam_score_ata(
    penn$case_a$eta, dat_p$y_fine, dat_p$e_fine, dat_p$y, dat_p$C,
    aic = penn$case_a$aic, method = "CL-GAMM Case A s(smoking)"
  ),
  clgam_score_ata(
    penn$case_b$eta, dat_p$y_fine, dat_p$e_fine, dat_p$y, dat_p$C,
    aic = penn$case_b$aic, method = "CL-GAMM Case B north_agg"
  ),
  clgam_score_ata(
    penn$case_c$eta, dat_p$y_fine, dat_p$e_fine, dat_p$y, dat_p$C,
    aic = penn$case_c$aic, method = "CL-GAMM Case C"
  ),
  clgam_score_ata(
    fit_mal_p$eta, dat_p$y_fine, dat_p$e_fine, dat_p$y, dat_p$C,
    aic = tryCatch(AIC(fit_mal_p$model), error = function(e) NA_real_),
    method = "Malone Method 1"
  ),
  clgam_score_ata(
    eta_base_p, dat_p$y_fine, dat_p$e_fine, dat_p$y, dat_p$C,
    method = "Baseline: region SMR expanded"
  )
)

tab_penn <- do.call(rbind, penn_rows)
tab_penn$dataset <- "pennLC"
tab_penn <- tab_penn[order(tab_penn$mse_eta), ]
rownames(tab_penn) <- NULL
message("pennLC (ranked by mse_eta):")
print(tab_penn[, c("method", "mse_eta", "cor_eta", "rmse_counts", "mass_err", "aic")])

# --------------------------------------------------------------------------
# Save + plot
# --------------------------------------------------------------------------
tab_all <- rbind(tab_madrid, tab_penn)
write.csv(tab_all, file.path(CLGAM_OUTPUT, "compare_metrics.csv"), row.names = FALSE)
saveRDS(
  list(
    madrid = tab_madrid,
    pennLC = tab_penn,
    primary_metric = "mse_eta = mean((eta - log((y_fine+0.5)/e_fine))^2)",
    notes = c(
      "Lower mse_eta is better (fine-scale vs oracle).",
      "AIC comparable only among pois_SAP rows on the same dataset.",
      "Malone AIC is from final heuristic GAM (not composite-link).",
      "Baseline expands coarse SMR as constant within each coarse unit.",
      "disaggregation (Madrid raster) omitted here: cell predictions need CT aggregation bridge."
    )
  ),
  file.path(CLGAM_OUTPUT, "compare_metrics.rds")
)

plot_rank <- function(tab, title) {
  tab$method <- factor(tab$method, levels = rev(tab$method))
  ggplot(tab, aes(method, mse_eta)) +
    geom_col(fill = "#4C78A8", width = 0.7) +
    coord_flip() +
    labs(
      title = title,
      subtitle = "Primary metric: MSE of fine eta vs oracle log-rate",
      x = NULL, y = "mse_eta"
    ) +
    theme_bw(base_size = 11)
}

png(file.path(CLGAM_OUTPUT, "compare_metrics.png"), width = 1400, height = 900, res = 140)
if (requireNamespace("patchwork", quietly = TRUE)) {
  print(
    plot_rank(tab_madrid, "Madrid MEDEA") / plot_rank(tab_penn, "pennLC (20 regions -> 67 counties)") +
      patchwork::plot_annotation(
        title = "Competitor comparison by mse_eta",
        subtitle = "Lower is better; oracle uses fine counts not seen by ATA fits"
      )
  )
} else {
  print(plot_rank(tab_madrid, "Madrid"))
}
dev.off()

message("Wrote ", file.path(CLGAM_OUTPUT, "compare_metrics.csv"))
message("Wrote ", file.path(CLGAM_OUTPUT, "compare_metrics.rds"))
message("Wrote ", file.path(CLGAM_OUTPUT, "compare_metrics.png"))

#!/usr/bin/env Rscript
# Figures for Madrid nested ATA / competitors (maps + partial effects ± SE).
#
# From experiments/:  Rscript scripts/12_plot_madrid_comparison.R
# Maps reuse saved RDS. Smooth bands always use paper bases (ndx=20, ndxnl=12);
# do not use CLGAM_FAST for partials (FAST warps Case A s(unemployed)).

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
})

dat <- clgam_load_madrid(load_maps = TRUE)
mun <- dat$map_mun$sf
ct <- dat$map_ct$sf
# align tidy order: mun codes via C like disaggregation helper
source(file.path(root, "R/04_disaggregation.R"))
mun_ord <- clgam_mun_sf_for_disag(dat)
ct_ord <- dat$map_ct$sf
stopifnot(nrow(ct_ord) == length(dat$yc))

ata <- readRDS(file.path(CLGAM_OUTPUT, "ATA_spatial_fit.rds"))
casea <- readRDS(file.path(CLGAM_OUTPUT, "CaseA_model25_fit.rds"))
stp <- readRDS(file.path(CLGAM_OUTPUT, "st_pclm_ata_nested.rds"))
mal <- readRDS(file.path(CLGAM_OUTPUT, "malone_corrected_fit.rds"))

eta_ata_spatial <- ata$model1$eta
eta_casea <- casea$fit$eta
eta_stp <- stp$fit_ata$eta
eta_mal <- mal$eta
logSMR_mun <- dat$logSMR_mun
logSMR_ct <- log((dat$yc + 1) / dat$ec)

# shared colour breaks (deciles of ATA spatial eta)
brks <- quantile(eta_ata_spatial, probs = seq(0, 1, 0.1), na.rm = TRUE)
brks[1] <- brks[1] - 1e-6
brks[length(brks)] <- brks[length(brks)] + 1e-6

cut_eta <- function(x) cut(x, breaks = brks, include.lowest = TRUE)

theme_map <- theme_void(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey30")
  )

pal <- scales::viridis_pal(option = "C")(10)

map_ct <- function(values, title, subtitle = NULL) {
  d <- ct_ord
  d$val <- values
  d$bin <- cut_eta(values)
  ggplot(d) +
    geom_sf(aes(fill = bin), colour = NA) +
    scale_fill_manual(values = pal, drop = FALSE, name = "log-rate") +
    labs(title = title, subtitle = subtitle) +
    theme_map
}

map_mun <- function(values, title, subtitle = NULL) {
  d <- mun_ord
  d$val <- values
  d$bin <- cut_eta(values)
  ggplot(d) +
    geom_sf(aes(fill = bin), colour = "grey40", linewidth = 0.05) +
    scale_fill_manual(values = pal, drop = FALSE, name = "log-rate") +
    labs(title = title, subtitle = subtitle) +
    theme_map
}

p1 <- map_mun(logSMR_mun, "Raw log(SMR) municipalities", "ym observed / em")
p2 <- map_ct(eta_ata_spatial, "ATA spatial (no covariates)", sprintf("pois_SAP  AIC=%.1f", ata$model1$aic))
p3 <- map_ct(eta_casea, "CL-GAMM Case A", sprintf("ageing + unemployed  AIC=%.1f", casea$fit$aic))
p4 <- map_ct(eta_stp, "ST-PCLM spatial ATA nested", sprintf("same Case A structure  AIC=%.1f", stp$fit_ata$aic))
p5 <- map_ct(eta_mal, "Malone heuristic", sprintf("Method 1 / poisson+round  iters=%s", mal$n_iter))
p6 <- map_ct(logSMR_ct, "Oracle CT log(SMR)", "(yc+1)/ec - not used in mun-only fits")

# Difference maps vs oracle
diff_casea <- eta_casea - logSMR_ct
diff_mal <- eta_mal - logSMR_ct
lim <- max(abs(c(diff_casea, diff_mal)), na.rm = TRUE)

map_diff <- function(values, title) {
  d <- ct_ord
  d$val <- values
  ggplot(d) +
    geom_sf(aes(fill = val), colour = NA) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = c(-lim, lim), name = "delta log-rate"
    ) +
    labs(title = title, subtitle = "fit eta - oracle log(SMR)") +
    theme_map
}

p7 <- map_diff(diff_casea, "Case A - oracle")
p8 <- map_diff(diff_mal, "Malone - oracle")

# Scatter eta vs oracle
sc <- data.frame(
  oracle = logSMR_ct,
  CaseA = eta_casea,
  ATA = eta_ata_spatial,
  Malone = eta_mal
)
sc_long <- rbind(
  data.frame(oracle = sc$oracle, fit = sc$CaseA, method = "Case A"),
  data.frame(oracle = sc$oracle, fit = sc$ATA, method = "ATA spatial"),
  data.frame(oracle = sc$oracle, fit = sc$Malone, method = "Malone")
)
p_sc <- ggplot(sc_long, aes(oracle, fit)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) +
  geom_point(alpha = 0.15, size = 0.6) +
  facet_wrap(~method) +
  labs(
    title = "Fine-scale eta vs oracle CT log(SMR)",
    x = "oracle log((yc+1)/ec)", y = "fitted eta"
  ) +
  theme_bw(base_size = 11) +
  coord_equal()

# --- Univariate smooths with 95% bands ---
# Case A: pois_SAP nleffects ± 1.96 sdnleffects (elements=TRUE)
# Malone: mgcv predict(..., type="terms", se.fit=TRUE)

smooth_band_df <- function(x_pct, fit, se, name) {
  # centre for display (same constant on fit and bands)
  ctr <- mean(fit, na.rm = TRUE)
  o <- order(x_pct)
  data.frame(
    x = x_pct[o],
    fit = fit[o] - ctr,
    lo = fit[o] - 1.96 * se[o] - ctr,
    hi = fit[o] + 1.96 * se[o] - ctr,
    se = se[o],
    halfwidth = 1.96 * se[o],
    term = name
  )
}

agg_band <- function(df) {
  aggregate(cbind(fit, lo, hi, se, halfwidth) ~ x, data = df, FUN = mean)
}

cols_method <- c("CL-GAMM Case A" = "#2166ac", "Malone" = "#b2182b")

plot_smooth_band <- function(df, title, xlab, colour, ylim = NULL) {
  ag <- agg_band(df)
  p <- ggplot() +
    geom_ribbon(
      data = ag, aes(x = x, ymin = lo, ymax = hi),
      fill = colour, alpha = 0.28
    ) +
    geom_line(data = ag, aes(x = x, y = fit), colour = colour, linewidth = 0.95) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = 3, linewidth = 0.4) +
    labs(
      title = title,
      subtitle = sprintf(
        "centred +/- 1.96 SE  |  median band width = %.3f",
        median(ag$hi - ag$lo)
      ),
      x = xlab, y = "partial effect (log-rate)"
    ) +
    theme_bw(base_size = 11)
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# Overlay both methods on a shared y-scale (emphasises magnitude + band contrast)
plot_overlay_bands <- function(df_ca, df_mal, title, xlab) {
  ca <- agg_band(df_ca)
  ca$method <- "CL-GAMM Case A"
  mal <- agg_band(df_mal)
  mal$method <- "Malone"
  both <- rbind(ca, mal)
  both$method <- factor(both$method, levels = names(cols_method))
  med_w <- tapply(both$hi - both$lo, both$method, median)
  ggplot(both, aes(x = x, colour = method, fill = method)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, colour = NA) +
    geom_line(aes(y = fit), linewidth = 1) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = 3, linewidth = 0.4) +
    scale_colour_manual(values = cols_method) +
    scale_fill_manual(values = cols_method) +
    labs(
      title = title,
      subtitle = sprintf(
        "median 95%% width Case A=%.3f / Malone=%.3f (~%.0fx)",
        med_w[["CL-GAMM Case A"]], med_w[["Malone"]],
        med_w[["CL-GAMM Case A"]] / med_w[["Malone"]]
      ),
      x = xlab, y = "partial effect (log-rate)",
      colour = NULL, fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

# Partials MUST use paper bases (ndx=20, ndxnl=12). FAST (8/6) warps
# s(unemployed) into an artifactual spike (~±4) and inflates SE contrast.
message("Case A refit (elements) for smooths + SE — paper ndx/ndxnl (not FAST)...")
ndx <- c(20L, 20L)
ndxnl <- 12L
cov <- as.data.frame(dat$covariates_cen)
fit_pe <- pois_SAP(
  y = dat$ym,
  x1 = dat$xxc[, 1],
  x2 = dat$xxc[, 2],
  efine = dat$ec,
  nlcovfine = cbind(ageing = cov$ageing, unemployed = cov$unemployed),
  C = dat$C_m,
  ndx = ndx,
  ndxnl = c(ndxnl, ndxnl),
  elements = TRUE,
  trace = FALSE
)
stopifnot(!is.null(fit_pe$nleffects), !is.null(fit_pe$sdnleffects))
nl <- as.matrix(fit_pe$nleffects)
sdnl <- as.matrix(fit_pe$sdnleffects)
age_pct <- dat$covariates[, "ageing"] * 100
un_pct <- dat$covariates[, "unemployed"] * 100

ca_age <- smooth_band_df(age_pct, nl[, 1], sdnl[, 1], "ageing")
ca_un <- smooth_band_df(un_pct, nl[, 2], sdnl[, 2], "unemployed")

message("Malone Method 1 refit for smooths + SE...")
source(file.path(root, "R/03_malone.R"))
cov_use <- cov[, c("ageing", "unemployed"), drop = FALSE]
mal_fit <- clgam_malone_fit(
  y_mun = dat$ym,
  e_fine = dat$ec,
  C = dat$C_m,
  cov_fine = cov_use,
  lon = dat$xxc[, 1],
  lat = dat$xxc[, 2],
  cov_smooth = c("ageing", "unemployed"),
  spatial = "none",
  family = "poisson",
  round_counts = TRUE,
  init = "equal",
  deliver = "fitted",
  stop_on = "map",
  max_iter = 200L,
  tol = 1e-3,
  k_cov = ndxnl,
  trace = FALSE
)
tr <- predict(mal_fit$model, type = "terms", se.fit = TRUE)
cn <- colnames(tr$fit)
ix_age <- grep("ageing", cn)[1]
ix_un <- grep("unemployed", cn)[1]
mal_age <- smooth_band_df(
  age_pct, as.numeric(tr$fit[, ix_age]), as.numeric(tr$se.fit[, ix_age]), "ageing"
)
mal_un <- smooth_band_df(
  un_pct, as.numeric(tr$fit[, ix_un]), as.numeric(tr$se.fit[, ix_un]), "unemployed"
)

# Shared y-limits within each covariate (separate panels still comparable)
ylim_age <- range(c(ca_age$lo, ca_age$hi, mal_age$lo, mal_age$hi), finite = TRUE)
ylim_un <- range(c(ca_un$lo, ca_un$hi, mal_un$lo, mal_un$hi), finite = TRUE)
p_ca_age <- plot_smooth_band(ca_age, "CL-GAMM Case A: s(ageing)", "% ageing", cols_method[1], ylim_age)
p_ca_un <- plot_smooth_band(ca_un, "CL-GAMM Case A: s(unemployed)", "% unemployed", cols_method[1], ylim_un)
p_age <- plot_smooth_band(mal_age, "Malone: s(ageing)", "% ageing", cols_method[2], ylim_age)
p_un <- plot_smooth_band(mal_un, "Malone: s(unemployed)", "% unemployed", cols_method[2], ylim_un)

p_ov_age <- plot_overlay_bands(ca_age, mal_age, "s(ageing)", "% ageing")
p_ov_un <- plot_overlay_bands(ca_un, mal_un, "s(unemployed)", "% unemployed")

bw_summ <- data.frame(
  term = c("ageing", "ageing", "unemployed", "unemployed"),
  method = c("CL-GAMM Case A", "Malone", "CL-GAMM Case A", "Malone"),
  median_width = c(
    median(ca_age$hi - ca_age$lo),
    median(mal_age$hi - mal_age$lo),
    median(ca_un$hi - ca_un$lo),
    median(mal_un$hi - mal_un$lo)
  )
)

# save band data for reuse
saveRDS(
  list(
    case_a = list(ageing = ca_age, unemployed = ca_un),
    malone = list(ageing = mal_age, unemployed = mal_un),
    band_width_summary = bw_summ,
    note = paste(
      "Centred partial effects with pointwise 95% bands (fit ± 1.96 SE).",
      "Case A SE from pois_SAP composite-link; Malone SE from mgcv on last-iteration pseudo-counts."
    )
  ),
  file.path(CLGAM_OUTPUT, "madrid_smooth_bands.rds")
)

# MSE summary bar
mse <- data.frame(
  method = c("ATA spatial", "Case A", "Malone"),
  mse = c(
    mean((eta_ata_spatial - logSMR_ct)^2),
    mean((eta_casea - logSMR_ct)^2),
    mean((eta_mal - logSMR_ct)^2)
  )
)
p_mse <- ggplot(mse, aes(reorder(method, mse), mse)) +
  geom_col(fill = "#4C78A8", width = 0.7) +
  coord_flip() +
  labs(
    title = "Oracle MSE of eta vs CT log(SMR)",
    x = NULL, y = "mean (eta - log((yc+1)/ec))^2"
  ) +
  theme_bw(base_size = 11)

# Write multi-page PDF
out_pdf <- file.path(CLGAM_OUTPUT, "madrid_comparison_plots.pdf")
pdf(out_pdf, width = 11, height = 8)
if (requireNamespace("patchwork", quietly = TRUE)) {
  print((p1 + p2) / (p3 + p5) + patchwork::plot_annotation(
    title = "Madrid MEDEA - nested ATA maps (common colour breaks)"
  ))
  print((p4 + p6) / (p7 + p8) + patchwork::plot_annotation(
    title = "ST-PCLM ATA nested / oracle / residuals"
  ))
  print(p_sc)
  print((p_ca_age + p_ca_un) / (p_age + p_un) + patchwork::plot_annotation(
    title = "Separate panels, shared y within covariate (band contrast)",
    subtitle = "Case A SE from composite link; Malone SE from heuristic GAM on pseudo-counts"
  ))
  print(p_ov_age + p_ov_un + patchwork::plot_annotation(
    title = "Partial effects +/- 1.96 SE: CL-GAMM Case A vs Malone"
  ))
  print(p_mse)
} else {
  print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)
  print(p7); print(p8); print(p_sc)
  print(p_ca_age); print(p_ca_un); print(p_age); print(p_un)
  print(p_ov_age); print(p_ov_un); print(p_mse)
}
dev.off()

# Also save a single PNG overview
png(file.path(CLGAM_OUTPUT, "madrid_maps_overview.png"),
    width = 1600, height = 1000, res = 140)
if (requireNamespace("patchwork", quietly = TRUE)) {
  print((p2 + p3) / (p5 + p6) + patchwork::plot_annotation(
    title = "Madrid nested ATA - ATA / Case A / Malone / oracle"
  ))
} else {
  print(p3)
}
dev.off()

# Main partials PNG: overlay with shared y-scale (bands carry the contrast)
png(file.path(CLGAM_OUTPUT, "madrid_partials.png"),
    width = 1500, height = 650, res = 140)
if (requireNamespace("patchwork", quietly = TRUE)) {
  print(p_ov_age + p_ov_un + patchwork::plot_annotation(
    title = "Univariate smooths +/- 1.96 SE: CL-GAMM Case A vs Malone"
  ))
} else {
  print(p_ov_un)
}
dev.off()

png(file.path(CLGAM_OUTPUT, "madrid_partials_separate.png"),
    width = 1400, height = 900, res = 140)
if (requireNamespace("patchwork", quietly = TRUE)) {
  print((p_ca_age + p_ca_un) / (p_age + p_un) + patchwork::plot_annotation(
    title = "Separate panels with shared y-limits within each covariate"
  ))
} else {
  print(p_ca_un)
}
dev.off()

message("Wrote ", out_pdf)
message("Wrote ", file.path(CLGAM_OUTPUT, "madrid_maps_overview.png"))
message("Wrote ", file.path(CLGAM_OUTPUT, "madrid_partials.png"))
message("Wrote ", file.path(CLGAM_OUTPUT, "madrid_partials_separate.png"))
message("Wrote ", file.path(CLGAM_OUTPUT, "madrid_smooth_bands.rds"))
print(bw_summ)
print(mse)

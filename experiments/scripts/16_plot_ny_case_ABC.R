#!/usr/bin/env Rscript
# Plots for NY Case A/B/C (ASCII labels only -- PNG devices drop many Unicode glyphs).
#
# From experiments/:  Rscript scripts/16_plot_ny_case_ABC.R
# Requires: output/ny_case_ABC_fit.rds from scripts/15_ny_case_ABC.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/06_load_ny_leukemia.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(sp)
})

out <- readRDS(file.path(CLGAM_OUTPUT, "ny_case_ABC_fit.rds"))
dat <- clgam_load_ny_leukemia(n_coarse = "county")
stopifnot(identical(out$data$coarse_scheme, "county_FIPS"))

data(NY8, package = "DClusterm")
ny_sf <- st_as_sf(NY8)
ny_sf$AREAKEY <- as.character(ny_sf$AREAKEY)
ny_sf <- ny_sf[match(dat$fips, ny_sf$AREAKEY), ]
stopifnot(identical(dat$fips, ny_sf$AREAKEY))

log_smr <- log((dat$cases_tract + 0.5) / dat$e_tract)
ny_sf$eta_A <- out$case_a$eta
ny_sf$eta_B <- out$case_b$eta
ny_sf$eta_C <- out$case_c$eta
ny_sf$oracle <- log_smr

ny8df <- as.data.frame(NY8)
ny8df$AREAKEY <- as.character(ny8df$AREAKEY)
ord <- match(dat$fips, ny8df$AREAKEY)
ny_sf$ownhome_pct <- as.numeric(ny8df$PCTOWNHOME[ord]) * 100
pexp <- as.numeric(ny8df$PEXPOSURE[ord])
pexp_c <- as.numeric(tapply(pexp, factor(dat$group, levels = dat$group_levels), mean))
ny_sf$tce_county <- pexp_c[match(dat$group, dat$group_levels)]

brks <- quantile(
  c(ny_sf$eta_A, ny_sf$eta_B, ny_sf$eta_C, ny_sf$oracle),
  probs = seq(0, 1, 0.1), na.rm = TRUE
)
brks[1] <- brks[1] - 1e-6
brks[length(brks)] <- brks[length(brks)] + 1e-6
pal <- scales::viridis_pal(option = "C")(10)

map1 <- function(var, title) {
  d <- ny_sf
  d$bin <- cut(d[[var]], breaks = brks, include.lowest = TRUE)
  ggplot(d) +
    geom_sf(aes(fill = bin), colour = NA) +
    scale_fill_manual(values = pal, drop = FALSE, name = "log-rate") +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 12)
    )
}

p1 <- map1("oracle", "Oracle tract log((cases+0.5)/e)")
p2 <- map1("eta_A", "Case A: spatial + ownhome")
p3 <- map1("eta_B", "Case B: spatial + tce_county")
p4 <- map1("eta_C", "Case C: ownhome + tce_county")

# In-sample check at the scale of the likelihood (8 counties), not tract oracle.
# Tract scatter looks flat with n=8 + poor bases (heavy shrinkage) -- misleading.
county_fit_df <- function(eta, method) {
  mu <- as.numeric(dat$C %*% (dat$e_tract * exp(as.numeric(eta))))
  data.frame(
    observed = as.numeric(dat$y),
    fitted = mu,
    county = dat$group_levels,
    method = method,
    stringsAsFactors = FALSE
  )
}
sc_c <- rbind(
  county_fit_df(out$case_a$eta, "A"),
  county_fit_df(out$case_b$eta, "B"),
  county_fit_df(out$case_c$eta, "C")
)
p_sc <- ggplot(sc_c, aes(observed, fitted)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 3, colour = "#2166ac") +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~method) +
  labs(
    title = "County fit: observed y vs C %*% (e*exp(eta)) [log-log]",
    subtitle = "In-sample at likelihood scale (8 counties)",
    x = "observed county cases (log10)", y = "fitted county mean (log10)"
  ) +
  theme_bw(base_size = 11) +
  coord_equal()

# Tract oracle kept as secondary panel (expect weak match with n=8)
sc_tr <- rbind(
  data.frame(oracle = log_smr, fit = out$case_a$eta, method = "A"),
  data.frame(oracle = log_smr, fit = out$case_b$eta, method = "B"),
  data.frame(oracle = log_smr, fit = out$case_c$eta, method = "C")
)
p_tr <- ggplot(sc_tr, aes(oracle, fit)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.25, size = 0.7) +
  facet_wrap(~method) +
  labs(
    title = "Tract eta vs oracle (descriptive only)",
    subtitle = "Flat cloud expected: only 8 county totals inform the fine surface",
    x = "tract oracle log((cases+0.5)/e)", y = "fitted eta"
  ) +
  theme_bw(base_size = 11)

# Linear effects from saved fits (recompute leffects quickly)
own <- dat$covariates_fine[, "ownhome", drop = FALSE]
tce <- dat$covariates_agg_expanded
fit_a <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_tract,
  lcovfine = own, C = dat$C, ndx = c(3, 3), elements = TRUE, trace = FALSE
)
fit_b <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_tract,
  lcovfine = tce, C = dat$C, ndx = c(3, 3), elements = TRUE, trace = FALSE
)
fit_c <- pois_SAP(
  y = dat$y, x1 = dat$x1, x2 = dat$x2, efine = dat$e_tract,
  lcovfine = cbind(own, tce), C = dat$C, ndx = c(3, 3), elements = TRUE, trace = FALSE
)

lin_band_df <- function(x, fit, se, method) {
  fit <- as.numeric(fit)
  se <- as.numeric(se)
  ctr <- mean(fit, na.rm = TRUE)
  data.frame(
    x = x,
    fit = fit - ctr,
    lo = fit - 1.96 * se - ctr,
    hi = fit + 1.96 * se - ctr,
    method = method,
    stringsAsFactors = FALSE
  )
}

agg_lin <- function(df) {
  aggregate(cbind(fit, lo, hi) ~ x + method, data = df, FUN = mean)
}

plot_lin_band <- function(df, title, xlab) {
  ag <- agg_lin(df)
  ggplot(ag, aes(x = x, colour = method, fill = method)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, colour = NA) +
    geom_line(aes(y = fit), linewidth = 0.95) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = 3, linewidth = 0.4) +
    labs(
      title = title,
      subtitle = "centred partial +/- 1.96 SE (pois_SAP)",
      x = xlab, y = "partial (centred)", colour = NULL, fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

pa <- lin_band_df(
  ny_sf$ownhome_pct, fit_a$leffects, fit_a$sdleffects, "A: ownhome"
)
pb <- lin_band_df(
  ny_sf$tce_county, fit_b$leffects, fit_b$sdleffects, "B: tce_county"
)
pc1 <- lin_band_df(
  ny_sf$ownhome_pct, fit_c$leffects[, 1], fit_c$sdleffects[, 1], "C: ownhome"
)
pc2 <- lin_band_df(
  ny_sf$tce_county, fit_c$leffects[, 2], fit_c$sdleffects[, 2], "C: tce_county"
)

p_own <- plot_lin_band(rbind(pa, pc1), "Linear effect: % own home", "% own home")
p_tce <- plot_lin_band(
  rbind(pb, pc2), "Linear effect: county TCE exposure", "county mean PEXPOSURE"
)

p_aic <- ggplot(out$table, aes(case, aic, fill = case)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", aic)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(A = "#2166ac", B = "#b2182b", C = "#4d9221")) +
  labs(title = "AIC (8 counties, ndx=3,3, linear)", y = "AIC", x = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(min(out$table$aic) - 1, max(out$table$aic) + 1))

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Need patchwork for multi-panel PNGs")
}

png(file.path(CLGAM_OUTPUT, "ny_case_ABC_maps.png"), width = 1600, height = 1100, res = 140)
print((p1 + p2) / (p3 + p4) + patchwork::plot_annotation(
  title = "NY leukemia - nested ATA (8 counties -> 281 tracts)",
  subtitle = "Poor bases ndx=(3,3); linear Case A/B/C"
))
dev.off()

png(file.path(CLGAM_OUTPUT, "ny_case_ABC_diagnostics.png"), width = 1500, height = 1200, res = 140)
print(
  (p_sc) / (p_own + p_tce) / (p_tr) +
    patchwork::plot_annotation(
      title = "NY Case A/B/C diagnostics",
      subtitle = "Top: county in-sample fit | Middle: linear partials +/- 1.96 SE | Bottom: tract oracle (weak by design)"
    ) +
    patchwork::plot_layout(heights = c(1.1, 1, 0.9))
)
dev.off()

png(file.path(CLGAM_OUTPUT, "ny_case_ABC_aic.png"), width = 700, height = 500, res = 140)
print(p_aic)
dev.off()

message("Wrote ", file.path(CLGAM_OUTPUT, "ny_case_ABC_maps.png"))
message("Wrote ", file.path(CLGAM_OUTPUT, "ny_case_ABC_diagnostics.png"))
message("Wrote ", file.path(CLGAM_OUTPUT, "ny_case_ABC_aic.png"))

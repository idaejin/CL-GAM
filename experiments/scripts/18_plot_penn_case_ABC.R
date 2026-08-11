#!/usr/bin/env Rscript
# Plots for pennLC Case A/B/C (ASCII labels).
# From experiments/:  Rscript scripts/18_plot_penn_case_ABC.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/07_load_pennLC.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(sp)
})

out <- readRDS(file.path(CLGAM_OUTPUT, "penn_case_ABC_fit.rds"))
dat <- clgam_load_pennLC(n_coarse = 20L, seed = 1L)

library(SpatialEpi)
data(pennLC, package = "SpatialEpi")
pa_sf <- st_as_sf(pennLC$spatial.polygon)
# spatial.polygon order may differ; match by county name from row.names / IDs
# pennLC spatial polygons often named by county
ids <- tolower(gsub("^.*\\.", "", row.names(pennLC$spatial.polygon)))
# fallback: match geo order
if (!all(dat$county %in% ids)) {
  # try polygon ID slot
  ids <- tolower(sapply(pennLC$spatial.polygon@polygons, function(p) p@ID))
}
ord <- match(dat$county, ids)
if (anyNA(ord)) {
  # build sf from geo only for point maps if polygon match fails
  pa_sf <- st_as_sf(
    data.frame(county = dat$county, x = dat$x1, y = dat$x2),
    coords = c("x", "y"), crs = 4326
  )
  use_poly <- FALSE
} else {
  pa_sf <- pa_sf[ord, ]
  pa_sf$county <- dat$county
  use_poly <- TRUE
}

log_smr <- log((dat$y_fine + 0.5) / dat$e_fine)
pa_sf$eta_A <- out$case_a$eta
pa_sf$eta_B <- out$case_b$eta
pa_sf$eta_C <- out$case_c$eta
pa_sf$oracle <- log_smr
pa_sf$smoking_pct <- dat$smoking * 100

brks <- quantile(
  c(pa_sf$eta_A, pa_sf$eta_B, pa_sf$eta_C, pa_sf$oracle),
  probs = seq(0, 1, 0.1), na.rm = TRUE
)
brks[1] <- brks[1] - 1e-6
brks[length(brks)] <- brks[length(brks)] + 1e-6
pal <- scales::viridis_pal(option = "C")(10)

map1 <- function(var, title) {
  d <- pa_sf
  d$bin <- cut(d[[var]], breaks = brks, include.lowest = TRUE)
  p <- ggplot(d)
  if (use_poly) {
    p <- p + geom_sf(aes(fill = bin), colour = "grey40", linewidth = 0.1)
  } else {
    p <- p + geom_sf(aes(colour = bin), size = 2.5) +
      scale_colour_manual(values = pal, drop = FALSE, name = "log-rate")
  }
  if (use_poly) {
    p <- p + scale_fill_manual(values = pal, drop = FALSE, name = "log-rate")
  }
  p + labs(title = title) +
    theme_void(base_size = 11) +
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 11))
}

p1 <- map1("oracle", "Oracle county log((Y+0.5)/E)")
p2 <- map1("eta_A", "Case A: spatial + s(smoking)")
p3 <- map1("eta_B", "Case B: spatial + north_agg")
p4 <- map1("eta_C", "Case C: smoking + north_agg")

# Coarse in-sample
county_fit_df <- function(eta, method) {
  mu <- as.numeric(dat$C %*% (dat$e_fine * exp(as.numeric(eta))))
  data.frame(
    observed = as.numeric(dat$y), fitted = mu, method = method,
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
  geom_point(size = 2.5, colour = "#2166ac") +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~method) +
  labs(
    title = "Region fit: observed y vs C %*% (e*exp(eta)) [log-log]",
    subtitle = sprintf("In-sample (%d k-means regions)", dat$n_coarse),
    x = "observed region cases (log10)", y = "fitted region mean (log10)"
  ) +
  theme_bw(base_size = 11) +
  coord_equal()

# Fine oracle scatter
sc_f <- rbind(
  data.frame(oracle = log_smr, fit = out$case_a$eta, method = "A"),
  data.frame(oracle = log_smr, fit = out$case_b$eta, method = "B"),
  data.frame(oracle = log_smr, fit = out$case_c$eta, method = "C")
)
p_f <- ggplot(sc_f, aes(oracle, fit)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.55, size = 1.6) +
  facet_wrap(~method) +
  labs(
    title = "County eta vs oracle",
    subtitle = sprintf(
      "cor A=%.2f B=%.2f C=%.2f",
      out$table$cor_vs_county[1], out$table$cor_vs_county[2], out$table$cor_vs_county[3]
    ),
    x = "county oracle log((Y+0.5)/E)", y = "fitted eta"
  ) +
  theme_bw(base_size = 11)

# Smoking smooth / bands from Case A and C
lin_band <- function(x, fit, se, method) {
  fit <- as.numeric(fit); se <- as.numeric(se)
  ctr <- mean(fit)
  o <- order(x)
  data.frame(
    x = x[o], fit = fit[o] - ctr,
    lo = fit[o] - 1.96 * se[o] - ctr,
    hi = fit[o] + 1.96 * se[o] - ctr,
    method = method
  )
}
# Case A: nleffects smoking; Case C: nleffects smoking + leffects north
pa <- lin_band(
  pa_sf$smoking_pct, out$case_a$nleffects, out$case_a$sdnleffects, "A: s(smoking)"
)
pc <- lin_band(
  pa_sf$smoking_pct, out$case_c$nleffects, out$case_c$sdnleffects, "C: s(smoking)"
)
ag <- aggregate(cbind(fit, lo, hi) ~ x + method, rbind(pa, pc), mean)
p_sm <- ggplot(ag, aes(x, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, colour = NA) +
  geom_line(aes(y = fit), linewidth = 0.95) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  labs(
    title = "s(smoking) +/- 1.96 SE",
    x = "% smoking", y = "partial (centred)", colour = NULL, fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_aic <- ggplot(out$table, aes(case, aic, fill = case)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", aic)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(A = "#2166ac", B = "#b2182b", C = "#4d9221")) +
  labs(title = "AIC (pennLC, 20 regions -> 67 counties)", y = "AIC", x = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

stopifnot(requireNamespace("patchwork", quietly = TRUE))

png(file.path(CLGAM_OUTPUT, "penn_case_ABC_maps.png"), width = 1600, height = 1100, res = 140)
print((p1 + p2) / (p3 + p4) + patchwork::plot_annotation(
  title = "PA lung cancer (pennLC) - nested ATA (20 regions -> 67 counties)",
  subtitle = "SpatialEpi::pennLC; Case A/B/C"
))
dev.off()

png(file.path(CLGAM_OUTPUT, "penn_case_ABC_diagnostics.png"), width = 1500, height = 1200, res = 140)
print(
  (p_sc) / (p_f + p_sm) / (p_aic) +
    patchwork::plot_annotation(
      title = "pennLC Case A/B/C diagnostics",
      subtitle = "Top: region in-sample | Middle: county oracle + smoking smooth | Bottom: AIC"
    ) +
    patchwork::plot_layout(heights = c(1, 1.1, 0.7))
)
dev.off()

message("Wrote ", file.path(CLGAM_OUTPUT, "penn_case_ABC_maps.png"))
message("Wrote ", file.path(CLGAM_OUTPUT, "penn_case_ABC_diagnostics.png"))
print(out$table)

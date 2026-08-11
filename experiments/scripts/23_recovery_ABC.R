# Shared geometry + scales: Case A, B, C recovery figures (2×4)
# Writes:
#   experiments/output/clgam_recovery_caseA.png
#   experiments/output/clgam_recovery_caseB.png
#   experiments/output/clgam_recovery_caseC.png
#
# From project root: Rscript experiments/scripts/23_recovery_ABC.R

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
  library(sf)
})
suppressMessages(sf_use_s2(FALSE))

out_dir <- "experiments/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Same Voronoi for A/B/C (seed + sizes fixed; geometry drawn before covariates)
SEED <- 6L
n_coarse <- 40L
n_fine_per <- 10L
spatial_amp_ab <- 0.8
nl_amp_ab <- 1.2
nl_fun_sin <- function(z) sin(2 * pi * z)
# Case C (easy demo): two distinct smooths that remain identifiable under ATA.
# Full-period fine sine is washed out by aggregation; use descending + ascending.
nl_fun_desc <- function(z) -tanh(2.5 * (z - 0.5))
nl_fun_asc <- function(z) tanh(2.5 * (z - 0.5))
spatial_amp_c <- 0.75
nl_amp_c <- c(1.2, 1.2)

# ---- simulate + fit ----------------------------------------------------------
dat_A <- simulate_ata(
  n_coarse = n_coarse, n_fine_per = n_fine_per, seed = SEED,
  include_covariate = TRUE, covariate_level = "fine",
  spatial_amp = spatial_amp_ab, nl_amp = nl_amp_ab, nl_fun = nl_fun_sin
)
fit_A <- clgam(
  dat_A$y, dat_A$x1, dat_A$x2, dat_A$C,
  exposure = dat_A$efine, smooth = dat_A$nlcovfine,
  knots = c(12, 12), knots_nl = 12, elements = TRUE, trace = FALSE
)

dat_B <- simulate_ata(
  n_coarse = n_coarse, n_fine_per = n_fine_per, seed = SEED,
  include_covariate = TRUE, covariate_level = "coarse",
  spatial_amp = spatial_amp_ab, nl_amp = nl_amp_ab, nl_fun = nl_fun_sin
)
fit_B <- clgam(
  dat_B$y, dat_B$x1, dat_B$x2, dat_B$C,
  exposure = dat_B$efine, smooth = dat_B$nlcovfine,
  knots = c(10, 10), knots_nl = 6, elements = TRUE, trace = FALSE
)

dat_C <- simulate_ata(
  n_coarse = n_coarse, n_fine_per = n_fine_per, seed = SEED,
  include_covariate = TRUE, covariate_level = "both",
  spatial_amp = spatial_amp_c, nl_amp = nl_amp_c,
  nl_fun = list(nl_fun_desc, nl_fun_asc)
)
fit_C <- clgam(
  dat_C$y, dat_C$x1, dat_C$x2, dat_C$C,
  exposure = dat_C$efine, smooth = dat_C$nlcovfine,
  knots = c(10, 10), knots_nl = c(8, 6), elements = TRUE, trace = FALSE
)

# Geometry check
stopifnot(all.equal(
  as.numeric(sf::st_area(dat_A$sf_coarse)),
  as.numeric(sf::st_area(dat_B$sf_coarse)),
  tolerance = 1e-12
))
stopifnot(all.equal(
  as.numeric(sf::st_area(dat_A$sf_fine)),
  as.numeric(sf::st_area(dat_C$sf_fine)),
  tolerance = 1e-12
))
cat("OK: A/B/C share the same Voronoi geometry (seed=", SEED, ")\n", sep = "")

# ---- helpers -----------------------------------------------------------------
prep1 <- function(dat, fit) {
  nk <- ncol(fit$nleffects)
  nl_hat <- fit$nleffects - matrix(colMeans(fit$nleffects), nrow(fit$nleffects), nk, byrow = TRUE)
  f_hat <- fit$eta - rowSums(fit$nleffects)
  f_hat <- f_hat - mean(f_hat)
  f_true <- dat$eta_spatial_true - mean(dat$eta_spatial_true)
  list(
    dat = dat, fit = fit, case = dat$case, nk = nk,
    nl_hat = nl_hat, se = fit$sdnleffects,
    f_hat = f_hat, f_true = f_true,
    eta_hat = fit$eta - mean(fit$eta),
    eta_true = dat$eta_true - mean(dat$eta_true),
    mse_eta = mean((fit$eta - dat$eta_true)^2)
  )
}
A <- prep1(dat_A, fit_A)
B <- prep1(dat_B, fit_B)
C <- prep1(dat_C, fit_C)

# Single-smooth curve helpers (A/B)
curve_one <- function(obj, kn, amp, fun, coarse = FALSE) {
  dat <- obj$dat
  fit <- obj$fit
  z_obs <- fit$nlcovfine[, 1]
  if (coarse) {
    idx1 <- match(seq_len(dat$n_coarse), dat$coarse_id)
    za <- dat$z_a
    oa <- order(za)
    h_hat_a <- obj$nl_hat[idx1, 1]
    se_a <- obj$se[idx1, 1]
    z_lo <- min(za); z_hi <- max(za)
    z_grid <- seq(z_lo, z_hi, length.out = 300)
    sf_hat <- stats::splinefun(za[oa], h_hat_a[oa], method = "natural")
    sf_se <- stats::splinefun(za[oa], se_a[oa], method = "natural")
    list(
      z = z_grid, hat = sf_hat(z_grid), se = pmax(sf_se(z_grid), 0),
      true = amp * fun(z_grid) - mean(amp * fun(za)),
      pts_z = za, pts_true = dat$h_true_coarse, pts_hat = h_hat_a,
      cor_nl = cor(h_hat_a, dat$h_true_coarse)
    )
  } else {
    oa <- order(z_obs)
    z_grid <- seq(min(z_obs), max(z_obs), length.out = 300)
    hat <- approx(z_obs[oa], obj$nl_hat[oa, 1], xout = z_grid, rule = 2)$y
    se <- approx(z_obs[oa], obj$se[oa, 1], xout = z_grid, rule = 2)$y
    zt <- if (!is.null(dat$z_f)) dat$z_f else dat$z
    list(
      z = z_grid, hat = hat, se = pmax(se, 0),
      true = amp * fun(z_grid) - mean(amp * fun(zt)),
      pts_z = NULL, pts_true = NULL, pts_hat = NULL,
      cor_nl = cor(obj$nl_hat[, 1], dat$g_true - mean(dat$g_true))
    )
  }
}

cur_A <- curve_one(A, 12L, nl_amp_ab, nl_fun_sin, coarse = FALSE)
cur_B <- curve_one(B, 6L, nl_amp_ab, nl_fun_sin, coarse = TRUE)

# Case C dual curves
amp_g <- nl_amp_c[1]; amp_h <- nl_amp_c[2]
z_grid <- seq(0, 1, length.out = 300)
g_true_c <- amp_g * nl_fun_desc(z_grid) - mean(amp_g * nl_fun_desc(dat_C$z_f))
og <- order(dat_C$z_f)
g_hat_c <- approx(dat_C$z_f[og], C$nl_hat[og, 1], xout = z_grid, rule = 2)$y
se_g <- approx(dat_C$z_f[og], C$se[og, 1], xout = z_grid, rule = 2)$y
idx1 <- match(seq_len(dat_C$n_coarse), dat_C$coarse_id)
za <- dat_C$z_a
oa <- order(za)
z_h <- seq(min(za), max(za), length.out = 300)
h_hat_a <- C$nl_hat[idx1, 2]
h_true_a <- dat_C$h_true_coarse
h_hat_c <- stats::splinefun(za[oa], h_hat_a[oa], method = "natural")(z_h)
h_true_h <- amp_h * nl_fun_asc(z_h) - mean(amp_h * nl_fun_asc(za))
se_h <- stats::splinefun(za[oa], C$se[idx1, 2][oa], method = "natural")(z_h)
cor_g_C <- cor(C$nl_hat[, 1], dat_C$g_true - mean(dat_C$g_true))
cor_h_C <- cor(h_hat_a, h_true_a)

cat(sprintf("A: cor_f=%.3f cor_g=%.3f MSE=%.3f\n",
            cor(A$f_hat, A$f_true), cur_A$cor_nl, A$mse_eta))
cat(sprintf("B: cor_f=%.3f cor_h=%.3f MSE=%.3f\n",
            cor(B$f_hat, B$f_true), cur_B$cor_nl, B$mse_eta))
cat(sprintf("C: cor_f=%.3f cor_g=%.3f cor_h=%.3f MSE=%.3f\n",
            cor(C$f_hat, C$f_true), cor_g_C, cor_h_C, C$mse_eta))

# ---- shared scales -----------------------------------------------------------
pal_div <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))
pal_seq <- grDevices::colorRampPalette(c("#fff7bc", "#fec44f", "#d95f0e", "#7f0000"))
n_pal <- 101L

lim_f <- {
  m <- max(abs(c(A$f_true, A$f_hat, B$f_true, B$f_hat, C$f_true, C$f_hat)))
  c(-m, m)
}
# counts: panels 3–4 share one scale (fine + coarse) within each figure
lim_mu <- range(c(dat_A$y, fit_A$mu, dat_B$y, fit_B$mu, dat_C$y, fit_C$mu))
lim_eta <- {
  m <- max(abs(c(A$eta_true, A$eta_hat, B$eta_true, B$eta_hat, C$eta_true, C$eta_hat)))
  c(-m, m)
}
lim_nl <- range(c(
  cur_A$true, cur_A$hat, cur_B$true, cur_B$hat,
  g_true_c, g_hat_c, h_true_h, h_hat_c, h_true_a, h_hat_a
), finite = TRUE)
lim_nl <- lim_nl + c(-1, 1) * 0.08 * diff(lim_nl)
xlim_z <- c(0, 1)

cols_from <- function(v, lim, pal, n = n_pal) {
  vv <- pmin(pmax(v, lim[1]), lim[2])
  br <- seq(lim[1], lim[2], length.out = n)
  pal(n)[cut(vv, breaks = br, include.lowest = TRUE)]
}

# Vertical colour bar with numeric ticks (panels 1/5 and 3/4)
add_colorbar <- function(lim, pal, n = n_pal, digits = NULL) {
  usr <- par("usr")
  dx <- diff(usr[1:2])
  x0 <- usr[2] + 0.03 * dx
  x1 <- usr[2] + 0.08 * dx
  cols <- pal(n)
  ys <- seq(usr[3], usr[4], length.out = n + 1L)
  for (i in seq_len(n)) {
    rect(x0, ys[i], x1, ys[i + 1L], col = cols[i], border = NA, xpd = NA)
  }
  rect(x0, usr[3], x1, usr[4], border = "grey30", xpd = NA)
  labs <- pretty(lim, n = 5)
  labs <- labs[labs >= lim[1] - 1e-9 & labs <= lim[2] + 1e-9]
  ylab <- usr[3] + (labs - lim[1]) / diff(lim) * diff(usr[3:4])
  if (is.null(digits)) {
    lab_txt <- format(labs, trim = TRUE, scientific = FALSE)
  } else {
    lab_txt <- format(round(labs, digits), nsmall = digits, trim = TRUE)
  }
  text(x1 + 0.015 * dx, ylab, labels = lab_txt, adj = 0, cex = 0.62, xpd = NA)
}

mar_map <- c(3.6, 3.2, 2.8, 3.6)
mar_other <- c(3.6, 3.6, 2.8, 1.0)

plot_chor <- function(sf_f, sf_c, v, lim, pal, main) {
  par(mar = mar_map)
  bb <- sf::st_bbox(sf_f)
  plot(NA, type = "n",
       xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
       xlab = "", ylab = "", main = main, xaxs = "i", yaxs = "i", asp = NA)
  plot(sf::st_geometry(sf_f), col = cols_from(v, lim, pal), border = NA, add = TRUE)
  plot(sf::st_geometry(sf_c), border = "grey20", lwd = 0.55, add = TRUE)
  box(col = "grey35", lwd = 0.9)
  add_colorbar(lim, pal, digits = 2)
}

plot_smooth_ab <- function(cur, ylab, main, show_hat = TRUE, is_B = FALSE) {
  par(mar = mar_other)
  plot(cur$z, cur$true, type = "n", xlim = xlim_z, ylim = lim_nl,
       xlab = if (is_B) expression(z[a] ~ "(coarse)") else expression(z ~ "(fine)"),
       ylab = ylab, main = main, xaxs = "i", yaxs = "i")
  if (isTRUE(show_hat)) {
    lo <- pmax(cur$hat - 1.96 * cur$se, lim_nl[1])
    hi <- pmin(cur$hat + 1.96 * cur$se, lim_nl[2])
    polygon(c(cur$z, rev(cur$z)), c(lo, rev(hi)),
            col = adjustcolor("#b2182b", 0.18), border = NA)
  }
  lines(cur$z, cur$true, lwd = 2.2, col = "grey20")
  if (isTRUE(show_hat)) lines(cur$z, cur$hat, lwd = 2.4, col = "#b2182b")
  if (is_B && !is.null(cur$pts_z)) {
    oa <- order(cur$pts_z)
    points(cur$pts_z[oa], cur$pts_true[oa], pch = 1, cex = 0.6, col = "grey20")
    if (isTRUE(show_hat))
      points(cur$pts_z[oa], cur$pts_hat[oa], pch = 16, cex = 0.5, col = "#b2182b")
  }
  if (isTRUE(show_hat)) {
    legend("topright",
           legend = if (is_B) c(expression(h[true]), expression(hat(h)))
           else c(expression(g[true]), expression(hat(g)), "95%"),
           col = c("grey20", "#b2182b", adjustcolor("#b2182b", 0.35)),
           lwd = c(2.2, 2.4, if (is_B) NA else 8)[seq_len(2 + !is_B)],
           bty = "n", cex = 0.72)
  } else {
    legend("topright",
           legend = if (is_B) expression(h[true](z)) else expression(g[true](z)),
           col = "grey20", lwd = 2.2, bty = "n", cex = 0.72)
  }
  box(col = "grey35", lwd = 0.8)
}

plot_counts <- function(dat, lim_counts) {
  sf_f <- dat$sf_fine; sf_c <- dat$sf_coarse
  par(mar = mar_map)
  bb <- sf::st_bbox(sf_f)
  plot(NA, type = "n", xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
       xlab = "", ylab = "", xaxs = "i", yaxs = "i", asp = NA,
       main = expression(paste("3. Intermediate  ", y[f], " ~ Poisson(", italic(f) + italic(s), ")")))
  plot(sf::st_geometry(sf_f), col = cols_from(dat$y_fine, lim_counts, pal_seq),
       border = NA, add = TRUE)
  plot(sf::st_geometry(sf_c), border = "grey20", lwd = 0.55, add = TRUE)
  box(col = "grey35", lwd = 0.9)
  add_colorbar(lim_counts, pal_seq, digits = 0)

  par(mar = mar_map)
  bb_c <- sf::st_bbox(sf_c)
  plot(NA, type = "n", xlim = c(bb_c["xmin"], bb_c["xmax"]), ylim = c(bb_c["ymin"], bb_c["ymax"]),
       xlab = "", ylab = "", xaxs = "i", yaxs = "i", asp = NA,
       main = expression(paste("4. Observed  ", y == C * y[f])))
  plot(sf::st_geometry(sf_c), col = cols_from(dat$y, lim_counts, pal_seq),
       border = "grey25", lwd = 0.8, add = TRUE)
  box(col = "grey35", lwd = 0.9)
  add_colorbar(lim_counts, pal_seq, digits = 0)
}

plot_fit_eta <- function(obj) {
  par(mar = mar_other)
  dat <- obj$dat; fit <- obj$fit
  plot(dat$y, fit$mu, pch = 16, cex = 0.85, col = adjustcolor("darkgreen", 0.75),
       xlim = lim_mu, ylim = lim_mu, xaxs = "i", yaxs = "i",
       xlab = "y (observed coarse)", ylab = expression(hat(mu)),
       main = sprintf("7. Aggregate fit  cor=%.2f", cor(fit$mu, dat$y)))
  abline(0, 1, lty = 2, col = "grey40")
  box(col = "grey35", lwd = 0.8)

  plot(obj$eta_true, obj$eta_hat, pch = 16, cex = 0.45, col = adjustcolor("steelblue", 0.5),
       xlim = lim_eta, ylim = lim_eta, xaxs = "i", yaxs = "i",
       xlab = expression(eta[true]), ylab = expression(hat(eta)),
       main = sprintf("8. eta vs truth  MSE=%.3f  cor=%.2f",
                      obj$mse_eta, cor(obj$eta_hat, obj$eta_true)))
  abline(0, 1, lty = 2, col = "grey40")
  box(col = "grey35", lwd = 0.8)
}

start_fig <- function(file) {
  png(file, width = 1900, height = 900, res = 140)
  layout(matrix(1:8, nrow = 2, byrow = TRUE), widths = rep(1, 4), heights = rep(1, 2))
  par(oma = c(0, 0, 0, 0), mgp = c(2.0, 0.6, 0), tcl = -0.3, font.main = 2)
}

# ---- Case A ------------------------------------------------------------------
fA <- file.path(out_dir, "clgam_recovery_caseA.png")
start_fig(fA)
lim_counts_A <- range(c(dat_A$y_fine, dat_A$y), finite = TRUE)
plot_chor(dat_A$sf_fine, dat_A$sf_coarse, A$f_true, lim_f, pal_div, "1. True spatial f")
plot_smooth_ab(cur_A, "g(z)", expression(paste("2. True smooth  ", g(z))), FALSE, FALSE)
plot_counts(dat_A, lim_counts_A)
plot_chor(dat_A$sf_fine, dat_A$sf_coarse, A$f_hat, lim_f, pal_div,
          sprintf("5. Spatial f hat  cor=%.2f", cor(A$f_hat, A$f_true)))
plot_smooth_ab(cur_A, "g(z)", sprintf("6. Additive smooth vs true  cor=%.2f", cur_A$cor_nl), TRUE, FALSE)
plot_fit_eta(A)
dev.off()

# ---- Case B ------------------------------------------------------------------
fB <- file.path(out_dir, "clgam_recovery_caseB.png")
start_fig(fB)
lim_counts_B <- range(c(dat_B$y_fine, dat_B$y), finite = TRUE)
plot_chor(dat_B$sf_fine, dat_B$sf_coarse, B$f_true, lim_f, pal_div, "1. True spatial f")
plot_smooth_ab(cur_B, "h(z)", expression(paste("2. True smooth  ", h(z[a]))), FALSE, TRUE)
plot_counts(dat_B, lim_counts_B)
plot_chor(dat_B$sf_fine, dat_B$sf_coarse, B$f_hat, lim_f, pal_div,
          sprintf("5. Spatial f hat  cor=%.2f", cor(B$f_hat, B$f_true)))
plot_smooth_ab(cur_B, "h(z)", sprintf("6. Additive smooth vs true  cor=%.2f", cur_B$cor_nl), TRUE, TRUE)
plot_fit_eta(B)
dev.off()

# ---- Case C ------------------------------------------------------------------
fC <- file.path(out_dir, "clgam_recovery_caseC.png")
start_fig(fC)
lim_counts_C <- range(c(dat_C$y_fine, dat_C$y), finite = TRUE)
plot_chor(dat_C$sf_fine, dat_C$sf_coarse, C$f_true, lim_f, pal_div, "1. True spatial f")

par(mar = mar_other)
plot(z_grid, g_true_c, type = "n", xlim = xlim_z, ylim = lim_nl,
     xlab = "z", ylab = "effect", xaxs = "i", yaxs = "i",
     main = expression(paste("2. True smooths  ", g == "desc", " & ", h == "asc")))
lines(z_grid, g_true_c, lwd = 2.2, col = "grey20")
lines(z_h, h_true_h, lwd = 2.2, lty = 2, col = "grey20")
points(za[oa], h_true_a[oa], pch = 1, cex = 0.55, col = "grey40")
legend("topleft",
       legend = c(expression(g[true](z[f]) ~ "(desc)"),
                  expression(h[true](z[a]) ~ "(asc)"),
                  expression(z[a])),
       col = "grey20", lwd = c(2.2, 2.2, NA), lty = c(1, 2, NA),
       pch = c(NA, NA, 1), pt.cex = 0.55, bty = "n", cex = 0.68)
box(col = "grey35", lwd = 0.9)

plot_counts(dat_C, lim_counts_C)
plot_chor(dat_C$sf_fine, dat_C$sf_coarse, C$f_hat, lim_f, pal_div,
          sprintf("5. Spatial f hat  cor=%.2f", cor(C$f_hat, C$f_true)))

par(mar = mar_other)
plot(z_grid, g_true_c, type = "n", xlim = xlim_z, ylim = lim_nl,
     xlab = "z", ylab = "effect", xaxs = "i", yaxs = "i",
     main = sprintf("6. Smooths vs true  cor(g,h)=%.2f,%.2f", cor_g_C, cor_h_C))
polygon(c(z_grid, rev(z_grid)),
        c(pmax(g_hat_c - 1.96 * se_g, lim_nl[1]),
          rev(pmin(g_hat_c + 1.96 * se_g, lim_nl[2]))),
        col = adjustcolor("#b2182b", 0.15), border = NA)
polygon(c(z_h, rev(z_h)),
        c(pmax(h_hat_c - 1.96 * se_h, lim_nl[1]),
          rev(pmin(h_hat_c + 1.96 * se_h, lim_nl[2]))),
        col = adjustcolor("#2166ac", 0.15), border = NA)
lines(z_grid, g_true_c, lwd = 1.8, col = "grey30")
lines(z_grid, g_hat_c, lwd = 2.2, col = "#b2182b")
lines(z_h, h_true_h, lwd = 1.8, lty = 2, col = "grey30")
lines(z_h, h_hat_c, lwd = 2.2, lty = 2, col = "#2166ac")
points(za[oa], h_true_a[oa], pch = 1, cex = 0.5, col = "grey40")
points(za[oa], h_hat_a[oa], pch = 16, cex = 0.45, col = "#2166ac")
legend("topleft",
       legend = c(expression(g), expression(hat(g)), expression(h), expression(hat(h))),
       col = c("grey30", "#b2182b", "grey30", "#2166ac"),
       lwd = 2, lty = c(1, 1, 2, 2), bty = "n", cex = 0.68)
box(col = "grey35", lwd = 0.9)

plot_fit_eta(C)
dev.off()

cat("Shared scales:\n")
cat(sprintf("  lim_f=[%.3f,%.3f] (panels 1 & 5)\n", lim_f[1], lim_f[2]))
cat(sprintf("  lim_counts A/B/C = [%s] / [%s] / [%s] (panels 3 & 4)\n",
            paste(round(lim_counts_A), collapse = ","),
            paste(round(lim_counts_B), collapse = ","),
            paste(round(lim_counts_C), collapse = ",")))
cat("Wrote\n ", fA, "\n ", fB, "\n ", fC, "\n")

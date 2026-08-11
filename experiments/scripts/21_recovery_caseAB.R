# Shared-scale recovery figures: Case A + Case B
# Writes: experiments/output/clgam_recovery_caseA.png
#         experiments/output/clgam_recovery_caseB.png
#
# From experiments/:  Rscript scripts/21_recovery_caseAB.R
# From project root:  Rscript experiments/scripts/21_recovery_caseAB.R

args_root <- if (dir.exists("clgam")) {
  "."
} else if (dir.exists("../clgam")) {
  ".."
} else {
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

nl_amp <- 1.2
nl_fun <- function(z) sin(2 * pi * z)
n_coarse <- 40L
n_fine_per <- 10L

# ---- fit both cases ----------------------------------------------------------
# Case A: fine covariate
dat_A <- simulate_ata(
  n_coarse = n_coarse, n_fine_per = n_fine_per, seed = 99L,
  include_covariate = TRUE, covariate_level = "fine",
  spatial_amp = 0.8, nl_amp = nl_amp, nl_fun = nl_fun
)
fit_A <- clgam(
  dat_A$y, dat_A$x1, dat_A$x2, dat_A$C,
  exposure = dat_A$efine, smooth = dat_A$nlcovfine,
  knots = c(12L, 12L), knots_nl = 12L,
  elements = TRUE, trace = FALSE
)

# Case B: coarse covariate — seed/amp/knots from grid search (dp≈0, RMSE≈0.05)
dat_B <- simulate_ata(
  n_coarse = n_coarse, n_fine_per = n_fine_per, seed = 6L,
  include_covariate = TRUE, covariate_level = "coarse",
  spatial_amp = 0.8, nl_amp = nl_amp, nl_fun = nl_fun
)
fit_B <- clgam(
  dat_B$y, dat_B$x1, dat_B$x2, dat_B$C,
  exposure = dat_B$efine, smooth = dat_B$nlcovfine,
  knots = c(10L, 10L), knots_nl = 6L,
  elements = TRUE, trace = FALSE
)

prep <- function(dat, fit) {
  nl_hat <- fit$nleffects[, 1] - mean(fit$nleffects[, 1])
  f_hat <- fit$eta - fit$nleffects[, 1]
  f_hat <- f_hat - mean(f_hat)
  f_true <- dat$eta_spatial_true - mean(dat$eta_spatial_true)
  list(
    dat = dat, fit = fit, case = dat$case,
    nl_hat = nl_hat,
    nl_true_fine = if (identical(dat$case, "B")) dat$h_true else {
      dat$g_true - mean(dat$g_true)
    },
    se = fit$sdnleffects[, 1],
    f_hat = f_hat, f_true = f_true,
    eta_hat = fit$eta - mean(fit$eta),
    eta_true = dat$eta_true - mean(dat$eta_true),
    mse_eta = mean((fit$eta - dat$eta_true)^2)
  )
}

A <- prep(dat_A, fit_A)
B <- prep(dat_B, fit_B)

smooth_curve <- function(obj, kn) {
  dat <- obj$dat
  fit <- obj$fit
  z_obs <- fit$nlcovfine[, 1]
  xl <- min(z_obs) - 0.01
  xr <- max(z_obs) + 0.01
  MMo <- mm_basis(z_obs, xl, xr, kn, 3, 2, decom = 2)
  n_nl <- ncol(MMo$Z)
  n_sp <- ncol(fit$matlist$Z) - n_nl
  beta <- fit$b.fixed[length(fit$b.fixed)]
  u <- fit$b.random[(n_sp + 1):(n_sp + n_nl)]

  # Evaluate only on observed support (avoid P-spline edge extrapolation)
  if (identical(obj$case, "B")) {
    z_lo <- min(dat$z_a)
    z_hi <- max(dat$z_a)
  } else {
    z_lo <- min(z_obs)
    z_hi <- max(z_obs)
  }
  z_grid <- seq(z_lo, z_hi, length.out = 300)
  MMg <- mm_basis(z_grid, xl, xr, kn, 3, 2, decom = 2)
  h_hat <- as.numeric(MMg$X[, -1, drop = FALSE] %*% beta + MMg$Z %*% u)
  h_hat <- h_hat - mean(fit$nleffects[, 1])

  if (identical(obj$case, "B")) {
    idx1 <- match(seq_len(dat$n_coarse), dat$coarse_id)
    za <- dat$z_a
    oa <- order(za)
    h_hat_a <- obj$nl_hat[idx1]
    se_a <- obj$se[idx1]
    # Curve through fitted values at observed z_a (demo-relevant support).
    # Dense mm_basis evaluation between sparse unique z_a can look phase-shifted
    # even when pointwise recovery at z_a is excellent.
    sf_hat <- stats::splinefun(za[oa], h_hat_a[oa], method = "natural")
    sf_se <- stats::splinefun(za[oa], se_a[oa], method = "natural")
    h_hat <- sf_hat(z_grid)
    se_grid <- pmax(sf_se(z_grid), 0)
    center <- mean(nl_amp * nl_fun(za))
    h_true <- nl_amp * nl_fun(z_grid) - center
    pts_z <- za
    pts_true <- dat$h_true_coarse
    cor_nl <- cor(h_hat_a, pts_true)
  } else {
    oa <- order(z_obs)
    se_grid <- approx(z_obs[oa], obj$se[oa], xout = z_grid, rule = 2)$y
    h_true <- nl_amp * nl_fun(z_grid) - mean(nl_amp * nl_fun(z_obs))
    pts_z <- NULL
    pts_true <- NULL
    cor_nl <- cor(obj$nl_hat, obj$nl_true_fine)
  }
  list(
    z = z_grid, hat = h_hat, true = h_true, se = se_grid,
    pts_z = pts_z, pts_true = pts_true, cor_nl = cor_nl,
    xlim = c(0, 1)  # shared axis; curves only drawn on support
  )
}

cur_A <- smooth_curve(A, kn = 12L)
cur_B <- smooth_curve(B, kn = 6L)

# ---- shared scales / palettes (identical across both PNGs) -------------------
pal_div <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))
pal_seq <- grDevices::colorRampPalette(c("#fff7bc", "#fec44f", "#d95f0e", "#7f0000"))

lim_f <- {
  m <- max(abs(c(A$f_true, A$f_hat, B$f_true, B$f_hat)), na.rm = TRUE)
  c(-m, m)
}
lim_yf <- range(c(A$dat$y_fine, B$dat$y_fine), finite = TRUE)
lim_y <- range(c(A$dat$y, B$dat$y), finite = TRUE)
lim_mu <- range(c(A$dat$y, A$fit$mu, B$dat$y, B$fit$mu), finite = TRUE)
lim_eta <- {
  m <- max(abs(c(A$eta_true, A$eta_hat, B$eta_true, B$eta_hat)), na.rm = TRUE)
  c(-m, m)
}
lim_nl <- {
  # Curves only (not SE edges) — wide boundary SE inflated ylim and made ĥ look shifted
  r <- range(c(cur_A$true, cur_A$hat, cur_B$true, cur_B$hat), finite = TRUE)
  r + c(-1, 1) * 0.08 * diff(r)
}
# True additive effect map scale (unused if panel 2 is the curve; kept for reuse)
lim_gmap <- {
  m <- max(abs(c(A$nl_true_fine, B$nl_true_fine)), na.rm = TRUE)
  c(-m, m)
}
xlim_z <- c(0, 1)

cols_from <- function(v, lim, pal, n = 101L) {
  vv <- pmin(pmax(v, lim[1]), lim[2])
  br <- seq(lim[1], lim[2], length.out = n)
  pal(n)[cut(vv, breaks = br, include.lowest = TRUE)]
}

plot_chor_f <- function(sf_f, sf_c, v, lim, pal, main, border_fine = NA) {
  bb <- sf::st_bbox(sf_f)
  xlim <- c(bb[["xmin"]], bb[["xmax"]])
  ylim <- c(bb[["ymin"]], bb[["ymax"]])
  plot(
    NA, type = "n",
    xlim = xlim, ylim = ylim,
    xlab = "", ylab = "", main = main,
    xaxs = "i", yaxs = "i", asp = NA
  )
  plot(sf::st_geometry(sf_f), col = cols_from(v, lim, pal), border = border_fine, add = TRUE)
  plot(sf::st_geometry(sf_c), border = "grey20", lwd = 0.55, add = TRUE)
  box(col = "grey35", lwd = 0.9)
}

plot_smooth_panel <- function(obj, cur, ylab, main, show_hat = TRUE) {
  dat <- obj$dat
  xlab <- if (identical(obj$case, "B")) {
    expression(z[a] ~ "(coarse)")
  } else {
    expression(z ~ "(fine)")
  }
  plot(cur$z, cur$true, type = "n", xlim = xlim_z, ylim = lim_nl,
       xlab = xlab, ylab = ylab, main = main, xaxs = "i", yaxs = "i")
  if (isTRUE(show_hat)) {
    lo <- pmax(cur$hat - 1.96 * cur$se, lim_nl[1])
    hi <- pmin(cur$hat + 1.96 * cur$se, lim_nl[2])
    polygon(c(cur$z, rev(cur$z)), c(lo, rev(hi)),
            col = grDevices::adjustcolor("#b2182b", 0.18), border = NA)
  }
  lines(cur$z, cur$true, lwd = 2.2, col = "grey20")
  if (isTRUE(show_hat)) {
    lines(cur$z, cur$hat, lwd = 2.4, col = "#b2182b")
  }
  if (identical(obj$case, "B") && !is.null(cur$pts_z)) {
    idx1 <- match(seq_len(dat$n_coarse), dat$coarse_id)
    oa <- order(cur$pts_z)
    points(cur$pts_z[oa], cur$pts_true[oa], pch = 1, cex = 0.65, col = "grey20")
    if (isTRUE(show_hat)) {
      points(cur$pts_z[oa], obj$nl_hat[idx1][oa], pch = 16, cex = 0.5, col = "#b2182b")
    }
  }
  if (isTRUE(show_hat)) {
    if (identical(obj$case, "B")) {
      legend("topright",
             legend = c(expression(h[true](z)), expression(hat(h)(z))),
             col = c("grey20", "#b2182b"), lwd = c(2.2, 2.4),
             bty = "n", cex = 0.75)
    } else {
      legend("topright",
             legend = c(expression(g[true](z)), expression(hat(g)(z)), "95%"),
             col = c("grey20", "#b2182b", grDevices::adjustcolor("#b2182b", 0.35)),
             lwd = c(2.2, 2.4, 8), bty = "n", cex = 0.75)
    }
  } else {
    lab <- if (identical(obj$case, "B")) expression(h[true](z)) else expression(g[true](z))
    legend("topright", legend = lab, col = "grey20", lwd = 2.2, bty = "n", cex = 0.75)
  }
  box(col = "grey35", lwd = 0.8)
}

# Layout 2×4:
# 1 true spatial | 2 true smooth | 3 y_f ~ Pois(spatial+smooth) | 4 observed y
# 5 spatial hat  | 6 smooth vs true | 7 aggregate fit | 8 eta vs truth
plot_recovery <- function(obj, cur, file, ylab_smooth) {
  dat <- obj$dat
  fit <- obj$fit
  sf_f <- dat$sf_fine
  sf_c <- dat$sf_coarse
  is_B <- identical(obj$case, "B")

  png(file, width = 1800, height = 900, res = 140)
  on.exit(dev.off(), add = TRUE)
  # Equal cells; common margins → aligned plot boxes
  layout(matrix(1:8, nrow = 2, byrow = TRUE),
         widths = rep(1, 4), heights = rep(1, 2))
  # Slightly tighter left margin: maps have no ylab text but keep same frame
  mar_common <- c(3.6, 3.6, 2.8, 1.0)
  par(mar = mar_common, oma = c(0, 0, 0, 0),
      mgp = c(2.0, 0.6, 0), tcl = -0.3, font.main = 2)

  # ---- row 1: generative truth / data ----------------------------------------
  plot_chor_f(sf_f, sf_c, obj$f_true, lim_f, pal_div, "1. True spatial f")

  plot_smooth_panel(
    obj, cur, ylab_smooth,
    main = if (is_B) {
      expression(paste("2. True smooth  ", h(z[a])))
    } else {
      expression(paste("2. True smooth  ", g(z)))
    },
    show_hat = FALSE
  )

  # 3–4. counts (same mar / axes as other panels)
  bb <- sf::st_bbox(sf_f)
  plot(NA, type = "n",
       xlim = c(bb[["xmin"]], bb[["xmax"]]),
       ylim = c(bb[["ymin"]], bb[["ymax"]]),
       xlab = "", ylab = "", xaxs = "i", yaxs = "i", asp = NA,
       main = expression(paste(
         "3. Intermediate  ", y[f], " ~ Poisson(", italic(f) + italic(s), ")"
       )))
  plot(sf::st_geometry(sf_f), col = cols_from(dat$y_fine, lim_yf, pal_seq),
       border = NA, add = TRUE)
  plot(sf::st_geometry(sf_c), border = "grey20", lwd = 0.55, add = TRUE)
  box(col = "grey35", lwd = 0.9)

  bb_c <- sf::st_bbox(sf_c)
  plot(NA, type = "n",
       xlim = c(bb_c[["xmin"]], bb_c[["xmax"]]),
       ylim = c(bb_c[["ymin"]], bb_c[["ymax"]]),
       xlab = "", ylab = "", xaxs = "i", yaxs = "i", asp = NA,
       main = expression(paste("4. Observed  ", y == C * y[f])))
  plot(sf::st_geometry(sf_c), col = cols_from(dat$y, lim_y, pal_seq),
       border = "grey25", lwd = 0.8, add = TRUE)
  box(col = "grey35", lwd = 0.9)

  # ---- row 2: recovery -------------------------------------------------------
  plot_chor_f(
    sf_f, sf_c, obj$f_hat, lim_f, pal_div,
    sprintf("5. Spatial f hat  cor=%.2f", cor(obj$f_hat, obj$f_true))
  )

  plot_smooth_panel(
    obj, cur, ylab_smooth,
    main = sprintf("6. Additive smooth vs true  cor=%.2f", cur$cor_nl),
    show_hat = TRUE
  )

  plot(dat$y, fit$mu, pch = 16, cex = 0.85,
       col = grDevices::adjustcolor("darkgreen", 0.75),
       xlim = lim_mu, ylim = lim_mu, xaxs = "i", yaxs = "i",
       xlab = "y (observed coarse)", ylab = expression(hat(mu)),
       main = sprintf("7. Aggregate fit  cor=%.2f", cor(fit$mu, dat$y)))
  abline(0, 1, lty = 2, col = "grey40")
  box(col = "grey35", lwd = 0.8)

  plot(obj$eta_true, obj$eta_hat, pch = 16, cex = 0.45,
       col = grDevices::adjustcolor("steelblue", 0.5),
       xlim = lim_eta, ylim = lim_eta, xaxs = "i", yaxs = "i",
       xlab = expression(eta[true]), ylab = expression(hat(eta)),
       main = sprintf("8. eta vs truth  MSE=%.3f  cor=%.2f",
                      obj$mse_eta, cor(obj$eta_hat, obj$eta_true)))
  abline(0, 1, lty = 2, col = "grey40")
  box(col = "grey35", lwd = 0.8)
}

fA <- file.path(out_dir, "clgam_recovery_caseA.png")
fB <- file.path(out_dir, "clgam_recovery_caseB.png")
plot_recovery(A, cur_A, fA, "g(z)")
plot_recovery(B, cur_B, fB, "h(z)")

cat("Shared scales (identical in both figures):\n")
cat(sprintf("  lim_f    = [%.3f, %.3f]  (panels 1,5)\n", lim_f[1], lim_f[2]))
cat(sprintf("  lim_yf   = [%.0f, %.0f]  (panel 3)\n", lim_yf[1], lim_yf[2]))
cat(sprintf("  lim_y    = [%.0f, %.0f]  (panel 4)\n", lim_y[1], lim_y[2]))
cat(sprintf("  lim_nl   = [%.3f, %.3f]  (panels 2,6)\n", lim_nl[1], lim_nl[2]))
cat(sprintf("  lim_mu   = [%.1f, %.1f]  (panel 7)\n", lim_mu[1], lim_mu[2]))
cat(sprintf("  lim_eta  = [%.3f, %.3f]  (panel 8)\n", lim_eta[1], lim_eta[2]))
cat("Wrote\n ", fA, "\n ", fB, "\n")

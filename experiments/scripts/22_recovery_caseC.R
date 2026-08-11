# Case C recovery figure (2×4 layout, shared style with Case A/B)
# From project root: Rscript experiments/scripts/22_recovery_caseC.R

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

nl_funs <- list(
  function(z) -tanh(2.5 * (z - 0.5)),
  function(z) tanh(2.5 * (z - 0.5))
)

# Easy Case C demo: descending + ascending (sine is weakly identified under ATA)
dat <- simulate_ata(
  n_coarse = 40L, n_fine_per = 10L, seed = 6L,
  include_covariate = TRUE, covariate_level = "both",
  spatial_amp = 0.75, nl_amp = c(1.2, 1.2),
  nl_fun = nl_funs
)
fit <- clgam(
  dat$y, dat$x1, dat$x2, dat$C,
  exposure = dat$efine, smooth = dat$nlcovfine,
  knots = c(10L, 10L), knots_nl = c(8L, 6L),
  elements = TRUE, trace = FALSE
)

g_hat <- fit$nleffects[, 1] - mean(fit$nleffects[, 1])
h_hat <- fit$nleffects[, 2] - mean(fit$nleffects[, 2])
g_true <- dat$g_true - mean(dat$g_true)
h_true <- dat$h_true - mean(dat$h_true)
f_hat <- fit$eta - fit$nleffects[, 1] - fit$nleffects[, 2]
f_hat <- f_hat - mean(f_hat)
f_true <- dat$eta_spatial_true - mean(dat$eta_spatial_true)
eta_hat <- fit$eta - mean(fit$eta)
eta_true <- dat$eta_true - mean(dat$eta_true)
mse_eta <- mean((fit$eta - dat$eta_true)^2)

idx1 <- match(seq_len(dat$n_coarse), dat$coarse_id)
h_hat_a <- h_hat[idx1]
h_true_a <- dat$h_true_coarse
cor_g <- cor(g_hat, g_true)
cor_h <- cor(h_hat_a, h_true_a)
cor_f <- cor(f_hat, f_true)

cat(sprintf("Case C: cor_f=%.3f cor_g=%.3f cor_h=%.3f cor_eta=%.3f MSE=%.3f\n",
            cor_f, cor_g, cor_h, cor(eta_hat, eta_true), mse_eta))

# Dense curves for g (fine) and h (coarse support)
z_grid <- seq(0, 1, length.out = 300)
g_true_c <- 1.1 * sin(2 * pi * z_grid) - mean(1.1 * sin(2 * pi * dat$z_f))
h_true_c <- 1.1 * cos(2 * pi * z_grid) - mean(1.1 * cos(2 * pi * dat$z_a))

# Hat curves via spline through observed evaluations
og <- order(dat$z_f)
g_hat_c <- approx(dat$z_f[og], g_hat[og], xout = z_grid, rule = 2)$y
za <- dat$z_a
oa <- order(za)
z_h <- seq(min(za), max(za), length.out = 300)
h_hat_c <- stats::splinefun(za[oa], h_hat_a[oa], method = "natural")(z_h)
h_true_h <- 1.1 * cos(2 * pi * z_h) - mean(1.1 * cos(2 * pi * za))
se_g <- approx(dat$z_f[og], fit$sdnleffects[og, 1], xout = z_grid, rule = 2)$y
se_h <- stats::splinefun(za[oa], fit$sdnleffects[idx1, 2][oa], method = "natural")(z_h)

pal_div <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))
pal_seq <- grDevices::colorRampPalette(c("#fff7bc", "#fec44f", "#d95f0e", "#7f0000"))
lim_f <- { m <- max(abs(c(f_true, f_hat))); c(-m, m) }
lim_yf <- range(dat$y_fine)
lim_y <- range(dat$y)
lim_mu <- range(c(dat$y, fit$mu))
lim_eta <- { m <- max(abs(c(eta_true, eta_hat))); c(-m, m) }
lim_nl <- range(c(g_true_c, g_hat_c, h_true_c, h_hat_c, h_true_a, h_hat_a),
                finite = TRUE)
lim_nl <- lim_nl + c(-1, 1) * 0.08 * diff(lim_nl)

cols_from <- function(v, lim, pal, n = 101L) {
  vv <- pmin(pmax(v, lim[1]), lim[2])
  br <- seq(lim[1], lim[2], length.out = n)
  pal(n)[cut(vv, breaks = br, include.lowest = TRUE)]
}

plot_chor <- function(sf_f, sf_c, v, lim, pal, main) {
  bb <- sf::st_bbox(sf_f)
  plot(NA, type = "n",
       xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
       xlab = "", ylab = "", main = main, xaxs = "i", yaxs = "i", asp = NA)
  plot(sf::st_geometry(sf_f), col = cols_from(v, lim, pal), border = NA, add = TRUE)
  plot(sf::st_geometry(sf_c), border = "grey20", lwd = 0.55, add = TRUE)
  box(col = "grey35", lwd = 0.9)
}

sf_f <- dat$sf_fine
sf_c <- dat$sf_coarse
file <- file.path(out_dir, "clgam_recovery_caseC.png")

png(file, width = 1800, height = 900, res = 140)
layout(matrix(1:8, nrow = 2, byrow = TRUE), widths = rep(1, 4), heights = rep(1, 2))
par(mar = c(3.6, 3.6, 2.8, 1.0), oma = c(0, 0, 0, 0),
    mgp = c(2.0, 0.6, 0), tcl = -0.3, font.main = 2)

# 1 spatial true
plot_chor(sf_f, sf_c, f_true, lim_f, pal_div, "1. True spatial f")

# 2 true smooths g & h
plot(z_grid, g_true_c, type = "n", xlim = c(0, 1), ylim = lim_nl,
     xlab = "z", ylab = "effect", xaxs = "i", yaxs = "i",
     main = expression(paste("2. True smooths  ", g(z[f]), " & ", h(z[a]))))
lines(z_grid, g_true_c, lwd = 2.2, col = "grey20")
lines(z_h, h_true_h, lwd = 2.2, lty = 2, col = "grey20")
points(za[oa], h_true_a[oa], pch = 1, cex = 0.55, col = "grey40")
legend("topright",
       legend = c(expression(g[true](z[f])), expression(h[true](z[a])), expression(z[a])),
       col = "grey20", lwd = c(2.2, 2.2, NA), lty = c(1, 2, NA),
       pch = c(NA, NA, 1), pt.cex = 0.6, bty = "n", cex = 0.72)
box(col = "grey35", lwd = 0.9)

# 3–4 counts
bb <- sf::st_bbox(sf_f)
plot(NA, type = "n", xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
     xlab = "", ylab = "", xaxs = "i", yaxs = "i", asp = NA,
     main = expression(paste("3. Intermediate  ", y[f], " ~ Poisson(", italic(f) + italic(s), ")")))
plot(sf::st_geometry(sf_f), col = cols_from(dat$y_fine, lim_yf, pal_seq), border = NA, add = TRUE)
plot(sf::st_geometry(sf_c), border = "grey20", lwd = 0.55, add = TRUE)
box(col = "grey35", lwd = 0.9)

bb_c <- sf::st_bbox(sf_c)
plot(NA, type = "n", xlim = c(bb_c["xmin"], bb_c["xmax"]), ylim = c(bb_c["ymin"], bb_c["ymax"]),
     xlab = "", ylab = "", xaxs = "i", yaxs = "i", asp = NA,
     main = expression(paste("4. Observed  ", y == C * y[f])))
plot(sf::st_geometry(sf_c), col = cols_from(dat$y, lim_y, pal_seq),
     border = "grey25", lwd = 0.8, add = TRUE)
box(col = "grey35", lwd = 0.9)

# 5 spatial hat
plot_chor(sf_f, sf_c, f_hat, lim_f, pal_div,
          sprintf("5. Spatial f hat  cor=%.2f", cor_f))

# 6 smooth vs true (both)
plot(z_grid, g_true_c, type = "n", xlim = c(0, 1), ylim = lim_nl,
     xlab = "z", ylab = "effect", xaxs = "i", yaxs = "i",
     main = sprintf("6. Smooths vs true  cor(g,h)=%.2f,%.2f", cor_g, cor_h))
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
legend("topright",
       legend = c(expression(g), expression(hat(g)), expression(h), expression(hat(h))),
       col = c("grey30", "#b2182b", "grey30", "#2166ac"),
       lwd = 2, lty = c(1, 1, 2, 2), bty = "n", cex = 0.7)
box(col = "grey35", lwd = 0.9)

# 7–8
plot(dat$y, fit$mu, pch = 16, cex = 0.85,
     col = adjustcolor("darkgreen", 0.75),
     xlim = lim_mu, ylim = lim_mu, xaxs = "i", yaxs = "i",
     xlab = "y (observed coarse)", ylab = expression(hat(mu)),
     main = sprintf("7. Aggregate fit  cor=%.2f", cor(fit$mu, dat$y)))
abline(0, 1, lty = 2, col = "grey40")
box(col = "grey35", lwd = 0.8)

plot(eta_true, eta_hat, pch = 16, cex = 0.45,
     col = adjustcolor("steelblue", 0.5),
     xlim = lim_eta, ylim = lim_eta, xaxs = "i", yaxs = "i",
     xlab = expression(eta[true]), ylab = expression(hat(eta)),
     main = sprintf("8. eta vs truth  MSE=%.3f  cor=%.2f",
                    mse_eta, cor(eta_hat, eta_true)))
abline(0, 1, lty = 2, col = "grey40")
box(col = "grey35", lwd = 0.8)

dev.off()
cat("Wrote", file, "\n")

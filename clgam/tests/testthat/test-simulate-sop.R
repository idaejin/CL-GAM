test_that("simulate_ata nested polygons + clgam", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 12L, n_fine_per = 6L, seed = 42L,
    spatial_amp = 0.8
  )
  expect_equal(dat$n_fine, 12L * 6L)
  expect_equal(length(dat$y), dat$n_coarse)
  expect_equal(nrow(dat$sf_coarse), dat$n_coarse)
  expect_equal(nrow(dat$sf_fine), dat$n_fine)
  expect_equal(as.numeric(colSums(dat$C)), rep(1, dat$n_fine))
  expect_equal(as.numeric(rowSums(dat$C)), rep(6, dat$n_coarse))
  expect_equal(dat$y, as.numeric(dat$C %*% dat$y_fine))
  expect_equal(sum(dat$y), sum(dat$y_fine))

  fit <- clgam(
    y = dat$y,
    coords = cbind(dat$x1, dat$x2),
    C = dat$C,
    exposure = dat$efine,
    knots = c(8L, 8L),
    elements = TRUE,
    trace = FALSE
  )
  expect_s3_class(fit, "clgam")
  expect_true(is.finite(fit$aic))
  expect_gt(cor(fit$eta, dat$eta_spatial_true), 0.2)

  pdf(NULL)
  expect_silent(plot(fit, which = c(1, 2, 6), ask = FALSE,
                     sf_fine = dat$sf_fine, sf_coarse = dat$sf_coarse))
  dev.off()
})

test_that("smooth covariate recovers nonlinear g(z)", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 40L, n_fine_per = 10L, seed = 21L,
    include_covariate = TRUE, nl_amp = 1.5, spatial_amp = 0.4,
    nl_fun = function(z) sin(2 * pi * z)
  )
  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    knots = c(8L, 8L), knots_nl = 14L,
    elements = FALSE, trace = FALSE
  )
  expect_equal(fit$nl.basis, "pspline")
  expect_true(fit$orth.smooth)
  expect_equal(length(fit$var.comp), 3L)
  expect_true(all(fit$var.comp > 0))
  g <- fit$nleffects[, 1] - mean(fit$nleffects[, 1])
  gt <- dat$g_true - mean(dat$g_true)
  z <- dat$nlcovfine[, 1]
  orth <- cor(resid(lm(g ~ z)), resid(lm(gt ~ z)))
  expect_gt(cor(g, gt), 0.7)
  expect_gt(orth, 0.7)

  pdf(NULL)
  expect_silent(plot(fit, which = 4, ask = FALSE, g_true = dat$g_true))
  dev.off()
})

test_that("simulate_ata partition and Poisson aggregation", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 8L, n_fine_per = 5L, seed = 5L)
  expect_equal(dat$n_fine, 40L)
  expect_equal(dat$y, as.numeric(dat$C %*% dat$y_fine))
  expect_equal(dat$mu_coarse, as.numeric(dat$C %*% dat$gamma), tolerance = 1e-10)
  expect_equal(dat$gamma, dat$efine * exp(dat$eta_true), tolerance = 1e-10)
  # nested: each fine centroid in its coarse polygon
  for (j in c(1L, 10L, 20L, 40L)) {
    expect_true(as.logical(sf::st_within(
      sf::st_centroid(sf::st_geometry(dat$sf_fine)[j]),
      sf::st_geometry(dat$sf_coarse)[dat$coarse_id[j]],
      sparse = FALSE
    )[1, 1]))
  }
})

test_that("Case B coarse covariate factors through C", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 40L, n_fine_per = 6L, seed = 42L,
    include_covariate = TRUE, covariate_level = "coarse",
    nl_amp = 0.8, spatial_amp = 1.0,
    nl_fun = function(z) sin(2 * pi * z)
  )
  expect_equal(dat$case, "B")
  expect_equal(length(dat$z_a), dat$n_coarse)
  expect_equal(dat$z, dat$z_a[dat$coarse_id])
  # piecewise constant within coarse units
  expect_true(all(vapply(split(dat$z, dat$coarse_id),
                         function(v) length(unique(v)) == 1L, logical(1))))

  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    smooth_level = dat$nl_level,
    knots = c(8L, 8L), knots_nl = 10L,
    elements = FALSE, trace = FALSE
  )
  expect_equal(fit$nl.level, "coarse")
  expect_false(isTRUE(fit$orth.info$applied)) # coarse-only → no A_f from smooths
  h <- fit$nleffects[, 1] - mean(fit$nleffects[, 1])
  ht <- dat$h_true - mean(dat$h_true)
  expect_gt(cor(h, ht), 0.7)
})

test_that("Case C fine+coarse covariates", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 40L, n_fine_per = 6L, seed = 11L,
    include_covariate = TRUE, covariate_level = "both",
    nl_amp = c(1.0, 1.0), spatial_amp = 0.6,
    nl_fun = list(
      function(z) sin(2 * pi * z),
      function(z) cos(2 * pi * z)
    )
  )
  expect_equal(dat$case, "C")
  expect_equal(ncol(dat$nlcovfine), 2L)
  expect_equal(dat$nlcovfine[, 2], dat$z_a[dat$coarse_id])
  expect_true(all(vapply(split(dat$nlcovfine[, 2], dat$coarse_id),
                         function(v) length(unique(v)) == 1L, logical(1))))

  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    smooth_level = dat$nl_level,
    knots = c(8L, 8L), knots_nl = c(8L, 6L),
    elements = FALSE, trace = FALSE
  )
  expect_equal(fit$nl.level, c("fine", "coarse"))
  expect_true(fit$orth.info$applied)
  expect_lt(fit$orth.info$max_abs_QZ, 1e-8)
  expect_equal(ncol(fit$nleffects), 2L)
  g <- fit$nleffects[, 1] - mean(fit$nleffects[, 1])
  h <- fit$nleffects[, 2] - mean(fit$nleffects[, 2])
  gt <- dat$g_true - mean(dat$g_true)
  ht <- dat$h_true - mean(dat$h_true)
  expect_gt(cor(g, gt), 0.5)
  expect_gt(cor(h, ht), 0.5)
})

test_that("S3 methods still work", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 9L, n_fine_per = 4L, seed = 3L)
  fit <- clgam(dat$y, dat$x1, dat$x2, dat$C, exposure = dat$efine,
               knots = c(5L, 5L), elements = TRUE, trace = FALSE)
  expect_equal(fitted(fit, type = "mu"), fit$mu)
  expect_equal(length(residuals(fit)), length(fit$y))
  expect_true(is.finite(AIC(fit)))
  expect_error(predict(fit, newdata = list(x = 1)), "newdata")
})

test_that("nl_fun presets sine/tanh match the paper curves", {
  skip_if_not_installed("sf")
  dat_sin <- simulate_ata(
    n_coarse = 8L, n_fine_per = 4L, seed = 4L,
    include_covariate = TRUE, nl_fun = "sine",
    spatial_amp = 0.4, nl_amp = 1.0, family = "none"
  )
  z <- dat_sin$z_f
  g <- dat_sin$g_true - mean(dat_sin$g_true)
  gs <- sin(2 * pi * z)
  gs <- gs - mean(gs)
  expect_gt(cor(g, gs), 0.999)

  dat_c <- simulate_ata(
    n_coarse = 8L, n_fine_per = 4L, seed = 5L,
    include_covariate = TRUE, covariate_level = "both",
    nl_fun = "tanh", nl_amp = c(1.2, 1.2),
    spatial_amp = 0.5, family = "none"
  )
  expect_equal(dat_c$case, "C")
  expect_equal(dat_c$nl_level, c("fine", "coarse"))
  gf <- dat_c$g_true - mean(dat_c$g_true)
  gt <- -tanh(2.5 * (dat_c$z_f - 0.5))
  gt <- gt - mean(gt)
  expect_gt(cor(gf, gt), 0.999)
  ht <- tanh(2.5 * (dat_c$z_a - 0.5))
  ht <- ht - mean(ht)
  h_coarse <- vapply(split(dat_c$h_true, dat_c$coarse_id), mean, 0)
  expect_gt(cor(h_coarse, ht), 0.999)
})

test_that("scenario presets fill paper A/C defaults but stay overridable", {
  skip_if_not_installed("sf")
  dat_a <- simulate_ata(scenario = "A", n_coarse = 8L, n_fine_per = 4L,
                        seed = 6L, family = "none")
  expect_equal(dat_a$scenario, "A")
  expect_equal(dat_a$case, "A")
  expect_equal(dat_a$n_fine, 32L)
  expect_equal(colnames(dat_a$nlcovfine), "z_f")
  z <- dat_a$z_f
  g <- dat_a$g_true - mean(dat_a$g_true)
  gs <- 1.2 * sin(2 * pi * z)
  gs <- gs - mean(gs)
  expect_gt(cor(g, gs), 0.999)

  dat_c <- simulate_ata(scenario = "C", n_coarse = 8L, n_fine_per = 4L,
                        seed = 7L, family = "none")
  expect_equal(dat_c$case, "C")
  expect_equal(ncol(dat_c$nlcovfine), 2L)
})

test_that("rho confounding raises kappa and stores f_perp", {
  skip_if_not_installed("sf")
  d0 <- simulate_ata(
    n_coarse = 10L, n_fine_per = 4L, seed = 11L,
    spatial_amp = 0.8, rho = 0, family = "none"
  )
  d9 <- simulate_ata(
    n_coarse = 10L, n_fine_per = 4L, seed = 11L,
    spatial_amp = 0.8, rho = 0.9, family = "none"
  )
  expect_equal(d0$n_fine, d9$n_fine)
  expect_equal(d0$x1, d9$x1)
  expect_true(is.finite(d0$kappa))
  expect_gt(d9$kappa, d0$kappa)
  expect_equal(length(d9$f_perp), d9$n_fine)
  expect_equal(
    as.numeric(d9$eta_true),
    as.numeric(d9$f_perp + d9$g_true),
    tolerance = 1e-10
  )
  expect_lt(abs(stats::sd(d0$g_true) - 0.45), 0.15)
})

test_that("Poisson DGP intercept, exposure_scale, and family='none'", {
  skip_if_not_installed("sf")
  d_mean <- simulate_ata(
    n_coarse = 8L, n_fine_per = 4L, seed = 2L, family = "none"
  )
  expect_equal(d_mean$family, "none")
  expect_equal(d_mean$y, d_mean$mu_coarse, tolerance = 1e-10)

  d0 <- simulate_ata(
    n_coarse = 8L, n_fine_per = 4L, seed = 2L,
    family = "none", intercept = 0
  )
  d1 <- simulate_ata(
    n_coarse = 8L, n_fine_per = 4L, seed = 2L,
    family = "none", intercept = 1
  )
  expect_equal(d1$gamma / d0$gamma, rep(exp(1), d0$n_fine), tolerance = 1e-10)

  de <- simulate_ata(
    n_coarse = 8L, n_fine_per = 4L, seed = 2L,
    family = "none", exposure_scale = 2
  )
  expect_equal(de$efine / d0$efine, rep(2, d0$n_fine), tolerance = 1e-10)

  d_pois <- simulate_ata(
    n_coarse = 8L, n_fine_per = 4L, seed = 2L, family = "poisson"
  )
  expect_true(all(d_pois$y_fine == as.integer(d_pois$y_fine)))
  expect_equal(d_pois$y, as.numeric(d_pois$C %*% d_pois$y_fine))
})

test_that("Matérn covariance is 1 at 0 and decreases with distance", {
  expect_equal(clgam:::.clgam_matern_cov(0, nu = 1, range = 0.3), 1)
  expect_equal(clgam:::.clgam_matern_cov(0, nu = 2, range = 0.3), 1)
  h <- seq(0, 1, length.out = 25)
  c1 <- clgam:::.clgam_matern_cov(h, nu = 1, range = 0.3)
  c2 <- clgam:::.clgam_matern_cov(h, nu = 2, range = 0.3)
  expect_true(all(diff(c1) <= 1e-10))
  expect_true(all(diff(c2) <= 1e-10))
  expect_gt(c2[3], c1[3])
})

test_that("matern1/matern2 truths are GP draws, not the P-spline field", {
  skip_if_not_installed("sf")
  d_ps <- simulate_ata(
    n_coarse = 8L, n_fine_per = 4L, seed = 8L,
    spatial_amp = 0.8, family = "none"
  )
  d1 <- simulate_ata(
    scenario = "matern1", n_coarse = 8L, n_fine_per = 4L,
    seed = 8L, family = "none"
  )
  d2 <- simulate_ata(
    scenario = "matern2", n_coarse = 8L, n_fine_per = 4L,
    seed = 8L, family = "none"
  )
  expect_equal(d1$spatial_truth, "matern")
  expect_equal(d1$matern_nu, 1)
  expect_equal(d2$matern_nu, 2)
  expect_true(is.null(d1$nlcovfine))
  expect_gt(
    mean((d1$eta_spatial_true - d_ps$eta_spatial_true)^2),
    1e-4
  )
  expect_gt(
    mean((d1$eta_spatial_true - d2$eta_spatial_true)^2),
    1e-6
  )
  expect_true(is.finite(stats::sd(d1$eta_spatial_true)))
})

test_that("jump truth is piecewise constant on two coarse groups", {
  skip_if_not_installed("sf")
  d <- simulate_ata(
    scenario = "jump", n_coarse = 10L, n_fine_per = 4L,
    seed = 9L, family = "none", spatial_amp = 0.8
  )
  expect_equal(d$spatial_truth, "jump")
  expect_equal(d$jump_amp, 1.6)
  expect_equal(length(unique(d$eta_spatial_true)), 2L)
  expect_true(all(vapply(
    split(d$eta_spatial_true, d$coarse_id),
    function(v) length(unique(round(v, 10))) == 1L,
    logical(1)
  )))
  expect_equal(
    diff(sort(unique(d$eta_spatial_true))),
    d$jump_amp,
    tolerance = 1e-10
  )
  expect_equal(length(unique(d$jump_group)), 2L)
  expect_true(all(d$jump_group[d$coarse_id] == d$sf_fine$jump_group))
})



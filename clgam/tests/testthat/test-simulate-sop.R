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

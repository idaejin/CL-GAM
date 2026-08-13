test_that("family=poisson is default and stores phi", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 16L, n_fine_per = 4L, seed = 7L,
    spatial_amp = 0.6
  )
  fit <- clgam(
    y = dat$y, coords = cbind(dat$x1, dat$x2), C = dat$C,
    exposure = dat$efine, knots = c(6L, 6L),
    elements = TRUE, trace = FALSE
  )
  expect_identical(fit$family, "poisson")
  expect_identical(fit$type, "spatial")
  expect_true(is.finite(fit$phi))
  expect_true(fit$phi > 0)
  expect_true(is.finite(fit$deviance_df))
})

test_that("quasipoisson matches poisson means and scales SEs by sqrt(phi)", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 20L, n_fine_per = 5L, seed = 11L,
    spatial_amp = 0.7
  )
  fit_p <- clgam(
    y = dat$y, coords = cbind(dat$x1, dat$x2), C = dat$C,
    exposure = dat$efine, knots = c(6L, 6L),
    family = poisson(), elements = TRUE, trace = FALSE
  )
  fit_q <- clgam(
    y = dat$y, coords = cbind(dat$x1, dat$x2), C = dat$C,
    exposure = dat$efine, knots = c(6L, 6L),
    family = quasipoisson(), elements = TRUE, trace = FALSE
  )
  expect_identical(fit_q$family, "quasipoisson")
  expect_equal(fit_q$eta, fit_p$eta, tolerance = 1e-10)
  expect_equal(fit_q$mu, fit_p$mu, tolerance = 1e-10)
  expect_equal(fit_q$var.comp, fit_p$var.comp, tolerance = 1e-10)
  expect_equal(fit_q$phi, fit_p$phi, tolerance = 1e-10)
  expect_gt(fit_q$phi, 0)
  expect_equal(fit_q$sd.eta, fit_p$sd.eta * sqrt(fit_q$phi), tolerance = 1e-10)
})

test_that("family string and invalid family work", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 10L, n_fine_per = 3L, seed = 3L)
  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, knots = c(5L, 5L),
    family = "quasipoisson", elements = FALSE, trace = FALSE
  )
  expect_identical(fit$family, "quasipoisson")
  expect_error(
    clgam(
      dat$y, dat$x1, dat$x2, dat$C,
      exposure = dat$efine, knots = c(5L, 5L),
      family = "negbin", elements = FALSE, trace = FALSE
    ),
    "poisson"
  )
})

test_that("C = I quasipoisson SEs align with glm scale inflation", {
  # Identity composition: coarse = fine. Compare SE inflation to glm().
  set.seed(99)
  n <- 80L
  x <- runif(n)
  eta <- 0.3 + 0.8 * sin(2 * pi * x)
  mu <- exp(eta)
  # Overdispersed counts via NB (mean mu, variance ~ 2*mu)
  y <- rnbinom(n, mu = mu, size = mu) # Var = mu + mu^2/mu = 2 mu
  C <- diag(n)
  ef <- rep(1, n)
  # Spatial coords on a line embedded in 2D (weak spatial field)
  fit_p <- pois_SOP(
    y = y, x1 = x, x2 = rep(0.5, n), efine = ef, C = C,
    ndx = c(8L, 3L), elements = TRUE, trace = FALSE,
    family = poisson(), maxit = c(40L, 40L)
  )
  fit_q <- pois_SOP(
    y = y, x1 = x, x2 = rep(0.5, n), efine = ef, C = C,
    ndx = c(8L, 3L), elements = TRUE, trace = FALSE,
    family = quasipoisson(), maxit = c(40L, 40L)
  )
  expect_equal(fit_q$eta, fit_p$eta, tolerance = 1e-10)
  expect_equal(
    mean(fit_q$sd.eta / fit_p$sd.eta),
    sqrt(fit_q$phi),
    tolerance = 1e-8
  )
  # Pearson phi should detect overdispersion
  expect_gt(fit_q$phi, 1.2)
})

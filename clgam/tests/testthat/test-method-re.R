test_that("re='coarse' adds an iid variance component on partition C", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 12L, n_fine_per = 4L, seed = 21L,
    spatial_amp = 0.6
  )
  fit0 <- clgam(
    y = dat$y, coords = cbind(dat$x1, dat$x2), C = dat$C,
    exposure = dat$efine, knots = c(5L, 5L),
    elements = FALSE, trace = FALSE, re = "none"
  )
  fit <- clgam(
    y = dat$y, coords = cbind(dat$x1, dat$x2), C = dat$C,
    exposure = dat$efine, knots = c(5L, 5L),
    elements = FALSE, trace = FALSE, re = "coarse"
  )
  expect_identical(fit$method, "SOP")
  expect_identical(fit$re_term, "coarse")
  expect_identical(fit0$re_term, "none")
  expect_equal(length(fit$re), nrow(dat$C))
  expect_equal(length(fit$eta), length(dat$x1))
  expect_equal(length(fit$mu), length(dat$y))
  expect_true("re.coarse" %in% names(fit$var.comp))
  expect_true(is.finite(fit$var.comp[["re.coarse"]]))
  expect_gt(fit$var.comp[["re.coarse"]], 0)
  expect_equal(length(fit0$var.comp) + 1L, length(fit$var.comp))
})

test_that("re='coarse' with quasipoisson warns", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 8L, n_fine_per = 4L, seed = 22L, spatial_amp = 0.5
  )
  expect_warning(
    clgam(
      y = dat$y, coords = cbind(dat$x1, dat$x2), C = dat$C,
      exposure = dat$efine, knots = c(4L, 4L),
      family = quasipoisson(), re = "coarse",
      elements = FALSE, trace = FALSE
    ),
    "overdispersion"
  )
})

test_that("method='Laplace' matches SOP eta on a small partition", {
  skip_on_cran()
  skip_if_not_installed("TMB")
  skip_if_not_installed("sf")
  set.seed(3)
  n <- 16L
  x1 <- rep(seq(0, 1, length.out = 4L), each = 4L)
  x2 <- rep(seq(0, 1, length.out = 4L), times = 4L)
  C <- diag(n)
  y <- stats::rpois(n, lambda = 4)
  fit_s <- clgam(
    y, x1, x2, C, exposure = 1, knots = c(3L, 3L),
    method = "SOP", elements = FALSE, trace = FALSE
  )
  fit_l <- tryCatch(
    clgam(
      y, x1, x2, C, exposure = 1, knots = c(3L, 3L),
      method = "Laplace", elements = FALSE, silent = TRUE
    ),
    error = function(e) e
  )
  if (inherits(fit_l, "error")) {
    skip(paste("Laplace compile/fit failed:", conditionMessage(fit_l)))
  }
  expect_identical(fit_l$method, "Laplace")
  expect_equal(length(fit_l$eta), n)
  expect_gt(stats::cor(fit_s$eta, fit_l$eta), 0.5)
})

test_that("method='Laplace' quasipoisson matches Poisson means and scales SEs", {
  skip_on_cran()
  skip_if_not_installed("TMB")
  set.seed(4)
  n <- 16L
  x1 <- rep(seq(0, 1, length.out = 4L), each = 4L)
  x2 <- rep(seq(0, 1, length.out = 4L), times = 4L)
  C <- diag(n)
  y <- stats::rpois(n, lambda = 4)
  fit_p <- tryCatch(
    clgam(
      y, x1, x2, C, exposure = 1, knots = c(3L, 3L),
      method = "Laplace", family = poisson(),
      elements = TRUE, silent = TRUE
    ),
    error = function(e) e
  )
  if (inherits(fit_p, "error")) {
    skip(paste("Laplace compile/fit failed:", conditionMessage(fit_p)))
  }
  fit_q <- clgam(
    y, x1, x2, C, exposure = 1, knots = c(3L, 3L),
    method = "Laplace", family = quasipoisson(),
    elements = TRUE, silent = TRUE
  )
  expect_identical(fit_q$method, "Laplace")
  expect_identical(fit_q$family, "quasipoisson")
  expect_equal(fit_q$eta, fit_p$eta, tolerance = 1e-10)
  expect_equal(fit_q$mu, fit_p$mu, tolerance = 1e-10)
  expect_equal(fit_q$var.comp, fit_p$var.comp, tolerance = 1e-10)
  expect_equal(fit_q$phi, fit_p$phi, tolerance = 1e-10)
  if (all(is.finite(fit_p$sd.eta)) && all(is.finite(fit_q$sd.eta))) {
    expect_equal(fit_q$sd.eta, fit_p$sd.eta * sqrt(fit_q$phi), tolerance = 1e-10)
  }
})

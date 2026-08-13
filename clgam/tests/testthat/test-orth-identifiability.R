test_that(".orth_cols and Af builders enforce fine-only space", {
  set.seed(1)
  m <- 40L
  A <- cbind(rnorm(m), rnorm(m))
  B <- cbind(A[, 1] + rnorm(m, sd = 0.01), rnorm(m))
  Ap <- clgam:::.orth_cols(A, B)
  expect_lt(clgam:::.orth_cross_max(Ap, B), 1e-10)

  Bk <- list(matrix(rnorm(m * 3), m, 3), matrix(rnorm(m * 3), m, 3))
  Af_all <- clgam:::.build_orth_Af(Bk, c("fine", "fine"), NULL)
  Af_fine <- clgam:::.build_orth_Af(Bk, c("fine", "coarse"), NULL)
  expect_equal(ncol(Af_all), 6L)
  expect_equal(ncol(Af_fine), 3L)

  lin <- matrix(rnorm(m), m, 1)
  Af_lin <- clgam:::.build_orth_Af(NULL, NULL, lin)
  expect_equal(ncol(Af_lin), 1L)
  expect_null(clgam:::.build_orth_Af(Bk, c("coarse", "coarse"), NULL))

  parts <- clgam:::.build_orth_Af_parts(
    Bk, c("fine", "coarse"), lin,
    smooth_names = c("z_f", "z_a"), linear_names = "xlin"
  )
  expect_equal(names(parts), c("z_f", "xlin"))
  expect_equal(ncol(parts$z_f), 3L)
})

test_that("Case A orth.smooth yields near-zero Q'Xs and Q'Zs", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 20L, n_fine_per = 6L, seed = 7L,
    include_covariate = TRUE, nl_amp = 1.0, spatial_amp = 0.5,
    nl_fun = function(z) sin(2 * pi * z)
  )
  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    smooth_level = dat$nl_level,
    knots = c(6L, 6L), knots_nl = 8L,
    elements = FALSE, trace = FALSE
  )
  expect_true(fit$orth.smooth)
  expect_true(fit$orth.info$applied)
  expect_lt(fit$orth.info$max_abs_QX, 1e-8)
  expect_lt(fit$orth.info$max_abs_QZ, 1e-8)
})

test_that("linear fine covariates enter the orth space", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 16L, n_fine_per = 5L, seed = 9L,
    include_covariate = FALSE, spatial_amp = 0.7
  )
  z <- as.numeric(scale(dat$eta_spatial_true))
  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, linear = cbind(z = z),
    knots = c(6L, 6L), elements = FALSE, trace = FALSE,
    orth.smooth = TRUE
  )
  expect_true(fit$orth.info$applied)
  expect_lt(fit$orth.info$max_abs_QZ, 1e-8)
})

test_that("Case C excludes coarse smooth from orth space", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 20L, n_fine_per = 6L, seed = 11L,
    include_covariate = TRUE, covariate_level = "both",
    nl_amp = c(1.0, 1.0), spatial_amp = 0.5,
    nl_fun = list(function(z) sin(2 * pi * z), function(z) cos(2 * pi * z))
  )
  expect_equal(dat$nl_level, c("fine", "coarse"))
  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    smooth_level = dat$nl_level,
    knots = c(6L, 6L), knots_nl = c(6L, 6L),
    elements = FALSE, trace = FALSE
  )
  expect_equal(fit$nl.level, c("fine", "coarse"))
  expect_true(fit$orth.info$applied)
  expect_lt(fit$orth.info$max_abs_QZ, 1e-8)

  # Wrongly marking both as fine still fits, but uses a larger A_f
  fit_both <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    smooth_level = c("fine", "fine"),
    knots = c(6L, 6L), knots_nl = c(6L, 6L),
    elements = FALSE, trace = FALSE
  )
  expect_true(fit_both$orth.info$applied)
  # Identified spatial fields need not match when A_f changes
  f1 <- fit$eta - rowSums(fit$nleffects)
  f2 <- fit_both$eta - rowSums(fit_both$nleffects)
  expect_gt(mean((f1 - f2)^2), 0)
})

test_that("orth.smooth=FALSE leaves spatial design unrestricted", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 12L, n_fine_per = 5L, seed = 3L,
    include_covariate = TRUE, nl_amp = 0.8, spatial_amp = 0.6
  )
  fit0 <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    knots = c(5L, 5L), knots_nl = 6L,
    elements = FALSE, trace = FALSE, orth.smooth = FALSE
  )
  expect_false(fit0$orth.smooth)
  expect_false(isTRUE(fit0$orth.info$applied))
})

test_that("kappa_diagnostic matches QR formula for eq. 21", {
  set.seed(1)
  A <- cbind(runif(40), runif(40))
  f <- rnorm(40)
  kap <- kappa_diagnostic(f, A)
  f0 <- f - mean(f)
  fp <- clgam:::.orth_cols(cbind(f0), A)[, 1]
  expect_equal(kap, sum((f0 - fp)^2) / sum(f0^2), tolerance = 1e-12)
  expect_gte(kap, 0)
  expect_lte(kap, 1)
})

test_that("kappa is ~1 when f lies in span(A) and ~0 when orthogonal", {
  set.seed(2)
  A <- cbind(runif(50), runif(50))
  f_in <- as.numeric(A %*% c(1.1, -0.4))
  expect_gt(kappa_diagnostic(f_in, A, center = FALSE), 0.99)
  Q <- qr.Q(qr(A))
  f_orth <- rnorm(50)
  f_orth <- f_orth - as.numeric(Q %*% crossprod(Q, f_orth))
  expect_lt(kappa_diagnostic(f_orth, A, center = FALSE), 1e-10)
})

test_that("clgam method matches stored kappa and accepts f_raw (eq. 21)", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 16L, n_fine_per = 5L, seed = 11L,
    include_covariate = TRUE, nl_amp = 1.0, spatial_amp = 0.5,
    nl_fun = function(z) sin(2 * pi * z)
  )
  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    knots = c(6L, 6L), knots_nl = 8L,
    elements = FALSE, trace = FALSE, orth.smooth = TRUE
  )
  expect_true(is.finite(fit$orth.info$kappa))
  expect_gte(fit$orth.info$kappa, 0)
  expect_lte(fit$orth.info$kappa, 1)
  expect_equal(
    kappa_diagnostic(fit),
    fit$orth.info$kappa,
    tolerance = 1e-10
  )
  f_raw <- as.numeric(dat$eta_spatial_true)
  Af <- clgam:::.clgam_Af(fit)
  expect_equal(
    kappa_diagnostic(fit, f = f_raw),
    kappa_diagnostic(f_raw, Af),
    tolerance = 1e-10
  )
})

test_that("Rcpp Schur path matches pure R and is used when compiled", {
  skip_if_not(exists("sop_solve_schur_cpp", mode = "function"),
              "Rcpp kernels not compiled")

  set.seed(42)
  p <- 8L
  q <- 12L
  XtX <- crossprod(matrix(rnorm(30 * p), 30, p))
  ZtX <- matrix(rnorm(q * p), q, p)
  ZtZ <- crossprod(matrix(rnorm(25 * q), 25, q))
  ZtXtZ <- rbind(t(ZtX), ZtZ)
  u <- rnorm(p + q)
  G <- runif(q, 0.5, 2)

  r_out <- clgam:::.sop_solve_schur_R(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache = NULL)
  options(clgam.use_rcpp = TRUE)
  cpp_out <- clgam:::.sop_solve_schur(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache = NULL)

  expect_equal(cpp_out$b.fixed, r_out$b.fixed, tolerance = 1e-10)
  expect_equal(cpp_out$b.random, r_out$b.random, tolerance = 1e-10)
  expect_equal(cpp_out$dZtNZ, r_out$dZtNZ, tolerance = 1e-10)

  # Cached second call (inner SOP iteration reuse)
  cpp2 <- clgam:::.sop_solve_schur(
    XtX, ZtX, ZtZ, ZtXtZ, u, G * 1.05, cache = cpp_out$cache
  )
  r2 <- clgam:::.sop_solve_schur_R(
    XtX, ZtX, ZtZ, ZtXtZ, u, G * 1.05, cache = r_out$cache
  )
  expect_equal(cpp2$b.fixed, r2$b.fixed, tolerance = 1e-10)
  expect_equal(cpp2$b.random, r2$b.random, tolerance = 1e-10)
})

test_that("clgam fit is unchanged by Rcpp vs R Schur solver", {
  skip_if_not_installed("sf")
  skip_if_not(exists("sop_solve_schur_cpp", mode = "function"),
              "Rcpp kernels not compiled")

  dat <- simulate_ata(
    n_coarse = 24L, n_fine_per = 6L, seed = 17L,
    include_covariate = TRUE, spatial_amp = 0.6
  )
  args <- list(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, C = dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    knots = c(8L, 8L), knots_nl = 10L,
    elements = TRUE, trace = FALSE
  )
  options(clgam.use_rcpp = TRUE)
  fit_cpp <- do.call(clgam, args)
  options(clgam.use_rcpp = FALSE)
  fit_r <- do.call(clgam, args)

  expect_equal(fit_cpp$eta, fit_r$eta, tolerance = 1e-10)
  expect_equal(fit_cpp$mu, fit_r$mu, tolerance = 1e-10)
  expect_equal(fit_cpp$var.comp, fit_r$var.comp, tolerance = 1e-8)
  expect_equal(fit_cpp$aic, fit_r$aic, tolerance = 1e-6)
})

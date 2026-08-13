test_that("formula interface matches positional (C = I)", {
  set.seed(1)
  n <- 36L
  x1 <- rep(seq(0, 1, length.out = 6L), each = 6L)
  x2 <- rep(seq(0, 1, length.out = 6L), times = 6L)
  z <- stats::runif(n)
  C <- diag(n)
  ef <- rep(1, n)
  eta <- 0.4 * sin(2 * pi * x1) * cos(2 * pi * x2) + 0.7 * sin(2 * pi * z)
  y <- stats::rpois(n, lambda = exp(eta))

  fit_pos <- clgam(
    y, x1, x2, C,
    exposure = ef, smooth = cbind(z = z),
    knots = c(4L, 4L), knots_nl = 5L,
    elements = FALSE, trace = FALSE, orth.smooth = TRUE
  )
  fit_f <- clgam(
    y ~ s(x1, x2) + s(z),
    C = C, exposure = ef,
    knots = c(4L, 4L), knots_nl = 5L,
    elements = FALSE, trace = FALSE, orth.smooth = TRUE
  )
  expect_equal(fit_pos$eta, fit_f$eta, tolerance = 1e-8)
  expect_equal(fit_f$formula, y ~ s(x1, x2) + s(z))
  expect_true(isTRUE(fit_f$orth.smooth))
})

test_that("formula data= list supplies C and exposure", {
  set.seed(2)
  n <- 36L
  x1 <- rep(seq(0, 1, length.out = 6L), each = 6L)
  x2 <- rep(seq(0, 1, length.out = 6L), times = 6L)
  z <- stats::runif(n)
  C <- diag(n)
  efine <- rep(1, n)
  y <- stats::rpois(n, lambda = exp(0.3 * sin(2 * pi * x1)))
  dat <- list(y = y, x1 = x1, x2 = x2, z = z, C = C, efine = efine)

  fit <- clgam(
    y ~ s(x1, x2) + s(z),
    data = dat,
    knots = c(4L, 4L), knots_nl = 5L,
    elements = FALSE, trace = FALSE, orth.smooth = TRUE
  )
  expect_s3_class(fit, "clgam")
  expect_equal(length(fit$eta), n)
})

test_that("s(..., k=) overrides knots when knots is omitted", {
  set.seed(3)
  n <- 25L
  x1 <- rep(seq(0, 1, length.out = 5L), each = 5L)
  x2 <- rep(seq(0, 1, length.out = 5L), times = 5L)
  C <- diag(n)
  y <- stats::rpois(n, lambda = 2)
  fit <- clgam(
    y ~ s(x1, x2, k = c(4, 4)),
    C = C, exposure = 1,
    elements = FALSE, trace = FALSE
  )
  expect_equal(as.integer(fit$ndx), c(4L, 4L))
})

test_that("linear formula terms are rejected", {
  expect_error(
    clgam(y ~ s(x1, x2) + z, C = diag(4)),
    "Linear terms"
  )
})

test_that("Case B formula expands z_a and sets level='coarse'", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    scenario = "B", n_coarse = 8L, n_fine_per = 4L,
    seed = 1L, family = "none"
  )
  fit_f <- clgam(
    y ~ s(x1, x2) + s(z_a, level = "coarse"),
    data = dat,
    knots = c(5L, 5L), knots_nl = 6L,
    elements = FALSE, trace = FALSE, orth.smooth = FALSE
  )
  fit_pos <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    smooth_level = dat$nl_level,
    knots = c(5L, 5L), knots_nl = 6L,
    elements = FALSE, trace = FALSE, orth.smooth = FALSE
  )
  expect_equal(fit_f$nl.level, "coarse")
  expect_equal(fit_f$eta, fit_pos$eta, tolerance = 1e-8)
})

test_that("simulate_ata_scenarios lists every preset", {
  sc <- simulate_ata_scenarios()
  expect_true(is.data.frame(sc))
  expect_equal(
    sc$scenario,
    c("default", "A", "B", "C", "confounding", "matern1", "matern2", "jump")
  )
  expect_equal(sc$orth.smooth[sc$scenario == "A"], FALSE)
  expect_equal(sc$orth.smooth[sc$scenario == "B"], FALSE)
  expect_equal(sc$orth.smooth[sc$scenario == "C"], FALSE)
  expect_equal(sc$orth.smooth[sc$scenario == "confounding"], TRUE)
  expect_match(sc$formula[sc$scenario == "B"], "level = \"coarse\"")
  expect_equal(sc$formula[sc$scenario == "matern1"], "y ~ s(x1, x2)")
  expect_equal(sc$n_fine[sc$scenario == "A"], 400L)
  expect_equal(sc$n_fine[sc$scenario == "default"], 96L)
})

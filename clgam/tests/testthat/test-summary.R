test_that("summary.clgam reports phi, labelled ED, AIC, and convergence", {
  # Identity C: no sf needed (same construction as the quasi-Poisson tests).
  set.seed(20)
  n <- 60L
  x1 <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  eta <- 0.2 + 0.6 * sin(2 * pi * x1) + 0.5 * (z - 0.5)^2
  y <- rpois(n, lambda = exp(eta))
  C <- diag(n)
  nlcov <- cbind(z_f = z)
  fit <- pois_SOP(
    y = y, x1 = x1, x2 = x2, efine = 1, C = C,
    nlcovfine = nlcov, ndx = c(6L, 3L), ndxnl = 6L,
    elements = TRUE, trace = FALSE, family = poisson(),
    maxit = c(40L, 40L)
  )
  expect_true(isTRUE(fit$converged) || isFALSE(fit$converged))
  expect_true(is.finite(fit$tol))
  expect_equal(names(fit$edf), c("spatial.x1", "spatial.x2", "z_f"))
  expect_true("z_f" %in% names(fit$orth.info$kappa_by))
  expect_equal(
    unname(fit$orth.info$kappa_by["z_f"]),
    fit$orth.info$kappa,
    tolerance = 1e-10
  )

  s <- summary(fit)
  expect_s3_class(s, "summary.clgam")
  expect_true(is.finite(s$phi))
  expect_true(is.finite(s$aic))
  expect_true(is.finite(s$bic))
  expect_true("z_f" %in% s$terms$term)
  expect_true(all(c("spatial.x1", "spatial.x2") %in% s$terms$term))
  expect_true(is.finite(s$terms$ed[s$terms$term == "z_f"]))
  expect_true(is.finite(s$kappa_by["z_f"]))

  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "Converged")
  expect_match(out, "Pearson phi")
  expect_match(out, "kappa")
  expect_match(out, "AIC")
  expect_match(out, "BIC")
  expect_match(out, "z_f")
  expect_match(out, "spatial.x1")
  expect_match(out, "Rel\\. change in eta")
})

test_that("summary notes quasi-Poisson phi SE scaling", {
  set.seed(21)
  n <- 50L
  x1 <- runif(n)
  x2 <- runif(n)
  y <- rnbinom(n, mu = exp(0.3 + sin(2 * pi * x1)), size = 4)
  fit <- pois_SOP(
    y = y, x1 = x1, x2 = x2, efine = 1, C = diag(n),
    ndx = c(6L, 3L), elements = TRUE, trace = FALSE,
    family = quasipoisson(), maxit = c(40L, 40L)
  )
  s <- summary(fit)
  expect_identical(s$family, "quasipoisson")
  expect_gt(s$phi, 1)
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "Quasi-Poisson")
  expect_match(out, "sqrt\\(phi\\)")
  expect_match(out, "Pearson phi")
  expect_match(out, "AIC")
})

test_that("summary rebuilds labelled terms for older un-named fits", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 12L, n_fine_per = 4L, seed = 9L,
    include_covariate = TRUE, spatial_amp = 0.5
  )
  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    knots = c(5L, 5L), knots_nl = 6L,
    elements = TRUE, trace = FALSE
  )
  fit$orth.info$kappa_by <- NULL
  names(fit$var.comp) <- NULL
  names(fit$edf) <- NULL
  fit$converged <- NULL
  fit$tol <- NULL
  s <- summary(fit)
  expect_true("z_f" %in% s$terms$term)
  expect_true(is.finite(s$kappa_by["z_f"]))
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "z_f")
  expect_match(out, "AIC")
  expect_match(out, "Converged")
})

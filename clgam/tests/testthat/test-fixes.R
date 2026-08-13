test_that("AIC.clgam/BIC.clgam do not double-count ed for clgam_contrast fits", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 10L, n_fine_per = 5L, seed = 5L, spatial_amp = 0.6)

  # Two-population stack (e.g. sex): same fine geometry/composition for both
  # groups, independent Poisson draws, mirroring the manuscript's Section 4.4
  # setup (clgam_contrast()/pois_incat_SOP()).
  set.seed(101)
  y1 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma))
  y2 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma))

  fit <- clgam_contrast(
    y = c(y1, y2),
    x1 = c(dat$x1, dat$x1), x2 = c(dat$x2, dat$x2),
    exposure = c(dat$efine, dat$efine),
    C1 = dat$C, C2 = dat$C,
    knots = c(6L, 6L), elements = TRUE
  )

  expect_s3_class(fit, "clgam_contrast")
  expect_length(fit$aic, 2L)
  expect_false(is.null(fit$aic_total))
  expect_false(is.null(fit$bic_total))

  # The correct joint total uses ONE shared ed, not one per category.
  expect_equal(fit$aic_total, fit$dev[1] + fit$dev[2] + 2 * fit$ed)
  expect_equal(fit$bic_total, fit$dev[1] + fit$dev[2] + log(length(fit$y)) * fit$ed)

  # AIC.clgam()/BIC.clgam() must return the joint total, NOT sum(fit$aic)
  # (which double-counts ed as 4*ed instead of 2*ed -- the bug this test
  # guards against).
  expect_equal(AIC(fit), fit$aic_total)
  expect_equal(BIC(fit), fit$bic_total)
  expect_false(isTRUE(all.equal(AIC(fit), sum(fit$aic))))
})

test_that("structure='contrast' uses shared+difference (sum-to-zero)", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 10L, n_fine_per = 4L, seed = 8L, spatial_amp = 0.5)
  set.seed(102)
  y1 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma))
  y2 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma * 1.2))
  fit <- clgam_contrast(
    y = c(y1, y2),
    x1 = c(dat$x1, dat$x1), x2 = c(dat$x2, dat$x2),
    exposure = c(dat$efine, dat$efine),
    C1 = dat$C, C2 = dat$C,
    knots = c(5L, 5L), elements = TRUE,
    structure = "contrast"
  )
  expect_identical(fit$structure, "contrast")
  expect_equal(
    names(fit$var.comp),
    c("spatial.shared.x1", "spatial.shared.x2",
      "spatial.contrast.x1", "spatial.contrast.x2")
  )
  n2 <- length(dat$x1)
  expect_equal(fit$eta.diff, fit$eta[seq_len(n2)] - fit$eta[(n2 + 1):length(fit$eta)])
  expect_equal(
    fit$eta.shared,
    0.5 * (fit$eta[seq_len(n2)] + fit$eta[(n2 + 1):length(fit$eta)])
  )
  # Sum-to-zero on the stacked spatial deviation: mean contrast over groups
  # is the stored difference; shared is the average surface.
  expect_equal(length(fit$eta.diff), n2)
  expect_true(all(is.finite(fit$sd.dif)))
})

test_that("structure='independent' remains the default contrast model", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 8L, n_fine_per = 4L, seed = 9L, spatial_amp = 0.4)
  set.seed(103)
  y1 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma))
  y2 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma))
  fit <- clgam_contrast(
    y = c(y1, y2),
    x1 = c(dat$x1, dat$x1), x2 = c(dat$x2, dat$x2),
    exposure = c(dat$efine, dat$efine),
    C1 = dat$C, C2 = dat$C,
    knots = c(4L, 4L), elements = FALSE
  )
  expect_identical(fit$structure, "independent")
  expect_true(grepl("^spatial\\.g1", names(fit$var.comp)[1]))
})

test_that("pois_incat_SOP validates equal-sized stacked groups", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 6L, n_fine_per = 4L, seed = 2L, spatial_amp = 0.4)
  y1 <- as.numeric(dat$C %*% dat$y_fine)

  # Mismatched Ccat1/Ccat2 column counts must error, not silently misalign
  # the difference-surface computation (previously unchecked; only the
  # clgam_contrast() wrapper enforced n_fine %% 2 == 0).
  bad_C2 <- dat$C[, seq_len(ncol(dat$C) - 2L), drop = FALSE]
  expect_error(
    pois_incat_SOP(
      y = c(y1, y1), x1 = c(dat$x1, dat$x1), x2 = c(dat$x2, dat$x2),
      efine = c(dat$efine, dat$efine),
      cat = factor(rep(c("g1", "g2"), each = dat$n_fine)),
      Ccat1 = dat$C, Ccat2 = bad_C2,
      ndx = c(5L, 5L), elements = FALSE
    ),
    "ncol\\(Ccat1\\)"
  )
})

test_that("elements=TRUE with BOTH linear and smooth fine covariates: no error, consistent SE source", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 15L, n_fine_per = 6L, seed = 13L,
    include_covariate = TRUE, nl_amp = 1.0, spatial_amp = 0.5,
    nl_fun = function(z) sin(2 * pi * z)
  )
  set.seed(7)
  lin <- matrix(rnorm(dat$n_fine), dat$n_fine, 1)

  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, linear = lin, smooth = dat$nlcovfine,
    knots = c(6L, 6L), knots_nl = 8L,
    elements = TRUE, trace = FALSE
  )

  expect_true(is.finite(fit$aic))
  expect_false(is.null(fit$sdleffects))
  expect_false(is.null(fit$sdnleffects))
  expect_true(all(is.finite(fit$sdleffects)))
  expect_true(all(is.finite(fit$sdnleffects)))
  expect_true(all(fit$sdleffects >= 0))
  expect_true(all(fit$sdnleffects >= 0))

  # sdleffects should now come from the same PIRLS-weighted covariance (M1$S)
  # as sdnleffects when elements=TRUE, not the cruder unweighted-Gram
  # fallback -- regression check: refitting with elements=FALSE (fallback
  # covariance only) should generally give different (typically larger, less
  # precise) SEs for the linear effect.
  fit_fallback <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, linear = lin, smooth = dat$nlcovfine,
    knots = c(6L, 6L), knots_nl = 8L,
    elements = FALSE, trace = FALSE
  )
  expect_false(is.null(fit_fallback$sdleffects))
  expect_false(isTRUE(all.equal(fit$sdleffects, fit_fallback$sdleffects)))
})

test_that("diverged flag exists and is FALSE for a healthy fit", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 10L, n_fine_per = 5L, seed = 4L, spatial_amp = 0.5)
  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, knots = c(6L, 6L),
    elements = FALSE, trace = FALSE
  )
  expect_false(is.null(fit$diverged))
  expect_false(fit$diverged)
})

test_that("mm_basis(): eigen-based non-null eigenvalues match the old svd-based ones", {
  x <- seq(0, 1, length.out = 60)
  mm <- mm_basis(x = x, xl = -0.01, xr = 1.01, ndx = 12, bdeg = 3, pord = 2, decom = 1)
  # Regression check against direct svd() of the same penalty, independent of
  # which internal routine mm_basis() now uses.
  D <- diff(diag(ncol(mm$B)), differences = 2)
  ref <- svd(crossprod(D))
  m <- ncol(mm$B)
  expect_equal(sort(mm$d), sort(ref$d[1:(m - 2)]), tolerance = 1e-8)
  expect_equal(ncol(mm$Z), m - 2)
})

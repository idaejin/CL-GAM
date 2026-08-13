test_that("cached PIRLS elements match fresh inv_bblock2 (pois_SOP)", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 32L, n_fine_per = 8L, seed = 99L,
    include_covariate = TRUE, spatial_amp = 0.7
  )
  fit <- clgam(
    dat$y, dat$x1, dat$x2, dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    knots = c(10L, 10L), knots_nl = 12L,
    elements = TRUE, trace = FALSE
  )
  expect_true(is.finite(fit$aic))
  expect_true(is.finite(fit$ed))
  expect_equal(length(fit$sd.eta), dat$n_fine)
  expect_true(all(fit$sd.eta > 0))
  expect_false(is.null(fit$sdnleffects))
})

test_that("elements cache does not change AIC or SEs vs reference fit", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 24L, n_fine_per = 6L, seed = 55L,
    include_covariate = TRUE
  )
  args <- list(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, C = dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    knots = c(8L, 8L), knots_nl = 10L,
    elements = TRUE, trace = FALSE,
    parold = c(0.8, 1.2), paroldnl = c(1.1)
  )
  fit1 <- do.call(clgam, args)
  fit2 <- do.call(clgam, args)
  expect_equal(fit1$eta, fit2$eta, tolerance = 1e-12)
  expect_equal(fit1$aic, fit2$aic, tolerance = 1e-8)
  expect_equal(fit1$sd.eta, fit2$sd.eta, tolerance = 1e-10)
  expect_equal(fit1$sdnleffects, fit2$sdnleffects, tolerance = 1e-10)
})

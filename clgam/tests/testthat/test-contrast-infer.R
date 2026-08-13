test_that("contrast inference: joint SE vs independence approx", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 10L, n_fine_per = 4L, seed = 21L, spatial_amp = 0.6)
  set.seed(104)
  y1 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma))
  y2 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma * 1.3))
  args <- list(
    y = c(y1, y2),
    x1 = c(dat$x1, dat$x1), x2 = c(dat$x2, dat$x2),
    exposure = c(dat$efine, dat$efine),
    C1 = dat$C, C2 = dat$C,
    knots = c(5L, 5L), elements = TRUE
  )
  fit_ind <- do.call(clgam_contrast, args)
  fit_con <- do.call(clgam_contrast, c(args, list(structure = "contrast")))

  expect_true(all(is.finite(fit_ind$sd.dif)))
  expect_true(all(is.finite(fit_ind$sd.shared)))
  expect_true(all(is.finite(fit_con$sd.dif)))
  expect_true(all(is.finite(fit_con$sd.shared)))
  expect_equal(length(fit_ind$sd.shared), length(dat$x1))
  expect_equal(length(fit_con$sd.shared), length(dat$x1))

  # Independent random fields: joint se(diff) may match the independence
  # approx (posterior Cov(eta1, eta2) ~ 0). A shared field induces
  # Corr > 0, so sd.dif2 overstates se(eta1-eta2).
  ratio_ind <- mean(fit_ind$sd.dif2 / pmax(fit_ind$sd.dif, .Machine$double.xmin))
  ratio_con <- mean(fit_con$sd.dif2 / pmax(fit_con$sd.dif, .Machine$double.xmin))
  expect_gte(ratio_ind, 1 - 1e-8)
  expect_gt(ratio_con, 1.02)
  expect_gt(ratio_con, ratio_ind)

  # Shared random columns cancel in the difference on a shared grid.
  n2 <- length(dat$x1)
  B <- cbind(fit_con$matlist$X, fit_con$matlist$Z)
  Dmat <- B[seq_len(n2), , drop = FALSE] - B[(n2 + 1L):nrow(B), , drop = FALSE]
  n_shared_z <- ncol(fit_con$matlist$Z) / 2
  z_shared <- seq.int(ncol(fit_con$matlist$X) + 1L,
                      ncol(fit_con$matlist$X) + n_shared_z)
  expect_lt(max(abs(Dmat[, z_shared])), 1e-10)

  pred_d <- predict(fit_con, type = "diff", se.fit = TRUE)
  expect_equal(pred_d$fit, fit_con$eta.diff)
  expect_equal(pred_d$se.fit, fit_con$sd.dif)
  pred_s <- predict(fit_con, type = "shared", se.fit = TRUE)
  expect_equal(pred_s$fit, fit_con$eta.shared)
  expect_equal(pred_s$se.fit, fit_con$sd.shared)

  inf <- clgam_contrast_infer(fit_ind, uncond = FALSE)
  expect_equal(nrow(inf), n2)
  expect_true(all(inf$sign %in% c(-1L, 0L, 1L)))
  expect_equal(attr(inf, "structure"), "independent")
  expect_equal(attr(inf, "se_type"), "bayes")
  expect_equal(inf$se, as.numeric(fit_ind$sd.dif))
  expect_equal(inf$p.value, 2 * pnorm(-abs(inf$z)))

  inf_s <- clgam_contrast_infer(fit_con, type = "shared", uncond = FALSE)
  expect_equal(inf_s$se, as.numeric(fit_con$sd.shared))
})

test_that("clgam_contrast_infer(uncond=TRUE) inflates SEs and respects structure", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(n_coarse = 8L, n_fine_per = 4L, seed = 22L, spatial_amp = 0.5)
  set.seed(105)
  y1 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma))
  y2 <- as.numeric(dat$C %*% rpois(dat$n_fine, dat$gamma))
  fit <- clgam_contrast(
    y = c(y1, y2),
    x1 = c(dat$x1, dat$x1), x2 = c(dat$x2, dat$x2),
    exposure = c(dat$efine, dat$efine),
    C1 = dat$C, C2 = dat$C,
    knots = c(4L, 4L), elements = TRUE,
    structure = "contrast"
  )
  inf_b <- clgam_contrast_infer(fit, uncond = FALSE)
  inf_u <- clgam_contrast_infer(fit)
  expect_s3_class(fit, "clgam")
  expect_equal(attr(inf_u, "structure"), "contrast")
  expect_equal(attr(inf_u, "se_type"), "uncond")
  expect_true(all(inf_u$se >= inf_b$se - 1e-10))
  expect_gt(mean(inf_u$se / pmax(inf_b$se, .Machine$double.xmin)), 1)
  expect_error(clgam_contrast_infer(list()), "clgam fit")
})

# W1: spatial B3 indices vs dense N block. Does not change the Schur solve.

.kron_fixture <- function(n = 40L, ndx = c(6L, 7L), pord = c(2L, 2L),
                          n_nl_random = 0L, seed = 11L) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  MM1 <- mm_basis(x1, min(x1) - 0.01, max(x1) + 0.01,
                  ndx = ndx[1], bdeg = 3, pord = pord[1], decom = 1)
  MM2 <- mm_basis(x2, min(x2) - 0.01, max(x2) + 0.01,
                  ndx = ndx[2], bdeg = 3, pord = pord[2], decom = 1)
  meta <- clgam:::.sop_kron_meta(MM1, MM2, pord = pord, n_nl_random = n_nl_random)
  xz <- clgam:::.sop_spatial_XZ(MM1, MM2)
  list(MM1 = MM1, MM2 = MM2, meta = meta, xz = xz, x1 = x1, x2 = x2)
}

test_that("B1/B2/B3 sizes match rten2 columns and pois_SOP np[2:4]", {
  fx <- .kron_fixture()
  meta <- fx$meta
  xz <- fx$xz

  expect_equal(ncol(xz$Z_B1), meta$n_b1)
  expect_equal(ncol(xz$Z_B2), meta$n_b2)
  expect_equal(ncol(xz$Z_B3), meta$n_b3)
  expect_equal(ncol(xz$Z), meta$q_spatial)
  expect_equal(ncol(xz$X), meta$n_fixed_spatial)

  # pois_SOP: np <- c(prod(pord), (c2-pord2)*pord1, (c1-pord1)*pord2, ...)
  pord <- meta$pord
  np <- c(
    prod(pord),
    (meta$c2 - pord[2]) * pord[1],
    (meta$c1 - pord[1]) * pord[2],
    (meta$c1 - pord[1]) * (meta$c2 - pord[2])
  )
  expect_equal(unname(np[2:4]), c(meta$n_b1, meta$n_b2, meta$n_b3))
})

test_that("spatial index sets partition 1:q_spatial without overlap", {
  fx <- .kron_fixture()
  meta <- fx$meta
  all_idx <- c(meta$idx_B1, meta$idx_B2, meta$idx_B3)
  expect_equal(sort(all_idx), seq_len(meta$q_spatial))
  expect_equal(length(unique(all_idx)), meta$q_spatial)
  expect_equal(meta$idx_spatial, seq_len(meta$q_spatial))
  expect_length(intersect(meta$idx_B1, meta$idx_B2), 0)
  expect_length(intersect(meta$idx_B1, meta$idx_B3), 0)
  expect_length(intersect(meta$idx_B2, meta$idx_B3), 0)
})

test_that("G1inv.n / G2inv.n masks match B1=tau2, B2=tau1, B3=both", {
  fx <- .kron_fixture()
  meta <- fx$meta
  eps <- 1e-12
  expect_true(all(abs(meta$G1inv.n[meta$idx_B1]) < eps))
  expect_true(all(meta$G1inv.n[meta$idx_B2] > 0))
  expect_true(all(meta$G1inv.n[meta$idx_B3] > 0))
  expect_true(all(meta$G2inv.n[meta$idx_B1] > 0))
  expect_true(all(abs(meta$G2inv.n[meta$idx_B2]) < eps))
  expect_true(all(meta$G2inv.n[meta$idx_B3] > 0))

  la <- c(0.4, 1.7)
  ginv <- clgam:::.sop_ginv_spatial(meta, la)
  expect_equal(
    ginv,
    c((1 / la[2]) * meta$g2u,
      (1 / la[1]) * meta$g1u,
      (1 / la[2]) * meta$g2b + (1 / la[1]) * meta$g1b)
  )
  expect_equal(ginv[meta$idx_B1], (1 / la[2]) * meta$G2inv.n[meta$idx_B1])
  expect_equal(ginv[meta$idx_B2], (1 / la[1]) * meta$G1inv.n[meta$idx_B2])
})

test_that("N[idx_B3, idx_B3] equals the dense Gram of Z_B3 after eliminating X", {
  fx <- .kron_fixture(n = 50L)
  meta <- fx$meta
  X <- fx$xz$X
  Z <- fx$xz$Z
  Z3 <- fx$xz$Z_B3
  n <- nrow(X)
  w <- runif(n, 0.4, 2.2)
  sw <- sqrt(w)
  Xw <- X * sw
  Zw <- Z * sw
  Z3w <- Z3 * sw
  XtX <- crossprod(Xw)
  ZtX <- crossprod(Zw, Xw)
  ZtZ <- crossprod(Zw)
  A11inv <- solve(XtX)
  N <- ZtZ - ZtX %*% (A11inv %*% t(ZtX))

  bl <- clgam:::.sop_N_blocks(N, meta)
  N33_direct <- crossprod(Z3w) - crossprod(Z3w, Xw) %*% (A11inv %*% crossprod(Xw, Z3w))
  expect_equal(bl$N_33, N33_direct, tolerance = 1e-10)
  expect_equal(bl$N_ss, N, tolerance = 1e-12)
})

test_that("Schur cache N B3 block matches helper extraction (solve unchanged)", {
  fx <- .kron_fixture(n = 45L)
  meta <- fx$meta
  X <- fx$xz$X
  Z <- fx$xz$Z
  n <- nrow(X)
  w <- runif(n, 0.5, 1.8)
  sw <- sqrt(w)
  Xw <- X * sw
  Zw <- Z * sw
  XtX <- crossprod(Xw)
  ZtX <- crossprod(Zw, Xw)
  ZtZ <- crossprod(Zw)
  ZtXtZ <- rbind(t(ZtX), ZtZ)
  u <- rnorm(ncol(X) + ncol(Z))
  G <- clgam:::.sop_ginv_spatial(meta, c(0.8, 1.3))
  G <- 1 / G

  sop <- clgam:::.sop_solve_schur_R(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache = NULL)
  bl <- clgam:::.sop_N_blocks(sop$cache$N, meta)
  expect_equal(bl$N_33, sop$cache$N[meta$idx_B3, meta$idx_B3, drop = FALSE])
  expect_equal(dim(bl$N_33), c(meta$n_b3, meta$n_b3))
})

test_that("nonlinear random coefficients sit after B3", {
  n_nl <- 5L
  fx <- .kron_fixture(n_nl_random = n_nl)
  meta <- fx$meta
  expect_equal(meta$q, meta$q_spatial + n_nl)
  expect_equal(meta$idx_nl, meta$q_spatial + seq_len(n_nl))
  expect_length(intersect(meta$idx_nl, meta$idx_spatial), 0)
  expect_true(all(meta$G1inv.n[meta$idx_nl] == 0))
  expect_true(all(meta$G2inv.n[meta$idx_nl] == 0))
  expect_equal(max(meta$idx_B3), meta$q_spatial)
})

# --- W2: isolated B3 solve -------------------------------------------------

.kron_b3_system <- function(fx, la = c(0.8, 1.3), Z = NULL) {
  meta <- fx$meta
  X <- fx$xz$X
  if (is.null(Z)) Z <- fx$xz$Z
  n <- nrow(X)
  w <- runif(n, 0.5, 1.8)
  sw <- sqrt(w)
  Xw <- X * sw
  Zw <- Z * sw
  XtX <- crossprod(Xw)
  ZtX <- crossprod(Zw, Xw)
  ZtZ <- crossprod(Zw)
  ZtXtZ <- rbind(t(ZtX), ZtZ)
  u <- rnorm(ncol(X) + ncol(Z))
  list(
    XtX = XtX, ZtX = ZtX, ZtZ = ZtZ, ZtXtZ = ZtXtZ, u = u,
    la = la, w = w
  )
}

test_that("PCG Jacobi matches dense B3 solve (same estimator)", {
  fx <- .kron_fixture(n = 50L)
  meta <- fx$meta
  sys <- .kron_b3_system(fx)
  ginv <- clgam:::.sop_ginv_spatial(meta, sys$la)
  G <- 1 / ginv
  sop <- clgam:::.sop_solve_schur_R(
    sys$XtX, sys$ZtX, sys$ZtZ, sys$ZtXtZ, sys$u, G, cache = NULL
  )
  N33 <- sop$cache$N[meta$idx_B3, meta$idx_B3, drop = FALSE]
  rhs <- sop$cache$rhs2[meta$idx_B3]
  g3 <- ginv[meta$idx_B3]
  dens <- clgam:::.sop_solve_b3(N33, g3, rhs, method = "dense")
  pcg <- clgam:::.sop_solve_b3(N33, g3, rhs, method = "pcg", tol = 1e-12)
  expect_true(pcg$converged)
  expect_equal(pcg$theta, dens$theta, tolerance = 1e-8)
})

test_that("B3-only model: dense B3 solve matches full Schur b.random", {
  fx <- .kron_fixture(n = 48L)
  meta <- fx$meta
  sys <- .kron_b3_system(fx, Z = fx$xz$Z_B3)
  g3 <- clgam:::.sop_ginv_b3(meta, sys$la)
  G <- 1 / g3
  sop <- clgam:::.sop_solve_schur_R(
    sys$XtX, sys$ZtX, sys$ZtZ, sys$ZtXtZ, sys$u, G, cache = NULL
  )
  dens <- clgam:::.sop_solve_b3(
    sop$cache$N, g3, sop$cache$rhs2, method = "dense"
  )
  pcg <- clgam:::.sop_solve_b3(
    sop$cache$N, g3, sop$cache$rhs2, method = "pcg", tol = 1e-12
  )
  expect_equal(dens$theta, sop$b.random, tolerance = 1e-10)
  expect_true(pcg$converged)
  expect_equal(pcg$theta, sop$b.random, tolerance = 1e-8)
})

test_that("diag B3 method is the Jacobi one-step (not the dense solve)", {
  fx <- .kron_fixture(n = 40L)
  meta <- fx$meta
  sys <- .kron_b3_system(fx, Z = fx$xz$Z_B3)
  g3 <- clgam:::.sop_ginv_b3(meta, sys$la)
  sop <- clgam:::.sop_solve_schur_R(
    sys$XtX, sys$ZtX, sys$ZtZ, sys$ZtXtZ, sys$u, 1 / g3, cache = NULL
  )
  dens <- clgam:::.sop_solve_b3(sop$cache$N, g3, sop$cache$rhs2, method = "dense")
  diag_sol <- clgam:::.sop_solve_b3(sop$cache$N, g3, sop$cache$rhs2, method = "diag")
  Mdiag <- diag(sop$cache$N) + g3
  expect_equal(diag_sol$theta, sop$cache$rhs2 / Mdiag, tolerance = 1e-12)
  expect_gt(max(abs(diag_sol$theta - dens$theta)), 1e-6)
})

test_that("sparse_P2D + diagonal weights matches dense solve (C = I diagnostic)", {
  m1 <- 8L
  m2 <- 9L
  set.seed(3)
  w <- runif(m1 * m2, 0.4, 2)
  b <- rnorm(m1 * m2)
  lam1 <- 1 / 0.5
  lam2 <- 1 / 1.2
  th <- clgam:::.sop_solve_P2D_diag(m1, m2, lam1, lam2, w, b, pord = 2L)
  P <- as.matrix(sparse_P2D(m1, m2, lam1, lam2, pord = 2L))
  th_d <- as.numeric(solve(P + diag(w), b))
  expect_equal(th, th_d, tolerance = 1e-8)
})

# --- W3: exact B3-vs-rest inverse of S (same estimator) --------------------

test_that("Sinv kron_hybrid matches dense solve(S)", {
  fx <- .kron_fixture(n = 50L)
  meta <- fx$meta
  sys <- .kron_b3_system(fx)
  G <- 1 / clgam:::.sop_ginv_spatial(meta, sys$la)
  sop <- clgam:::.sop_solve_schur_R(
    sys$XtX, sys$ZtX, sys$ZtZ, sys$ZtXtZ, sys$u, G, cache = NULL
  )
  Sinv_d <- clgam:::.sop_Sinv_dense(sop$cache$N, G)
  Sinv_k <- clgam:::.sop_Sinv_kron(sop$cache$N, G, meta)
  expect_equal(Sinv_k, Sinv_d, tolerance = 1e-8)
})

test_that("kron_hybrid Schur matches dense b.fixed / b.random / dZtNZ", {
  fx <- .kron_fixture(n = 48L)
  meta <- fx$meta
  sys <- .kron_b3_system(fx)
  G <- 1 / clgam:::.sop_ginv_spatial(meta, sys$la)
  old <- getOption("clgam.sop.backend")
  on.exit(options(clgam.sop.backend = old), add = TRUE)

  options(clgam.sop.backend = "dense")
  dens <- clgam:::.sop_solve_schur_R(
    sys$XtX, sys$ZtX, sys$ZtZ, sys$ZtXtZ, sys$u, G,
    cache = NULL, meta = meta
  )
  options(clgam.sop.backend = "kron_hybrid")
  kron <- clgam:::.sop_solve_schur_R(
    sys$XtX, sys$ZtX, sys$ZtZ, sys$ZtXtZ, sys$u, G,
    cache = dens$cache, meta = meta
  )
  expect_equal(kron$b.fixed, dens$b.fixed, tolerance = 1e-10)
  expect_equal(kron$b.random, dens$b.random, tolerance = 1e-10)
  expect_equal(kron$dZtNZ, dens$dZtNZ, tolerance = 1e-10)
})

test_that("kron_hybrid with nonlinear random effects still matches dense", {
  n_nl <- 4L
  fx <- .kron_fixture(n = 45L, n_nl_random = n_nl)
  meta <- fx$meta
  n <- nrow(fx$xz$X)
  set.seed(19)
  fx$xz$Z <- cbind(fx$xz$Z, matrix(rnorm(n * n_nl), n, n_nl))
  sys <- .kron_b3_system(fx)
  G <- c(1 / clgam:::.sop_ginv_spatial(meta, sys$la), runif(n_nl, 0.4, 1.6))
  old <- getOption("clgam.sop.backend")
  on.exit(options(clgam.sop.backend = old), add = TRUE)

  options(clgam.sop.backend = "dense")
  dens <- clgam:::.sop_solve_schur_R(
    sys$XtX, sys$ZtX, sys$ZtZ, sys$ZtXtZ, sys$u, G,
    cache = NULL, meta = meta
  )
  options(clgam.sop.backend = "kron_hybrid")
  kron <- clgam:::.sop_solve_schur_R(
    sys$XtX, sys$ZtX, sys$ZtZ, sys$ZtXtZ, sys$u, G,
    cache = dens$cache, meta = meta
  )
  expect_equal(kron$b.fixed, dens$b.fixed, tolerance = 1e-10)
  expect_equal(kron$b.random, dens$b.random, tolerance = 1e-10)
  expect_equal(kron$dZtNZ, dens$dZtNZ, tolerance = 1e-10)
})

test_that("default SOP backend remains dense", {
  expect_identical(getOption("clgam.sop.backend", "dense"), "dense")
})

test_that("clgam Case A: kron_hybrid matches dense eta / tau / AIC", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 16L, n_fine_per = 5L, seed = 13L,
    spatial_amp = 0.6
  )
  args <- list(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, C = dat$C,
    exposure = dat$efine, knots = c(6L, 6L),
    elements = TRUE, trace = FALSE
  )
  old <- getOption("clgam.sop.backend")
  on.exit(options(clgam.sop.backend = old), add = TRUE)

  options(clgam.sop.backend = "dense")
  fit_d <- do.call(clgam, args)
  options(clgam.sop.backend = "kron_hybrid")
  fit_k <- do.call(clgam, args)

  expect_lt(max(abs(fit_k$eta - fit_d$eta)), 1e-6)
  expect_equal(fit_k$var.comp, fit_d$var.comp, tolerance = 1e-6)
  expect_equal(fit_k$aic, fit_d$aic, tolerance = 1e-5)
})

test_that("clgam with a smooth covariate: kron_hybrid matches dense", {
  skip_if_not_installed("sf")
  dat <- simulate_ata(
    n_coarse = 16L, n_fine_per = 5L, seed = 15L,
    include_covariate = TRUE, nl_amp = 1.0, spatial_amp = 0.5,
    nl_fun = function(z) sin(2 * pi * z)
  )
  args <- list(
    y = dat$y, x1 = dat$x1, x2 = dat$x2, C = dat$C,
    exposure = dat$efine, smooth = dat$nlcovfine,
    knots = c(6L, 6L), knots_nl = 8L,
    elements = TRUE, trace = FALSE
  )
  old <- getOption("clgam.sop.backend")
  on.exit(options(clgam.sop.backend = old), add = TRUE)

  options(clgam.sop.backend = "dense")
  fit_d <- do.call(clgam, args)
  options(clgam.sop.backend = "kron_hybrid")
  fit_k <- do.call(clgam, args)

  expect_lt(max(abs(fit_k$eta - fit_d$eta)), 1e-6)
  expect_equal(fit_k$var.comp, fit_d$var.comp, tolerance = 1e-6)
  expect_equal(fit_k$aic, fit_d$aic, tolerance = 1e-5)
})

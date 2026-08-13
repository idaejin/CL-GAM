#' Spatial SOP block indices and Kronecker-Schur prototype (W1--W3)
#'
#' Index sets for the three spatial random blocks in \code{pois_SOP()}:
#' \code{Z_s = cbind(rten2(Z2, X1), rten2(X2, Z1), rten2(Z2, Z1))}.
#' W2 adds an isolated B3 solver (\code{.sop_solve_b3}); it does not replace
#' the default dense Schur path. W3 inverts \code{S} by an exact B3-vs-rest
#' block factorization when \code{options(clgam.sop.backend = "kron_hybrid")}.
#'
#' Layout (random coefficients only, after dropping the \code{np[1]} fixed
#' block) must stay in lockstep with \code{pois_SOP}:
#' \itemize{
#'   \item B1: \code{Z2} box \code{X1}, size \code{(c2 - pord[2]) * pord[1]},
#'     penalized by \eqn{\tau_2} only (\code{g2u}).
#'   \item B2: \code{X2} box \code{Z1}, size \code{(c1 - pord[1]) * pord[2]},
#'     penalized by \eqn{\tau_1} only (\code{g1u}).
#'   \item B3: \code{Z2} box \code{Z1}, size
#'     \code{(c1 - pord[1]) * (c2 - pord[2])}, penalized by both
#'     (\code{g1b + g2b}).
#' }
#' Nonlinear covariate random effects, if any, are appended after B3.
#'
#' @name sop_kron_spatial
#' @keywords internal
NULL

#' Build Kronecker-Schur metadata from \code{mm_basis()} margins
#'
#' @param MM1,MM2 output of \code{\link{mm_basis}} for \code{x1} and \code{x2}
#' @param pord length-2 difference orders (same as \code{pois_SOP})
#' @param n_nl_random number of nonlinear-covariate random coefficients
#'   appended after the spatial block (0 for spatial-only Case A)
#' @return list with block sizes, 1-based index vectors, and the spatial
#'   G-inverse masks used by \code{pois_SOP}
#' @keywords internal
.sop_kron_meta <- function(MM1, MM2, pord = c(2L, 2L), n_nl_random = 0L) {
  pord <- as.integer(pord)
  if (length(pord) == 1L) {
    pord <- c(pord, pord)
  }
  if (length(pord) != 2L) {
    stop("pord must have length 1 or 2.", call. = FALSE)
  }
  c1 <- as.integer(MM1$m)
  c2 <- as.integer(MM2$m)
  d1 <- as.numeric(MM1$d)
  d2 <- as.numeric(MM2$d)
  n_nl_random <- as.integer(n_nl_random)
  if (n_nl_random < 0L) {
    stop("n_nl_random must be >= 0.", call. = FALSE)
  }
  if (length(d1) != (c1 - pord[1L]) || length(d2) != (c2 - pord[2L])) {
    stop(
      "mm_basis non-null eigenvalues do not match (c - pord); ",
      "check decom/pord against pois_SOP.",
      call. = FALSE
    )
  }

  n_fixed_spatial <- pord[1L] * pord[2L]
  n_b1 <- (c2 - pord[2L]) * pord[1L]
  n_b2 <- (c1 - pord[1L]) * pord[2L]
  n_b3 <- (c1 - pord[1L]) * (c2 - pord[2L])
  q_spatial <- n_b1 + n_b2 + n_b3

  idx_B1 <- seq_len(n_b1)
  idx_B2 <- n_b1 + seq_len(n_b2)
  idx_B3 <- n_b1 + n_b2 + seq_len(n_b3)
  idx_nl <- if (n_nl_random > 0L) {
    q_spatial + seq_len(n_nl_random)
  } else {
    integer(0L)
  }

  # Same construction as pois_SOP (G1inv.n / G2inv.n / Ginv).
  g1u <- rep(d1, times = pord[2L])
  g2u <- rep(d2, each = pord[1L])
  g1b <- rep(d1, times = (c2 - pord[2L]))
  g2b <- rep(d2, each = (c1 - pord[1L]))

  G1inv.n <- c(rep(0, n_b1), g1u, g1b, rep(0, n_nl_random))
  G2inv.n <- c(g2u, rep(0, n_b2), g2b, rep(0, n_nl_random))

  list(
    c1 = c1,
    c2 = c2,
    pord = pord,
    n_fixed_spatial = n_fixed_spatial,
    n_b1 = n_b1,
    n_b2 = n_b2,
    n_b3 = n_b3,
    q_spatial = q_spatial,
    q = q_spatial + n_nl_random,
    idx_B1 = idx_B1,
    idx_B2 = idx_B2,
    idx_B3 = idx_B3,
    idx_spatial = seq_len(q_spatial),
    idx_nl = idx_nl,
    d1 = d1,
    d2 = d2,
    g1u = g1u,
    g2u = g2u,
    g1b = g1b,
    g2b = g2b,
    G1inv.n = G1inv.n,
    G2inv.n = G2inv.n
  )
}

#' Spatial G-inverse diagonal matching \code{pois_SOP} inner SOP
#'
#' \code{Ginv = c((1/la[2])*g2u, (1/la[1])*g1u, (1/la[2])*g2b + (1/la[1])*g1b)}.
#' Nonlinear blocks are not included; append them as in \code{pois_SOP}.
#'
#' @param meta output of \code{\link{.sop_kron_meta}}
#' @param la length-2 spatial variance components \code{c(tau1^2, tau2^2)}
#' @return numeric vector of length \code{meta$q_spatial}
#' @keywords internal
.sop_ginv_spatial <- function(meta, la) {
  la <- as.numeric(la)
  if (length(la) < 2L) {
    stop("la must contain at least two spatial variance components.", call. = FALSE)
  }
  c(
    (1 / la[2L]) * meta$g2u,
    (1 / la[1L]) * meta$g1u,
    (1 / la[2L]) * meta$g2b + (1 / la[1L]) * meta$g1b
  )
}

#' Extract spatial / B3 blocks of the cached Schur Gram \code{N}
#'
#' @param N \code{q x q} matrix from \code{.sop_solve_schur()$cache$N}
#' @param meta output of \code{\link{.sop_kron_meta}}
#' @return named list of submatrices
#' @keywords internal
.sop_N_blocks <- function(N, meta) {
  if (!is.matrix(N) || nrow(N) != ncol(N)) {
    stop("N must be a square matrix.", call. = FALSE)
  }
  if (nrow(N) < meta$q_spatial) {
    stop("N is smaller than the spatial random block.", call. = FALSE)
  }
  list(
    N_11 = N[meta$idx_B1, meta$idx_B1, drop = FALSE],
    N_22 = N[meta$idx_B2, meta$idx_B2, drop = FALSE],
    N_33 = N[meta$idx_B3, meta$idx_B3, drop = FALSE],
    N_ss = N[meta$idx_spatial, meta$idx_spatial, drop = FALSE]
  )
}

#' Form \eqn{S = N\,\mathrm{diag}(G)+I}
#' @keywords internal
.sop_form_S <- function(N, G) {
  S <- sweep(N, 2L, G, `*`)
  diag(S) <- diag(S) + 1
  S
}

#' Dense inverse of the Schur complement \code{S}
#' @keywords internal
.sop_Sinv_dense <- function(N, G) {
  S <- .sop_form_S(N, G)
  Sinv <- try(solve(S), silent = TRUE)
  if (inherits(Sinv, "try-error")) {
    Sinv <- MASS::ginv(S)
  }
  Sinv
}

#' Whether B3 block inverse applies to this \code{meta} / \code{q}
#' @keywords internal
.sop_kron_applicable <- function(meta, q) {
  if (is.null(meta) || is.null(meta$idx_B3)) {
    return(FALSE)
  }
  idx3 <- as.integer(meta$idx_B3)
  if (length(idx3) < 1L || max(idx3) > q || min(idx3) < 1L) {
    return(FALSE)
  }
  n_rest <- q - length(idx3)
  n_rest >= 1L && length(idx3) >= 1L
}

#' Exact block inverse of \code{S}, splitting B3 vs the remaining random effects
#'
#' Same matrix as \code{solve(N \%*\% diag(G) + I)}. B3 is inverted as the
#' Schur complement of the (typically small) B1/B2/nl block. Used by
#' \code{options(clgam.sop.backend = "kron_hybrid")}.
#' @keywords internal
.sop_Sinv_kron <- function(N, G, meta) {
  q <- length(G)
  if (!.sop_kron_applicable(meta, q)) {
    return(.sop_Sinv_dense(N, G))
  }
  idx3 <- as.integer(meta$idx_B3)
  idx_r <- seq_len(q)[-idx3]
  S <- .sop_form_S(N, G)
  P <- S[idx_r, idx_r, drop = FALSE]
  Q <- S[idx_r, idx3, drop = FALSE]
  Rmat <- S[idx3, idx_r, drop = FALSE]
  Tm <- S[idx3, idx3, drop = FALSE]
  Pinv <- try(solve(P), silent = TRUE)
  if (inherits(Pinv, "try-error")) {
    Pinv <- MASS::ginv(P)
  }
  RP <- Rmat %*% Pinv
  T_eff <- Tm - RP %*% Q
  U <- try(solve(T_eff), silent = TRUE)
  if (inherits(U, "try-error")) {
    U <- MASS::ginv(T_eff)
  }
  PinvQ <- Pinv %*% Q
  URP <- U %*% RP
  Sinv <- matrix(0, q, q)
  Sinv[idx3, idx3] <- U
  Sinv[idx3, idx_r] <- -URP
  Sinv[idx_r, idx3] <- -PinvQ %*% U
  Sinv[idx_r, idx_r] <- Pinv + PinvQ %*% URP
  Sinv
}

#' Spatial mixed-model matrices matching \code{pois_SOP} (no PIRLS)
#'
#' @param MM1,MM2 \code{\link{mm_basis}} margins
#' @return list with \code{X}, \code{Z}, and the three spatial random blocks
#' @keywords internal
.sop_spatial_XZ <- function(MM1, MM2) {
  X <- rten2(MM2$X, MM1$X)
  Z1b <- rten2(MM2$Z, MM1$X)
  Z2b <- rten2(MM2$X, MM1$Z)
  Z3b <- rten2(MM2$Z, MM1$Z)
  list(X = X, Z = cbind(Z1b, Z2b, Z3b), Z_B1 = Z1b, Z_B2 = Z2b, Z_B3 = Z3b)
}

#' B3 slice of the spatial G-inverse diagonal
#' @keywords internal
.sop_ginv_b3 <- function(meta, la) {
  .sop_ginv_spatial(meta, la)[meta$idx_B3]
}

#' Solve the isolated B3 SOP normal equations
#'
#' Solves \eqn{(N_{33}+G_3^{-1})\theta = \mathrm{rhs}} where \code{rhs} is the
#' Schur \code{rhs2} slice (the right-hand side for \code{b2} in
#' \code{.sop_solve_schur_R}). Then \eqn{\theta} is the B3 random-effect
#' subvector (\code{b.random[idx_B3]}), not \code{b2}.
#'
#' \describe{
#'   \item{\code{dense}}{Direct \code{solve()}; reference for the same estimator.}
#'   \item{\code{pcg}}{CG with Jacobi preconditioner
#'     \eqn{M=\mathrm{diag}(N_{33})+G_3^{-1}}. Same system; tight \code{tol}
#'     reproduces \code{dense}.}
#'   \item{\code{diag}}{One Jacobi step \eqn{\theta=\mathrm{rhs}/M}. Diagnostic
#'     only (changes the estimator unless \eqn{N_{33}} is diagonal).}
#' }
#'
#' Isolated B3 does not replace the production Schur solve: W3 uses
#' \code{.sop_Sinv_kron()} so that \code{dZtNZ} stays exact.
#'
#' @param N33 B3 Gram block (\code{n_b3 x n_b3})
#' @param ginv3 B3 G-inverse diagonal (length \code{n_b3})
#' @param rhs right-hand side (Schur \code{rhs2} on B3)
#' @param method \code{"dense"}, \code{"pcg"}, or \code{"diag"}
#' @param tol,maxit PCG stopping rule (relative residual)
#' @param x0 optional PCG start (default: Jacobi guess)
#' @return list with \code{theta}, \code{iter}, \code{converged}, \code{relres},
#'   \code{method}
#' @keywords internal
.sop_solve_b3 <- function(N33, ginv3, rhs,
                          method = c("dense", "pcg", "diag"),
                          tol = 1e-10, maxit = NULL, x0 = NULL) {
  method <- match.arg(method)
  rhs <- as.numeric(rhs)
  ginv3 <- as.numeric(ginv3)
  n <- length(rhs)
  if (!is.matrix(N33) || nrow(N33) != n || ncol(N33) != n) {
    stop("N33 must be an n x n matrix with n = length(rhs).", call. = FALSE)
  }
  if (length(ginv3) != n) {
    stop("ginv3 must have length n = length(rhs).", call. = FALSE)
  }
  Mdiag <- diag(N33) + ginv3
  Mdiag <- pmax(Mdiag, .Machine$double.eps)

  if (identical(method, "diag")) {
    theta <- rhs / Mdiag
    r <- as.numeric(N33 %*% theta) + ginv3 * theta - rhs
    bnorm <- sqrt(sum(rhs * rhs))
    relres <- if (bnorm > 0) sqrt(sum(r * r)) / bnorm else 0
    return(list(
      theta = theta, iter = 1L, converged = TRUE,
      relres = relres, method = method
    ))
  }

  if (identical(method, "dense")) {
    A <- N33
    diag(A) <- diag(A) + ginv3
    theta <- as.numeric(solve(A, rhs))
    return(list(
      theta = theta, iter = NA_integer_, converged = TRUE,
      relres = 0, method = method
    ))
  }

  .sop_pcg_jacobi(N33, ginv3, rhs, Mdiag = Mdiag,
                  tol = tol, maxit = maxit, x0 = x0)
}

#' Jacobi-preconditioned CG for \eqn{(N + \mathrm{diag}(g))\theta = \mathrm{rhs}}
#' @keywords internal
.sop_pcg_jacobi <- function(N, ginv, rhs, Mdiag = NULL,
                            tol = 1e-10, maxit = NULL, x0 = NULL) {
  n <- length(rhs)
  if (is.null(maxit)) {
    maxit <- as.integer(max(200L, 4L * n))
  }
  if (is.null(Mdiag)) {
    Mdiag <- diag(N) + as.numeric(ginv)
    Mdiag <- pmax(Mdiag, .Machine$double.eps)
  }
  ginv <- as.numeric(ginv)
  x <- if (is.null(x0)) rhs / Mdiag else as.numeric(x0)
  Ax <- as.numeric(N %*% x) + ginv * x
  r <- rhs - Ax
  z <- r / Mdiag
  p <- z
  rz <- sum(r * z)
  bnorm <- sqrt(sum(rhs * rhs))
  if (!is.finite(bnorm) || bnorm == 0) {
    return(list(
      theta = x, iter = 0L, converged = TRUE, relres = 0, method = "pcg"
    ))
  }
  relres <- sqrt(sum(r * r)) / bnorm
  if (relres < tol) {
    return(list(
      theta = x, iter = 0L, converged = TRUE, relres = relres, method = "pcg"
    ))
  }
  for (k in seq_len(maxit)) {
    Ap <- as.numeric(N %*% p) + ginv * p
    pAp <- sum(p * Ap)
    if (!is.finite(pAp) || pAp <= 0) {
      return(list(
        theta = x, iter = k, converged = FALSE, relres = relres, method = "pcg"
      ))
    }
    alpha <- rz / pAp
    x <- x + alpha * p
    r <- r - alpha * Ap
    relres <- sqrt(sum(r * r)) / bnorm
    if (relres < tol) {
      return(list(
        theta = x, iter = as.integer(k), converged = TRUE,
        relres = relres, method = "pcg"
      ))
    }
    z <- r / Mdiag
    rz_new <- sum(r * z)
    p <- z + (rz_new / rz) * p
    rz <- rz_new
  }
  list(
    theta = x, iter = as.integer(maxit), converged = FALSE,
    relres = relres, method = "pcg"
  )
}

#' Sparse Kronecker \eqn{P+W} solve on original tensor coefficients (diagnostic)
#'
#' Boer / LMMsolver path: \code{sparse_P2D()} + diagonal weights, size
#' \code{m1 * m2}. Not the SVD-SOP B3 system; use for \code{C = I} benches.
#'
#' @keywords internal
.sop_solve_P2D_diag <- function(m1, m2, lambda1, lambda2, w, rhs,
                                pord = 2L) {
  P <- sparse_P2D(m1, m2, lambda1 = lambda1, lambda2 = lambda2, pord = pord)
  sparse_chol_solve(P, w, rhs)
}

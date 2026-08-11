#' Penalized CLM via Camarda–Durbán latent working response
#'
#' Implements the iterative scheme of Camarda & Durbán
#' (\url{https://arxiv.org/abs/2412.04956}): redistribute aggregated counts to a
#' fine-scale working response
#' \eqn{\tilde y = \gamma \odot (C^\top (y \oslash \mu))}, then update
#' P-spline coefficients with Kronecker penalty
#' \eqn{P = \lambda_1 (I\otimes D_1'D_1)+\lambda_2(D_2'D_2\otimes I)} without
#' building the composite working design \eqn{\breve B = W^{-1}C\Gamma B}.
#'
#' Unlike \code{pois_SOP}, smoothing is controlled by fixed \code{lambda}
#' (not SOP variance components). For irregular spatial \code{C} (Madrid),
#' GLAM array identities for \eqn{C=C_2\otimes C_1} are not used; the latent
#' \eqn{\tilde y} step still applies.
#'
#' @param y coarse Poisson counts
#' @param x1,x2 fine-scale coordinates (length = ncol(C))
#' @param efine fine-scale offset / exposure (default 1)
#' @param C composition matrix (coarse × fine)
#' @param ndx,bdeg,pord P-spline settings (same spirit as \code{pois_SOP})
#' @param lambda length-2 smoothness weights for the Kronecker penalty
#' @param thr,maxit convergence on relative change in \eqn{\eta}
#' @param sparse.backend passed to \code{as_comp_C}
#' @param trace print iteration info
#' @return list with \code{eta}, \code{gamma}, \code{mu}, \code{alpha}, \code{dev}, …
#' @export
#' @references
#' Camarda, C. G. and Durbán, M. (2025).
#' Fast Estimation of the Composite Link Model for Multidimensional Grouped Counts.
#' \emph{arXiv}:2412.04956.
pois_PCLM_latent <- function(y, x1, x2, efine = NULL, C,
                             ndx = c(20, 20), bdeg = c(3, 3), pord = 2L,
                             lambda = c(1, 1),
                             x1lim = NULL, x2lim = NULL,
                             thr = 1e-6, maxit = 100L,
                             sparse.backend = "auto",
                             trace = FALSE) {
  start.all <- proc.time()[3]
  dimfine <- length(x1)
  if (is.null(efine)) efine <- rep(1, dimfine)
  if (length(efine) == 1L) efine <- rep(efine, dimfine)
  if (is.null(x1lim)) x1lim <- c(min(x1) - 0.01, max(x1) + 0.01)
  if (is.null(x2lim)) x2lim <- c(min(x2) - 0.01, max(x2) + 0.01)
  if (length(bdeg) == 1L) bdeg <- c(bdeg, bdeg)
  if (length(ndx) == 1L) ndx <- c(ndx, ndx)
  if (length(lambda) == 1L) lambda <- c(lambda, lambda)

  C <- .as_comp_C(C, backend = sparse.backend)
  C_groups <- if (.is_partition_C(C)) .partition_groups(C) else NULL

  B1 <- bbasis(x1, x1lim[1], x1lim[2], ndx[1], bdeg[1])$B
  B2 <- bbasis(x2, x2lim[1], x2lim[2], ndx[2], bdeg[2])$B
  # Row-wise tensor product at scattered fine locations (Ayma / Currie style)
  B <- rten2(B2, B1)
  c1 <- ncol(B1)
  c2 <- ncol(B2)
  P <- sparse_P2D(c1, c2, lambda1 = lambda[1], lambda2 = lambda[2],
                  pord = as.integer(pord), backend = "Matrix")

  alpha <- rep(0, ncol(B))
  eta <- as.numeric(B %*% alpha)
  gamma <- efine * exp(eta)
  mu <- .comp_agg(C, gamma, C_groups)

  for (i in seq_len(maxit)) {
    y_tilde <- gamma * .comp_tmul(C, y / mu, C_groups)
    # IWLS working vector for Poisson mean γ (Camarda–Durbán)
    z <- eta + (y_tilde - gamma) / pmax(gamma, 1e-12)
    w <- gamma
    # (B' W B + P) α = B' (W z)  with W = diag(γ)
    if (isTRUE(getOption("clgam.use_rcpp", TRUE)) &&
        exists("btWb_cpp", mode = "function")) {
      tmp <- btWb_cpp(B, w, z)
      BtWB <- tmp$BtWB
      rhs <- as.numeric(tmp$rhs)
    } else {
      sw <- sqrt(pmax(w, 0))
      BtWB <- crossprod(sw * B)
      rhs <- as.numeric(crossprod(B, w * z))
    }
    A <- BtWB + P
    alpha_new <- tryCatch(
      as.numeric(Matrix::solve(A, rhs)),
      error = function(e) as.numeric(Matrix::solve(Matrix::forceSymmetric(A), rhs))
    )
    eta_new <- as.numeric(B %*% alpha_new)
    gamma <- efine * exp(eta_new)
    mu <- .comp_agg(C, gamma, C_groups)
    tol <- sum((eta_new - eta)^2) / max(sum(eta_new^2), 1e-12)
    alpha <- alpha_new
    eta <- eta_new
    if (trace) {
      cat(sprintf("%3d  tol=%.3e  dev=%.3f\n", i, tol,
                  2 * sum(y * log(ifelse(y == 0, 1, y / mu)) - (y - mu))))
    }
    if (tol < thr) break
  }

  dev <- 2 * sum(y * log(ifelse(y == 0, 1, y / mu)) - (y - mu))
  elapsed <- proc.time()[3] - start.all
  cat("Number of iterations:", i, "\n")
  cat("Lambda:", lambda, "\n")
  cat("Convergence criterion value: ", tol, "\n")
  cat("Elapsed time of estimation procedure:", elapsed, "seconds\n")

  list(
    eta = eta,
    gamma = gamma,
    mu = mu,
    alpha = alpha,
    lambda = lambda,
    ndx = ndx,
    bdeg = bdeg,
    pord = pord,
    niter = i,
    tol = tol,
    elapsed.time = elapsed,
    dev = dev,
    method = "Camarda-Durban-latent",
    matlist = list(B = B, B1 = B1, B2 = B2, C = C, P = P)
  )
}

#' C %*% gamma (partition-aware)
#' @keywords internal
.comp_agg <- function(C, gamma, groups = NULL) {
  if (!is.null(groups)) {
    as.numeric(rowsum(gamma, group = factor(groups, levels = seq_len(nrow(C)))))
  } else {
    as.numeric(C %*% gamma)
  }
}

#' C' %*% v (partition-aware): fine gets coarse value of its cell
#' @keywords internal
.comp_tmul <- function(C, v, groups = NULL) {
  if (!is.null(groups)) {
    as.numeric(v[groups])
  } else if (inherits(C, "Matrix")) {
    as.numeric(Matrix::crossprod(C, v))
  } else if (inherits(C, "spam")) {
    as.numeric(t(C) %*% v)
  } else {
    as.numeric(crossprod(C, v))
  }
}

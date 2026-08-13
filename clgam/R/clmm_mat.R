#' Matrices involved in CLMM estimation (speed-tuned)
#'
#' Same algebra as the original Ayma implementation, but:
#' - uses sparse \code{C} (\pkg{Matrix} or \pkg{spam}) when beneficial,
#' - uses \code{rowsum} aggregation when \code{C} is a 0-1 partition,
#' - avoids redundant \code{1/w} scalings.
#'
#' @param C composition matrix (dense, \pkg{Matrix}, or \pkg{spam})
#' @param gamma fine-scale mean intensity
#' @param X fixed-effects design (fine scale)
#' @param Z random-effects design (fine scale)
#' @param z working response (coarse)
#' @param w working weights (= \code{mu})
#' @param groups optional partition group ids (from \code{.partition_groups})
#'
#' @return list of cross-product blocks
#' @export
clmm_mat <- function(C, gamma, X, Z, z, w, groups = NULL) {
  if (!is.null(groups) &&
      isTRUE(getOption("clgam.use_rcpp", TRUE)) &&
      exists("comp_mul_groups_cpp", mode = "function")) {
    CGX <- comp_mul_groups_cpp(gamma, X, as.integer(groups), nrow(C))
    CGZ <- comp_mul_groups_cpp(gamma, Z, as.integer(groups), nrow(C))
  } else {
    CGX <- .comp_mul(C, gamma, X, groups = groups)
    CGZ <- .comp_mul(C, gamma, Z, groups = groups)
  }

  if (isTRUE(getOption("clgam.use_rcpp", TRUE)) &&
      exists("clmm_crossprod_cpp", mode = "function")) {
    return(clmm_crossprod_cpp(CGX, CGZ, z, w))
  }

  sw <- sqrt(1 / w)
  XtX <- crossprod(sw * CGX)
  ZtZ <- crossprod(sw * CGZ)
  XtZ <- crossprod(CGX, CGZ / w)
  ZtX <- t(XtZ)
  Xty <- crossprod(CGX, z)
  Zty <- crossprod(CGZ, z)
  yty <- sum((z^2) * w)
  ZtXtZ <- rbind(XtZ, ZtZ)
  u <- c(Xty, Zty)

  list(
    XtX = XtX, XtZ = XtZ, ZtX = ZtX, ZtZ = ZtZ,
    Xty = Xty, Zty = Zty, yty = yty, ZtXtZ = ZtXtZ, u = u
  )
}

#' Diagonal of A \%*\% S \%*\% t(A) without forming the full Gram matrix
#'
#' Row-wise quadratic forms used for pointwise SEs.
#' @keywords internal
.quad_diag <- function(A, S) {
  rowSums((A %*% S) * A)
}

#' Working response for composite-link Poisson PIRLS
#' @keywords internal
.clmm_working_z <- function(C, gamma, eta, mu, y, groups = NULL) {
  Cg_eta <- if (!is.null(groups)) {
    n_coarse <- nrow(C)
    as.numeric(rowsum(gamma * eta, group = factor(groups, levels = seq_len(n_coarse))))
  } else {
    as.numeric(C %*% (gamma * eta))
  }
  Cg_eta / mu + (y - mu) / mu
}

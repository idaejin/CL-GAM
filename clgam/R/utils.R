# Shared helpers for CL-GAMM fitting (internal)

#' Default coordinate limits (small pad)
#' @keywords internal
.clgam_xlim <- function(x) {
  c(min(x) - 0.01, max(x) + 0.01)
}

#' Expand scalar / NULL exposure to length n
#' @keywords internal
.clgam_exposure <- function(efine, n) {
  if (is.null(efine)) {
    rep(1, n)
  } else if (length(efine) == 1L) {
    rep(efine, n)
  } else {
    stopifnot(length(efine) == n)
    as.numeric(efine)
  }
}

#' Column-wise residualize A against span(B) (QR). Drop near-zero columns of B.
#' @keywords internal
.orth_cols <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  if (!nrow(A) || !ncol(A) || !ncol(B)) return(A)
  keep <- apply(B, 2L, function(v) stats::sd(v) > 1e-12)
  if (!any(keep)) return(A)
  B <- B[, keep, drop = FALSE]
  qr.resid(qr(B), A)
}

#' Pointwise SE for an additive smooth under a sum-to-zero (centering) constraint.
#'
#' For design columns \code{Cg} and coefficient covariance \code{Sg}, returns
#' \(\mathrm{se}_i=\sqrt{\mathrm{Var}(g_i-\bar g)}\), matching centred plots of \(g\).
#' @keywords internal
.se_smooth_centered <- function(Cg, Sg) {
  Cg <- as.matrix(Cg)
  Sg <- as.matrix(Sg)
  if (!nrow(Cg) || !ncol(Cg)) {
    return(rep(NA_real_, nrow(Cg)))
  }
  Cg_c <- sweep(Cg, 2L, colMeans(Cg), "-")
  sqrt(pmax(rowSums((Cg_c %*% Sg) * Cg_c), 0))
}

#' Attach S3 class and optional quiet fit banner
#' @keywords internal
.as_clgam <- function(x, call = NULL, family = c("spatial", "contrast")) {
  family <- match.arg(family)
  x$call <- call
  x$family <- family
  class(x) <- c(if (family == "contrast") "clgam_contrast" else NULL, "clgam", "list")
  x
}

#' Print fit summary lines (honours \code{trace})
#' @keywords internal
.clgam_report <- function(trace, niter, la, tol, elapsed) {
  if (!isTRUE(trace)) return(invisible(NULL))
  cat("Number of iterations:", niter, "\n")
  cat("Optimal variance components:", la, "\n")
  cat("Convergence criterion value:", tol, "\n")
  cat("Elapsed time of estimation procedure:", elapsed, "seconds\n")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

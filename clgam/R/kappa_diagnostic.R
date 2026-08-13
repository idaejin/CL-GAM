#' Functional-space overlap \eqn{\kappa}
#'
#' Geometric overlap between a spatial field \eqn{f} and the fine-scale
#' covariate space \eqn{A_{\mathrm{f}}},
#' \deqn{\kappa = \|P_{A_{\mathrm{f}}} f\|^2 / \|f\|^2.}
#' This is the diagnostic in the CL-GAMM identifiability section (manuscript
#' eq. 21): the fraction of spatial signal that also lies in the fine
#' B-spline / linear covariate space. Marginal correlation of \eqn{z} with
#' \eqn{f} does not determine \eqn{\kappa}.
#'
#' The projection uses the same QR residualization as \code{orth.smooth}
#' (\code{\link{.orth_cols}}). By default \eqn{f} is centred first, matching
#' the Monte Carlo construction.
#'
#' For a \code{clgam} fit, the default \eqn{f} is the fitted spatial field
#' \eqn{\hat\eta - \hat g(z) - \hat x'\beta_{\mathrm{lin}}}. That is an
#' \emph{empirical} overlap of the estimated spatial component with
#' \eqn{A_{\mathrm{f}}}; after \code{orth.smooth=TRUE} it is not forced to 0
#' (the spatial intercept is not residualized). The paper's design
#' \eqn{\kappa} uses an unrestricted field or known \eqn{f_{\mathrm{raw}}}
#' --- pass \code{f} explicitly. Fits also store \code{orth.info$kappa_by},
#' one \eqn{\kappa_k} per fine-scale smooth / linear covariate; see
#' \code{\link{summary.clgam}}.
#'
#' @param object a numeric spatial field, or a \code{clgam} fit
#' @param ... passed to methods
#' @return a length-1 numeric in \eqn{[0, 1]}
#' @seealso \code{\link{clgam}}, \code{\link{pois_SOP}}
#' @export
#' @examples
#' set.seed(1)
#' A <- cbind(runif(40), runif(40))
#' f <- as.numeric(A %*% c(0.8, -0.3)) + rnorm(40, sd = 0.05)
#' kappa_diagnostic(f, A)
kappa_diagnostic <- function(object, ...) {
  UseMethod("kappa_diagnostic")
}

#' @rdname kappa_diagnostic
#' @param A fine-covariate design whose column span is \eqn{A_{\mathrm{f}}}
#' @param center if \code{TRUE} (default), centre \code{object} before the ratio
#' @export
kappa_diagnostic.default <- function(object, A, center = TRUE, ...) {
  if (missing(A) || is.null(A)) {
    stop("Provide the fine-covariate design A (column span of A_f).", call. = FALSE)
  }
  .kappa_overlap(object, A, center = center)
}

#' @rdname kappa_diagnostic
#' @param f optional spatial field (defaults to the fitted spatial component)
#' @export
kappa_diagnostic.clgam <- function(object, f = NULL, A = NULL, center = TRUE, ...) {
  if (is.null(A)) {
    A <- .clgam_Af(object)
  }
  if (is.null(A) || !ncol(as.matrix(A))) {
    stop("Fit has no fine-scale covariate space A_f; kappa is undefined.",
         call. = FALSE)
  }
  if (is.null(f)) {
    f <- .clgam_spatial_field(object)
  }
  .kappa_overlap(f, A, center = center)
}

#' \eqn{\kappa = \|P_A f\|^2 / \|f\|^2} via the same QR as \code{.orth_cols}
#' @keywords internal
.kappa_overlap <- function(f, A, center = TRUE) {
  f <- as.numeric(f)
  if (!length(f) || !all(is.finite(f))) {
    stop("f must be a finite numeric vector.", call. = FALSE)
  }
  if (isTRUE(center)) {
    f <- f - mean(f)
  }
  den <- sum(f * f)
  if (!is.finite(den) || den <= 1e-12) {
    return(0)
  }
  A <- as.matrix(A)
  if (nrow(A) != length(f)) {
    stop("nrow(A) must equal length(f).", call. = FALSE)
  }
  f_perp <- as.numeric(.orth_cols(cbind(f), A)[, 1L])
  sum((f - f_perp)^2) / den
}

#' Fitted spatial field: eta minus covariate effects
#' @keywords internal
.clgam_spatial_field <- function(fit) {
  if (is.null(fit$eta)) {
    stop("Fit does not store eta.", call. = FALSE)
  }
  .clgam_spatial_from_parts(fit$eta, fit$nleffects, fit$leffects)
}

#' @keywords internal
.clgam_spatial_from_parts <- function(eta, nleffects = NULL, leffects = NULL) {
  f <- as.numeric(eta)
  if (!is.null(nleffects)) {
    f <- f - as.numeric(rowSums(as.matrix(nleffects)))
  }
  if (!is.null(leffects)) {
    f <- f - as.numeric(rowSums(as.matrix(leffects)))
  }
  f
}

#' Rebuild univariate B-spline bases from a stored fit (same as \code{pois_SOP}).
#' @keywords internal
.clgam_rebuild_Bk <- function(fit) {
  nlcov <- fit$nlcovfine
  if (is.null(nlcov)) return(NULL)
  nlcov <- as.matrix(nlcov)
  nk <- ncol(nlcov)
  if (nk < 1L) return(NULL)
  ndxnl <- fit$ndxnl %||% 15
  bdegnl <- fit$bdegnl %||% 3
  pordnl <- fit$pordnl %||% 2
  if (length(ndxnl) == 1L) ndxnl <- rep(ndxnl, nk)
  if (length(bdegnl) == 1L) bdegnl <- rep(bdegnl, nk)
  if (length(pordnl) == 1L) pordnl <- rep(pordnl, nk)
  decom <- fit$decom %||% 1
  nl.basis <- fit$nl.basis %||% "legacy"
  Bk <- vector("list", nk)
  for (k in seq_len(nk)) {
    xk <- nlcov[, k]
    xklim <- c(min(xk) - 0.01, max(xk) + 0.01)
    decom_k <- if (identical(nl.basis, "pspline")) 2L else decom
    Bk[[k]] <- mm_basis(
      x = xk, xl = xklim[1], xr = xklim[2], ndx = ndxnl[k],
      bdeg = bdegnl[k], pord = pordnl[k], decom = decom_k
    )$B
  }
  Bk
}

#' Rebuild \eqn{A_{\mathrm{f}}} from a \code{clgam} fit (raw fine B-splines + linear)
#' @keywords internal
.clgam_Af <- function(fit) {
  nlcov <- fit$nlcovfine
  nk <- if (!is.null(nlcov)) ncol(as.matrix(nlcov)) else 0L
  nl.level <- fit$nl.level
  Bk <- .clgam_rebuild_Bk(fit)
  if (nk > 0L && (is.null(nl.level) || !length(nl.level))) {
    nl.level <- rep("fine", nk)
  }
  .build_orth_Af(Bk = Bk, nl.level = nl.level, lcovfine = fit$lcovfine)
}

#' Named per-covariate designs for \eqn{\kappa_k} (fine smooths + linear).
#' @keywords internal
.clgam_Af_parts <- function(fit) {
  nlcov <- fit$nlcovfine
  nk <- if (!is.null(nlcov)) ncol(as.matrix(nlcov)) else 0L
  nl.level <- fit$nl.level
  if (nk > 0L && (is.null(nl.level) || !length(nl.level))) {
    nl.level <- rep("fine", nk)
  }
  .build_orth_Af_parts(
    Bk = .clgam_rebuild_Bk(fit),
    nl.level = nl.level,
    lcovfine = fit$lcovfine,
    smooth_names = .clgam_mat_colnames(nlcov, "s"),
    linear_names = .clgam_mat_colnames(fit$lcovfine, "x")
  )
}

#' Per-covariate overlap \eqn{\kappa_k = \|P_{A_k} f\|^2 / \|f\|^2}.
#' @keywords internal
.clgam_kappa_by <- function(fit, f = NULL) {
  stored <- fit$orth.info$kappa_by
  if (!is.null(stored) && length(stored)) return(stored)
  parts <- .clgam_Af_parts(fit)
  if (!length(parts)) return(numeric(0))
  if (is.null(f)) f <- .clgam_spatial_field(fit)
  vapply(parts, function(A) .kappa_overlap(f, A), numeric(1))
}

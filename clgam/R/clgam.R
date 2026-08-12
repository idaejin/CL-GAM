#' Fit a CL-GAMM (spatial / Case A--C)
#'
#' User-facing wrapper around \code{\link{pois_SOP}}.
#'
#' @param y coarse counts
#' @param x1,x2 fine-scale spatial coordinates (or pass \code{coords})
#' @param coords optional \code{n_fine x 2} matrix
#' @param C composition matrix (coarse \(\times\) fine)
#' @param exposure fine-scale offset / expected counts
#' @param linear optional fine-scale linear covariates
#' @param smooth optional covariates for univariate P-spline smooths
#'   (Case A/B/C). Matrix with one column (A or B) or two (Case C: fine then
#'   coarse-expanded). Default basis \code{nl.basis="pspline"}.
#' @param smooth_level optional \code{"fine"} / \code{"coarse"} labels for
#'   \code{smooth} columns (passed as \code{nl.level} to \code{pois_SOP}).
#'   Use \code{"coarse"} for Case B and \code{c("fine","coarse")} for Case C.
#' @param knots spatial knot counts (\code{ndx})
#' @param knots_nl knot count(s) for each smooth covariate (\code{ndxnl})
#' @param elements compute AIC/BIC/SEs (default \code{TRUE})
#' @param ... further arguments to \code{pois_SOP} (e.g. \code{nl.basis="legacy"}
#'   for archival SMiMR behaviour)
#' @return An object of class \code{"clgam"}.
#' @export
clgam <- function(y,
                  x1 = NULL,
                  x2 = NULL,
                  C,
                  exposure = NULL,
                  linear = NULL,
                  smooth = NULL,
                  smooth_level = NULL,
                  coords = NULL,
                  knots = c(15L, 15L),
                  knots_nl = 10L,
                  elements = TRUE,
                  ...) {
  if (!is.null(coords)) {
    coords <- as.matrix(coords)
    stopifnot(ncol(coords) == 2L)
    x1 <- coords[, 1L]
    x2 <- coords[, 2L]
  }
  stopifnot(!is.null(x1), !is.null(x2), !missing(C))
  dots <- list(...)
  if (is.null(dots$nl.basis)) dots$nl.basis <- "pspline"
  if (is.null(dots$nl.level) && !is.null(smooth_level)) {
    dots$nl.level <- smooth_level
  }
  out <- do.call(pois_SOP, c(list(
    y = y, x1 = x1, x2 = x2, efine = exposure,
    lcovfine = linear, nlcovfine = smooth,
    C = C, ndx = knots, ndxnl = knots_nl, elements = elements
  ), dots))
  out$call <- match.call()
  out
}

#' Fit a two-group CL-GAMM contrast (e.g. sex)
#'
#' @inheritParams clgam
#' @param C1,C2 composition matrices per group
#' @param group unused placeholder
#' @export
clgam_contrast <- function(y,
                           x1 = NULL,
                           x2 = NULL,
                           exposure = NULL,
                           group = NULL,
                           C1,
                           C2,
                           coords = NULL,
                           knots = c(15L, 15L),
                           elements = TRUE,
                           ...) {
  if (!is.null(coords)) {
    coords <- as.matrix(coords)
    stopifnot(ncol(coords) == 2L)
    x1 <- coords[, 1L]
    x2 <- coords[, 2L]
  }
  stopifnot(!is.null(x1), !is.null(x2), !missing(C1), !missing(C2))
  n_fine <- length(x1)
  stopifnot(n_fine %% 2L == 0L)
  n2 <- n_fine %/% 2L
  cat_fac <- factor(c(rep("g1", n2), rep("g2", n2)))
  out <- pois_incat_SOP(
    y = y, x1 = x1, x2 = x2, efine = exposure,
    cat = cat_fac, Ccat1 = C1, Ccat2 = C2,
    ndx = knots, elements = elements, ...
  )
  out$call <- match.call()
  out
}

#' @export
#' @method print clgam
print.clgam <- function(x, ...) {
  fam <- if (!is.null(x$family)) x$family else "spatial"
  cat("CL-GAMM fit (", fam, ")\n", sep = "")
  if (!is.null(x$call)) {
    cat("Call: ")
    print(x$call)
  }
  cat("Iterations:", x$niter, "  elapsed:", round(x$elapsed.time, 3), "s\n")
  if (isTRUE(x$diverged)) {
    cat("WARNING: variance-component updates diverged; fit did not converge.\n")
  }
  cat("Variance components (tau^2):", paste(signif(x$var.comp, 4), collapse = ", "), "\n")
  if (!is.null(x$aic)) {
    cat("AIC:", paste(round(x$aic, 3), collapse = ", "))
    if (!is.null(x$ed)) cat("  ed:", round(x$ed, 3))
    cat("\n")
  }
  invisible(x)
}

#' @export
#' @method summary clgam
summary.clgam <- function(object, ...) {
  out <- list(
    call = object$call,
    family = object$family,
    niter = object$niter,
    elapsed = object$elapsed.time,
    var.comp = object$var.comp,
    edf = object$edf,
    ed = object$ed,
    aic = object$aic,
    bic = object$bic,
    dev = object$dev
  )
  class(out) <- "summary.clgam"
  out
}

#' @export
#' @method print summary.clgam
print.summary.clgam <- function(x, ...) {
  cat("Summary of CL-GAMM fit\n")
  if (!is.null(x$call)) {
    cat("Call: ")
    print(x$call)
  }
  cat("Family:", x$family %||% "spatial", "\n")
  cat("Iterations:", x$niter, "  elapsed:", round(x$elapsed, 3), "s\n")
  cat("tau^2:", paste(signif(x$var.comp, 5), collapse = ", "), "\n")
  if (!is.null(x$edf)) cat("edf:", paste(round(x$edf, 3), collapse = ", "), "\n")
  if (!is.null(x$ed)) cat("ed:", round(x$ed, 3), "\n")
  if (!is.null(x$aic)) cat("AIC:", paste(round(x$aic, 3), collapse = ", "), "\n")
  if (!is.null(x$bic)) cat("BIC:", paste(round(x$bic, 3), collapse = ", "), "\n")
  invisible(x)
}

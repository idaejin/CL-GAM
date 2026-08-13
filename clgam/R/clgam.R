#' Fit a CL-GAMM (spatial / Case A--C)
#'
#' User-facing wrapper around \code{\link{pois_SOP}}. The first argument may
#' be coarse counts (positional API) or an \pkg{mgcv}-style formula
#' \code{y ~ s(x1, x2) + s(z)}. The formula is parsed by \pkg{clgam}; it is
#' not evaluated by \pkg{mgcv} (the symbol \code{s} is matched even if
#' \code{mgcv::s} is attached; see \code{\link{s}}). Left-hand side \code{y}
#' is \emph{coarse}; \code{s()} covariates are on the \emph{fine} support
#' (a length-\code{n_coarse} vector such as Case B \code{z_a} is expanded
#' through partition \code{C}). A bare name \code{z} is linear; \code{s(z)}
#' is a P-spline.
#'
#' @param y coarse counts, or a two-sided formula such as
#'   \code{y ~ s(x1, x2) + s(z)} or \code{y ~ s(x1, x2) + z} (bare names
#'   are linear)
#' @param x1,x2 fine-scale spatial coordinates (or pass \code{coords}).
#'   Ignored when \code{y} is a formula.
#' @param coords optional \code{n_fine x 2} matrix
#' @param C composition matrix (coarse \eqn{\times} fine). Optional when
#'   \code{data} is a \code{\link{simulate_ata}} list that already contains
#'   \code{C}.
#' @param exposure fine-scale offset / expected counts. Optional when
#'   \code{data} supplies \code{efine}.
#' @param linear optional fine-scale linear covariates. With a formula,
#'   a bare name \code{z} (or \code{I(z)}) is equivalent; if both are
#'   supplied, the formula terms are used.
#' @param smooth optional covariates for univariate P-spline smooths
#'   (Case A/B/C). Matrix with one column (A or B) or two (Case C: fine then
#'   coarse-expanded). Default basis \code{nl.basis="pspline"}. Ignored when
#'   \code{y} is a formula (\code{s()} terms are used instead).
#' @param smooth_level optional \code{"fine"} / \code{"coarse"} labels for
#'   \code{smooth} columns (passed as \code{nl.level} to \code{pois_SOP}).
#'   Use \code{"coarse"} for Case B and \code{c("fine","coarse")} for Case C.
#'   With a formula, use \code{s(z, level = "coarse")} instead.
#' @param knots spatial knot counts (\code{ndx}). Overridden by
#'   \code{s(x1, x2, ndx = )} (or \code{k = }) when \code{knots} is not passed.
#' @param knots_nl knot count(s) for each smooth covariate (\code{ndxnl}).
#'   Overridden by \code{s(z, ndx = )} when \code{knots_nl} is not passed.
#' @param elements compute AIC/BIC/SEs (default \code{TRUE})
#' @param family \code{stats::poisson()} (default) or
#'   \code{stats::quasipoisson()}. Under quasi-Poisson, point estimates match
#'   Poisson; Pearson \eqn{\hat\phi} is stored in \code{fit$phi} and SEs are
#'   scaled by \eqn{\sqrt{\hat\phi}}.
#' @param data optional list or data frame in which to find formula variables.
#'   A \code{\link{simulate_ata}} result also supplies \code{C} and
#'   \code{exposure} (\code{efine}) when those arguments are omitted.
#' @param orth.smooth identifiability projection; see \code{\link{pois_SOP}}.
#'   Default \code{NULL} follows the basis (\code{TRUE} for P-splines).
#'   Identified Cases A--C use \code{TRUE} with DGP \eqn{f_\perp}; pass
#'   \code{identifying="unrestricted"} in \code{\link{simulate_ata}} to
#'   reproduce the older tables. See \code{\link{simulate_ata_scenarios}}.
#' @param method \code{"SOP"} (default; PIRLS + SOP) or \code{"Laplace"}
#'   (TMB Laplace on the same mixed-model design; requires \pkg{TMB}).
#' @param re \code{NULL}/\code{"none"} or \code{"coarse"} (iid random
#'   effect per coarse count). See \code{\link{pois_SOP}}.
#' @param ... further arguments to \code{pois_SOP} or \code{pois_Laplace}
#'   (e.g. \code{nl.basis="legacy"} for archival SMiMR behaviour)
#' @return An object of class \code{"clgam"}. Formula fits store
#'   \code{fit$formula}. \code{fit$method} is \code{"SOP"} or
#'   \code{"Laplace"}; \code{fit$re_term} is \code{"none"} or
#'   \code{"coarse"}.
#' @seealso \code{\link{s}}, \code{\link{pois_SOP}}, \code{\link{pois_Laplace}}, \code{\link{simulate_ata}},
#'   \code{\link{simulate_ata_scenarios}}, \code{\link{kappa_diagnostic}},
#'   \code{\link{summary.clgam}}
#' @examples
#' \donttest{
#' dat <- simulate_ata(scenario = "A", n_coarse = 12, n_fine_per = 4, seed = 1)
#' fit <- clgam(
#'   y ~ s(x1, x2, ndx = c(8, 8), bdeg = 3, pord = 2) + s(z_f, ndx = 10),
#'   data = dat,
#'   orth.smooth = TRUE
#' )
#' # Same call with C and exposure named (user-facing interface):
#' # clgam(y ~ s(x1, x2) + s(z), C = C, exposure = ef, orth.smooth = TRUE)
#' }
#' @export
clgam <- function(y,
                  x1 = NULL,
                  x2 = NULL,
                  C = NULL,
                  exposure = NULL,
                  linear = NULL,
                  smooth = NULL,
                  smooth_level = NULL,
                  coords = NULL,
                  knots = c(15L, 15L),
                  knots_nl = 10L,
                  elements = TRUE,
                  family = stats::poisson(),
                  data = NULL,
                  orth.smooth = NULL,
                  method = c("SOP", "Laplace"),
                  re = NULL,
                  ...) {
  formula <- NULL
  parsed <- NULL
  if (inherits(y, "formula")) {
    parsed <- .clgam_parse_formula(y, data = data, env = parent.frame())
    formula <- parsed$formula
    y <- parsed$y
    x1 <- parsed$x1
    x2 <- parsed$x2
    if (!is.null(parsed$smooth_vars)) {
      if (!is.null(smooth)) {
        warning("Ignoring smooth=; formula s() terms are used.", call. = FALSE)
      }
      if (!is.null(smooth_level)) {
        warning(
          "Ignoring smooth_level=; formula s(..., level=) is used.",
          call. = FALSE
        )
      }
    }
    if (!is.null(parsed$linear_vars) && !is.null(linear)) {
      warning("Ignoring linear=; formula linear terms are used.", call. = FALSE)
    }
    if (!is.null(parsed$knots) && missing(knots)) {
      knots <- parsed$knots
    }
    if (!is.null(parsed$knots_nl) && missing(knots_nl)) {
      kn <- parsed$knots_nl
      kn[is.na(kn)] <- 10
      knots_nl <- kn
    }
  } else if (!is.null(coords)) {
    coords <- as.matrix(coords)
    stopifnot(ncol(coords) == 2L)
    x1 <- coords[, 1L]
    x2 <- coords[, 2L]
  }
  pulled <- .clgam_pull_design(data, C, exposure)
  C <- pulled$C
  exposure <- pulled$exposure
  if (is.null(C)) {
    stop(
      "C is required (composition matrix), or pass data= from simulate_ata().",
      call. = FALSE
    )
  }
  if (inherits(formula, "formula") && !is.null(parsed$smooth_vars)) {
    smooth <- .clgam_bind_smooth(parsed$smooth_vars, C)
    smooth_level <- parsed$smooth_level
  }
  if (inherits(formula, "formula") && !is.null(parsed$linear_vars)) {
    linear <- .clgam_bind_smooth(parsed$linear_vars, C)
  }
  stopifnot(!is.null(x1), !is.null(x2))
  dots <- list(...)
  if (is.null(dots$nl.basis)) dots$nl.basis <- "pspline"
  if (is.null(dots$nl.level) && !is.null(smooth_level)) {
    dots$nl.level <- smooth_level
  }
  if (inherits(formula, "formula") && !is.null(parsed)) {
    if (is.null(dots$bdeg) && !is.null(parsed$bdeg)) dots$bdeg <- parsed$bdeg
    if (is.null(dots$pord) && !is.null(parsed$pord)) dots$pord <- parsed$pord
    if (is.null(dots$bdegnl) && !is.null(parsed$bdegnl)) {
      dots$bdegnl <- parsed$bdegnl
    }
    if (is.null(dots$pordnl) && !is.null(parsed$pordnl)) {
      dots$pordnl <- parsed$pordnl
    }
  }
  if (!is.null(dots$family)) {
    stop("Pass family= as a named argument of clgam(), not via ...", call. = FALSE)
  }
  if (!is.null(dots$orth.smooth)) {
    stop("Pass orth.smooth= as a named argument of clgam(), not via ...",
         call. = FALSE)
  }
  if (!is.null(dots$method)) {
    stop("Pass method= as a named argument of clgam(), not via ...",
         call. = FALSE)
  }
  if (!is.null(dots$re)) {
    stop("Pass re= as a named argument of clgam(), not via ...", call. = FALSE)
  }
  method <- match.arg(method)
  args <- c(list(
    y = y, x1 = x1, x2 = x2, efine = exposure,
    lcovfine = linear, nlcovfine = smooth,
    C = C, ndx = knots, ndxnl = knots_nl, elements = elements,
    family = family, orth.smooth = orth.smooth, re = re
  ), dots)
  if (identical(method, "Laplace")) {
    keep <- intersect(names(args), names(formals(pois_Laplace)))
    out <- do.call(pois_Laplace, args[keep])
  } else {
    out <- do.call(pois_SOP, args)
  }
  out$call <- match.call()
  if (!is.null(formula)) out$formula <- formula
  out
}

#' Fit a two-group CL-GAMM contrast (e.g. sex)
#'
#' @inheritParams clgam
#' @param y stacked coarse counts for both groups
#' @param x1,x2 stacked fine-scale spatial coordinates (or pass \code{coords})
#' @param exposure stacked fine-scale offset / expected counts
#' @param knots spatial knot counts (\code{ndx})
#' @param C1,C2 composition matrices per group
#' @param group unused placeholder
#' @param structure \code{"independent"} (default; two spatial fields, as
#'   in the manuscript) or \code{"contrast"} (shared field plus
#'   difference: \eqn{\eta_1=f+d}, \eqn{\eta_2=f-d}).
#' @details
#' Pointwise inference uses Bayesian SEs conditional on
#' \eqn{\widehat{\tau}^2} (the same linear functional of the joint
#' mixed-model covariance in both structures). With \code{elements=TRUE}
#' the fit stores \code{sd.eta} (stacked predictors), \code{sd.dif} for
#' \eqn{\widehat{\eta}_1-\widehat{\eta}_2} (accounts for
#' \eqn{\mathrm{Cov}(\eta_1,\eta_2)}), \code{sd.shared} for
#' \eqn{(\widehat{\eta}_1+\widehat{\eta}_2)/2}, and \code{sd.dif2}
#' (\eqn{\sqrt{\mathrm{se}_1^2+\mathrm{se}_2^2}}, which ignores that
#' covariance). Use \code{\link{clgam_contrast_infer}(fit)} for the
#' difference surface (unconditional SEs by default). Under
#' \code{structure="contrast"}
#' and a shared fine grid, the shared random field cancels in the
#' difference, so \code{sd.dif} is the SE of \eqn{2\widehat{d}} plus the
#' fixed group contrast.
#' @seealso \code{\link{clgam_contrast_infer}}, \code{\link{pois_incat_SOP}}
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
                           family = stats::poisson(),
                           structure = c("independent", "contrast"),
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
  structure <- match.arg(structure)
  dots <- list(...)
  if (!is.null(dots$method) || !is.null(dots$re)) {
    stop(
      "clgam_contrast is PIRLS+SOP only; method= and re= are not supported.",
      call. = FALSE
    )
  }
  out <- pois_incat_SOP(
    y = y, x1 = x1, x2 = x2, efine = exposure,
    cat = cat_fac, Ccat1 = C1, Ccat2 = C2,
    ndx = knots, elements = elements, family = family,
    structure = structure, ...
  )
  out$call <- match.call()
  out
}

#' @export
#' @method print clgam
print.clgam <- function(x, ...) {
  type <- x$type %||% {
    if (identical(x$family, "contrast") || inherits(x, "clgam_contrast")) {
      "contrast"
    } else {
      "spatial"
    }
  }
  fam <- x$family %||% "poisson"
  if (fam %in% c("spatial", "contrast")) fam <- "poisson"
  cat("CL-GAMM fit (", type, ", family=", fam,
      ", method=", x$method %||% "SOP",
      if (!is.null(x$re_term) && !identical(x$re_term, "none"))
        paste0(", re=", x$re_term) else "",
      if (identical(type, "contrast") && !is.null(x$structure))
        paste0(", structure=", x$structure) else "",
      ")\n", sep = "")
  if (!is.null(x$call)) {
    cat("Call: ")
    print(x$call)
  }
  if (!is.null(x$formula)) {
    cat("Formula: ")
    print(x$formula)
  }
  conv <- x$converged
  if (isTRUE(x$diverged)) {
    cat("Converged: FALSE (variance-component update diverged)\n")
  } else if (isTRUE(conv)) {
    cat("Converged: TRUE  (outer PIRLS ", x$niter, ", ",
        round(x$elapsed.time, 3), " s)\n", sep = "")
  } else if (identical(conv, FALSE)) {
    cat("Converged: FALSE  (outer PIRLS ", x$niter, "/",
        if (!is.null(x$maxit)) x$maxit[1] else "?", ", ",
        round(x$elapsed.time, 3), " s)\n", sep = "")
  } else {
    cat("Iterations:", x$niter, "  elapsed:", round(x$elapsed.time, 3), "s\n")
  }
  vc <- x$var.comp
  if (!is.null(names(vc)) && any(nzchar(names(vc)))) {
    cat("tau^2:", paste(sprintf("%s=%s", names(vc), signif(vc, 4)), collapse = "  "), "\n")
  } else {
    cat("Variance components (tau^2):", paste(signif(vc, 4), collapse = ", "), "\n")
  }
  if (!is.null(x$phi)) {
    cat("phi (Pearson):", signif(x$phi, 4), "\n")
  }
  aic_show <- if (!is.null(x$aic_total)) x$aic_total else x$aic
  bic_show <- if (!is.null(x$bic_total)) x$bic_total else x$bic
  if (!is.null(aic_show) && length(aic_show) == 1L) {
    cat("AIC:", round(aic_show, 3))
    if (!is.null(bic_show) && length(bic_show) == 1L) cat("  BIC:", round(bic_show, 3))
    if (!is.null(x$ed)) cat("  ED:", round(x$ed, 3))
    cat("\n")
  }
  invisible(x)
}

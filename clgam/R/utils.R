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
  Q <- qr.Q(qr(B))
  A - Q %*% crossprod(Q, A)
}

#' Resolve per-smooth support labels for identifiability projection.
#' @param nk number of nonlinear smooth columns
#' @param nl.level \code{NULL} (all \code{"fine"}), or length-\code{nk} labels
#'   \code{"fine"} / \code{"coarse"} (recycled if length 1).
#' @keywords internal
.resolve_nl_level <- function(nk, nl.level = NULL) {
  nk <- as.integer(nk)
  if (nk <= 0L) return(character(0))
  if (is.null(nl.level)) {
    return(rep("fine", nk))
  }
  lvl <- tolower(as.character(nl.level))
  if (length(lvl) == 1L) lvl <- rep(lvl, nk)
  if (length(lvl) != nk) {
    stop("nl.level must be NULL or length 1 / ncol(nlcovfine).", call. = FALSE)
  }
  bad <- setdiff(unique(lvl), c("fine", "coarse"))
  if (length(bad)) {
    stop("nl.level entries must be 'fine' or 'coarse'.", call. = FALSE)
  }
  lvl
}

#' Column names for a covariate matrix; fill missing with \code{prefix} + index.
#' @keywords internal
.clgam_mat_colnames <- function(M, prefix = "s") {
  if (is.null(M)) return(character(0))
  M <- as.matrix(M)
  n <- ncol(M)
  if (!n) return(character(0))
  cn <- colnames(M)
  if (is.null(cn) || !length(cn)) cn <- paste0(prefix, seq_len(n))
  empty <- is.na(cn) | !nzchar(cn)
  if (any(empty)) cn[empty] <- paste0(prefix, which(empty))
  make.unique(cn, sep = ".")
}

#' Design space used to orthogonalize spatial bases: fine B-splines + linear fine.
#' Coarse-scale smooths (Case B/C aggregated terms) are excluded by design.
#' @keywords internal
.build_orth_Af <- function(Bk = NULL, nl.level = NULL, lcovfine = NULL) {
  parts <- .build_orth_Af_parts(Bk = Bk, nl.level = nl.level, lcovfine = lcovfine)
  if (!length(parts)) return(NULL)
  do.call(cbind, unname(parts))
}

#' Named per-covariate designs whose union is \eqn{A_{\mathrm{f}}}.
#' Fine smooths (raw B-splines) and linear columns; coarse smooths omitted.
#' @keywords internal
.build_orth_Af_parts <- function(Bk = NULL, nl.level = NULL, lcovfine = NULL,
                                 smooth_names = NULL, linear_names = NULL) {
  parts <- list()
  if (!is.null(Bk) && length(Bk)) {
    nk <- length(Bk)
    if (is.null(nl.level) || !length(nl.level)) nl.level <- rep("fine", nk)
    if (is.null(smooth_names) || length(smooth_names) != nk) {
      smooth_names <- paste0("s", seq_len(nk))
    }
    for (k in seq_len(nk)) {
      if (!isTRUE(nl.level[k] == "fine")) next
      parts[[smooth_names[k]]] <- as.matrix(Bk[[k]])
    }
  }
  if (!is.null(lcovfine)) {
    L <- as.matrix(lcovfine)
    nj <- ncol(L)
    if (is.null(linear_names) || length(linear_names) != nj) {
      linear_names <- .clgam_mat_colnames(L, "x")
    }
    for (j in seq_len(nj)) {
      parts[[linear_names[j]]] <- L[, j, drop = FALSE]
    }
  }
  parts
}

#' Max absolute cross-product ||Q' A||_max after residualizing A perp span(B).
#' @keywords internal
.orth_cross_max <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  if (!nrow(A) || !ncol(A) || !ncol(B)) return(0)
  keep <- apply(B, 2L, function(v) stats::sd(v) > 1e-12)
  if (!any(keep)) return(0)
  B <- B[, keep, drop = FALSE]
  Q <- qr.Q(qr(B))
  max(abs(crossprod(Q, A)))
}

#' Pointwise SE for an additive smooth under a sum-to-zero (centering) constraint.
#'
#' For design columns \code{Cg} and coefficient covariance \code{Sg}, returns
#' \eqn{\mathrm{se}_i=\sqrt{\mathrm{Var}(g_i-\bar g)}}, matching centred plots of \eqn{g}.
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

#' Resolve GLM family for CL-GAM (poisson / quasipoisson only).
#'
#' Accepts a \code{stats::family} object, a family generator, or a string.
#' Quasi-Poisson uses the same PIRLS+SOP point estimates as Poisson; a Pearson
#' dispersion \eqn{\hat\phi} is estimated post-hoc and standard errors are
#' inflated by \eqn{\sqrt{\hat\phi}} (see \code{\link{.pearson_phi}}).
#'
#' @keywords internal
.resolve_clgam_family <- function(family = stats::poisson()) {
  if (is.null(family)) {
    family <- stats::poisson()
  } else if (is.character(family)) {
    family <- tolower(family[[1L]])
    family <- switch(
      family,
      poisson = stats::poisson(),
      quasipoisson = stats::quasipoisson(),
      stop(
        "clgam family must be poisson() or quasipoisson() (got '",
        family, "').",
        call. = FALSE
      )
    )
  } else if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, "family")) {
    stop("clgam 'family' must be a stats::family object or name.", call. = FALSE)
  }
  fam <- tolower(family$family)
  if (!fam %in% c("poisson", "quasipoisson")) {
    stop(
      "clgam supports only poisson and quasipoisson (got '",
      family$family, "'). Negative binomial is not implemented.",
      call. = FALSE
    )
  }
  list(
    family = family,
    name = fam,
    quasi = identical(fam, "quasipoisson")
  )
}

#' Pearson dispersion for (quasi-)Poisson composite-link fits.
#'
#' \deqn{\hat\phi = \frac{1}{n-p}\sum_i (y_i-\hat\mu_i)^2 / \hat\mu_i}
#' with \eqn{p} the effective dimension (\code{ed} when available).
#'
#' @keywords internal
.pearson_phi <- function(y, mu, ed) {
  y <- as.numeric(y)
  mu <- pmax(as.numeric(mu), 1e-12)
  stopifnot(length(y) == length(mu))
  n <- length(y)
  p <- if (length(ed) != 1L || is.null(ed) || !is.finite(ed[[1L]])) {
    0
  } else {
    as.numeric(ed[[1L]])
  }
  df <- max(n - p, 1)
  sum((y - mu)^2 / mu) / df
}

#' Inflate stored SE fields by \eqn{\sqrt{\phi}} (quasi-Poisson).
#' @keywords internal
.scale_se_phi <- function(out, phi) {
  if (!is.finite(phi) || abs(phi - 1) < .Machine$double.eps) {
    return(out)
  }
  s <- sqrt(phi)
  for (nm in c(
    "sd.eta", "sd.exp.eta", "sd.dif", "sd.dif2", "sd.shared",
    "sdleffects", "sdnleffects"
  )) {
    if (!is.null(out[[nm]])) {
      out[[nm]] <- out[[nm]] * s
    }
  }
  out
}

#' PIRLS-weighted Bayesian blocks (M1, M2) for AIC/BIC and SEs at convergence.
#'
#' Called once at the final outer PIRLS step when \code{elements=TRUE}, then
#' reused instead of recomputing \code{clmm_mat} + \code{inv_bblock2}.
#' @keywords internal
.pirls_bayes_blocks <- function(C, gamma, eta, mu, y, X, Z, ginvsp,
                                groups = NULL) {
  z <- .clmm_working_z(C, gamma, eta, mu, y, groups = groups)
  opt.mat <- clmm_mat(C, gamma, X, Z, z, mu, groups = groups)
  ZtZpen <- opt.mat$ZtZ
  diag(ZtZpen) <- diag(ZtZpen) + ginvsp
  M1 <- inv_bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, ZtZpen)
  M2 <- bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, opt.mat$ZtZ)
  list(mat = opt.mat, M1 = M1, M2 = M2)
}

#' Diagonal penalty inverse vector ginvsp for spatial + optional nl smooths
#' and optional coarse RE.
#' @keywords internal
.ginvsp_from_la <- function(la, g2u, g1u, g2b, g1b, dk = NULL, n_re = 0L) {
  ginvsp <- c((1 / la[2]) * g2u, (1 / la[1]) * g1u,
              (1 / la[2]) * g2b + (1 / la[1]) * g1b)
  if (!is.null(dk) && length(dk)) {
    for (k in seq_along(dk)) {
      ginvsp <- c(ginvsp, (1 / la[k + 2L]) * dk[[k]])
    }
  }
  n_re <- as.integer(n_re)
  if (n_re > 0L) {
    ginvsp <- c(ginvsp, rep(1 / la[length(la)], n_re))
  }
  ginvsp
}

#' @keywords internal
.ginvsp_from_la_incat <- function(la, g2u, g1u, g2b, g1b) {
  c(
    (1 / la[2]) * g2u, (1 / la[1]) * g1u,
    (1 / la[2]) * g2b + (1 / la[1]) * g1b,
    (1 / la[4]) * g2u, (1 / la[3]) * g1u,
    (1 / la[4]) * g2b + (1 / la[3]) * g1b
  )
}

#' Attach S3 class; \code{family} is the GLM family, \code{type} is spatial/contrast.
#' @keywords internal
.as_clgam <- function(x,
                      call = NULL,
                      type = c("spatial", "contrast"),
                      family = "poisson") {
  type <- match.arg(type)
  # Back-compat: older callers passed family = "spatial"|"contrast"
  if (is.character(family) && length(family) == 1L &&
      family %in% c("spatial", "contrast")) {
    type <- family
    family <- "poisson"
  }
  fam <- .resolve_clgam_family(family)
  x$call <- call
  x$type <- type
  x$family <- fam$name
  if (is.null(x$phi)) x$phi <- 1
  class(x) <- c(if (type == "contrast") "clgam_contrast" else NULL, "clgam", "list")
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

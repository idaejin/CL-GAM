#' Pointwise inference for a two-group CL-GAMM contrast
#'
#' @description
#' Call as \code{clgam_contrast_infer(fit)} where \code{fit} is a
#' \code{\link{clgam}} object returned by \code{\link{clgam_contrast}}.
#' Defaults are the manuscript contrast: difference surface
#' \eqn{\widehat{\eta}_1-\widehat{\eta}_2} with unconditional SEs.
#'
#' The difference uses the joint covariance \code{sd.dif}, not
#' \code{sd.dif2}. Unconditional SEs add the Wood–Pya–Säfken correction
#' for uncertainty in \eqn{\log\tau^2}. Pass \code{uncond = FALSE} for
#' the faster Bayesian SEs already stored on \code{fit}.
#'
#' @param fit a \code{\link{clgam}} fit from \code{\link{clgam_contrast}}
#'   (\code{elements = TRUE})
#' @param type \code{"diff"} (default; \eqn{\eta_1-\eta_2}),
#'   \code{"shared"} (\eqn{(\eta_1+\eta_2)/2}), or \code{"eta"}
#' @param level pointwise interval level
#' @param uncond \code{TRUE} (default): inflate SEs for variance-component
#'   uncertainty. \code{FALSE} uses \code{fit$sd.dif} / \code{sd.shared} /
#'   \code{sd.eta} as stored.
#' @param eps step on the log-variance scale when \code{uncond = TRUE}
#' @param ... unused; for S3
#' @return A data frame with \code{estimate}, \code{se}, \code{z},
#'   \code{p.value}, \code{lower}, \code{upper}. For \code{type = "diff"}
#'   also \code{sign} in \eqn{\{-1,0,1\}} (interval below / covers /
#'   above zero). Attributes record \code{structure}, \code{level},
#'   \code{uncond}, and \code{se_type}.
#' @seealso \code{\link{clgam_contrast}}
#' @export
clgam_contrast_infer <- function(fit, ...) {
  UseMethod("clgam_contrast_infer")
}

#' @rdname clgam_contrast_infer
#' @export
#' @method clgam_contrast_infer clgam
clgam_contrast_infer.clgam <- function(fit,
                                       type = c("diff", "shared", "eta"),
                                       level = 0.95,
                                       uncond = TRUE,
                                       eps = 1e-3,
                                       ...) {
  type <- match.arg(type)
  if (!identical(fit$type, "contrast") &&
      !inherits(fit, "clgam_contrast") &&
      is.null(fit$eta.diff)) {
    stop(
      "clgam_contrast_infer(): 'fit' is a clgam object but not a two-group ",
      "contrast. Pass the result of clgam_contrast().",
      call. = FALSE
    )
  }
  if (is.null(fit$eta) || is.null(fit$eta.diff)) {
    stop("clgam_contrast_infer(): two-group contrast fit required.", call. = FALSE)
  }
  zcrit <- stats::qnorm((1 + level) / 2)
  est <- switch(type,
    diff = as.numeric(fit$eta.diff),
    shared = as.numeric(fit$eta.shared),
    eta = as.numeric(fit$eta)
  )
  se <- switch(type,
    diff = fit$sd.dif,
    shared = fit$sd.shared,
    eta = fit$sd.eta
  )
  if (is.null(se)) {
    stop(
      "clgam_contrast_infer(): SEs missing; refit with elements=TRUE.",
      call. = FALSE
    )
  }
  se <- as.numeric(se)
  se_type <- "bayes"
  if (isTRUE(uncond)) {
    se <- .clgam_contrast_uncond_se(fit, type = type, eps = eps)
    se_type <- "uncond"
  }
  z <- est / se
  z[!is.finite(z)] <- NA_real_
  p <- 2 * stats::pnorm(-abs(z))
  out <- data.frame(
    estimate = est,
    se = se,
    z = z,
    p.value = p,
    lower = est - zcrit * se,
    upper = est + zcrit * se,
    stringsAsFactors = FALSE
  )
  if (identical(type, "diff")) {
    out$sign <- as.integer(
      as.integer(sign(out$lower) == sign(out$upper)) * sign(out$estimate)
    )
    out$sign[out$se == 0 & out$estimate == 0] <- 0L
  }
  attr(out, "type") <- type
  attr(out, "structure") <- fit$structure %||% "independent"
  attr(out, "level") <- level
  attr(out, "uncond") <- isTRUE(uncond)
  attr(out, "se_type") <- se_type
  out
}

#' @rdname clgam_contrast_infer
#' @export
#' @method clgam_contrast_infer default
clgam_contrast_infer.default <- function(fit, ...) {
  stop(
    "clgam_contrast_infer() needs a clgam fit from clgam_contrast().",
    call. = FALSE
  )
}

#' Recover C1/C2 from a contrast fit (stored or split from bdiag C).
#' @keywords internal
.clgam_split_contrast_C <- function(fit) {
  if (!is.null(fit$matlist$C1) && !is.null(fit$matlist$C2)) {
    return(list(C1 = fit$matlist$C1, C2 = fit$matlist$C2))
  }
  C <- fit$matlist$C
  if (is.null(C)) {
    stop("clgam_contrast_infer(): fit does not store C; cannot compute uncond SEs.",
         call. = FALSE)
  }
  n2 <- length(as.numeric(fit$eta)) %/% 2L
  C2_cols <- C[, seq.int(n2 + 1L, ncol(C)), drop = FALSE]
  rs2 <- as.numeric(Matrix::rowSums(abs(C2_cols)))
  n_g1 <- match(TRUE, rs2 > 0) - 1L
  if (is.na(n_g1) || n_g1 < 1L) {
    stop("clgam_contrast_infer(): could not split stacked C into C1/C2.",
         call. = FALSE)
  }
  list(
    C1 = C[seq_len(n_g1), seq_len(n2), drop = FALSE],
    C2 = C[seq.int(n_g1 + 1L, nrow(C)), seq.int(n2 + 1L, ncol(C)), drop = FALSE]
  )
}

#' Wood–Pya–Säfken inflate of a contrast linear functional.
#' @keywords internal
.clgam_contrast_uncond_se <- function(fit, type, eps = 1e-3) {
  if (is.null(fit$var.comp) || is.null(fit$edf) ||
      is.null(fit$b.fixed) || is.null(fit$b.random)) {
    stop("uncond=TRUE needs var.comp, edf, and coefficients on the fit.",
         call. = FALSE)
  }
  la <- as.numeric(fit$var.comp)
  edf <- as.numeric(fit$edf)
  if (length(la) != length(edf)) {
    stop("uncond=TRUE: length(var.comp) != length(edf).", call. = FALSE)
  }
  bold <- c(as.numeric(fit$b.fixed), as.numeric(fit$b.random))
  n_fine <- length(as.numeric(fit$eta))
  n2 <- n_fine %/% 2L
  CC <- .clgam_split_contrast_C(fit)
  fam_arg <- if (!is.null(fit$family) &&
                 fit$family %in% c("poisson", "quasipoisson")) {
    fit$family
  } else {
    "poisson"
  }
  structure <- fit$structure %||% "independent"
  knots <- if (!is.null(fit$ndx)) as.integer(fit$ndx) else c(15L, 15L)

  pick <- function(eta) {
    eta <- as.numeric(eta)
    switch(type,
      eta = eta,
      diff = eta[seq_len(n2)] - eta[(n2 + 1L):n_fine],
      shared = 0.5 * (eta[seq_len(n2)] + eta[(n2 + 1L):n_fine])
    )
  }

  fit_at_la <- function(la_use) {
    clgam_contrast(
      y = fit$y, x1 = fit$x1, x2 = fit$x2,
      exposure = fit$efine,
      C1 = CC$C1, C2 = CC$C2,
      knots = knots,
      elements = FALSE,
      trace = FALSE,
      parold = la_use,
      bold = bold,
      maxit = c(1L, 1L),
      family = fam_arg,
      structure = structure
    )
  }

  se_b <- switch(type,
    diff = as.numeric(fit$sd.dif),
    shared = as.numeric(fit$sd.shared),
    eta = as.numeric(fit$sd.eta)
  )
  d0 <- pick(fit_at_la(la)$eta)
  n_q <- length(d0)
  J <- matrix(0, nrow = n_q, ncol = length(la))
  for (j in seq_along(la)) {
    la_p <- la
    la_p[j] <- la[j] * exp(eps)
    J[, j] <- (pick(fit_at_la(la_p)$eta) - d0) / eps
  }
  Vrho_diag <- 2 / pmax(edf, 1e-8)
  extra <- as.numeric(J^2 %*% Vrho_diag)
  sqrt(pmax(se_b^2 + extra, 0))
}

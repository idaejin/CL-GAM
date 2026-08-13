#' Summarize a CL-GAMM fit
#'
#' Reviewer-facing summary of a \code{\link{clgam}} / \code{\link{pois_SOP}}
#' object. Reports family and type, PIRLS convergence, Pearson dispersion
#' \eqn{\phi} (and SE scaling under quasi-Poisson), SOP \eqn{\tau^2} and
#' effective df by smooth, overlap \eqn{\kappa} by fine-scale covariate
#' (manuscript eq. 21), and AIC/BIC when \code{elements=TRUE}.
#'
#' @param object a \code{clgam} fit
#' @param ... unused
#' @return An object of class \code{"summary.clgam"}.
#' @seealso \code{\link{clgam}}, \code{\link{kappa_diagnostic}}
#' @export
#' @method summary clgam
summary.clgam <- function(object, ...) {
  fam <- object$family %||% "poisson"
  if (fam %in% c("spatial", "contrast")) fam <- "poisson"
  type <- object$type %||% if (inherits(object, "clgam_contrast")) {
    "contrast"
  } else {
    "spatial"
  }
  kappa_by <- tryCatch(.clgam_kappa_by(object), error = function(e) numeric(0))
  terms <- .clgam_term_table(object, kappa_by = kappa_by)
  aic_joint <- if (!is.null(object$aic_total)) object$aic_total else {
    a <- object$aic
    if (!is.null(a) && length(a) == 1L) a else NULL
  }
  bic_joint <- if (!is.null(object$bic_total)) object$bic_total else {
    b <- object$bic
    if (!is.null(b) && length(b) == 1L) b else NULL
  }
  out <- list(
    call = object$call,
    type = type,
    family = fam,
    method = object$method %||% "SOP",
    re_term = object$re_term %||% "none",
    n_coarse = if (!is.null(object$y)) length(object$y) else NA_integer_,
    n_fine = if (!is.null(object$eta)) length(object$eta) else NA_integer_,
    niter = object$niter,
    maxit = object$maxit,
    thr = object$thr,
    tol = object$tol,
    elapsed = object$elapsed.time,
    diverged = isTRUE(object$diverged),
    converged = .clgam_infer_converged(object),
    var.comp = object$var.comp,
    edf = object$edf,
    terms = terms,
    ed = object$ed,
    aic = aic_joint,
    bic = bic_joint,
    aic_group = if (!is.null(object$aic) && length(object$aic) > 1L) object$aic else NULL,
    bic_group = if (!is.null(object$bic) && length(object$bic) > 1L) object$bic else NULL,
    dev = object$dev,
    phi = object$phi,
    deviance_df = object$deviance_df,
    df.residual = object$df.residual,
    kappa = object$orth.info$kappa,
    kappa_by = kappa_by,
    orth.smooth = object$orth.smooth,
    orth.applied = object$orth.info$applied,
    structure = object$structure,
    mean_se_dif = if (!is.null(object$sd.dif)) mean(object$sd.dif) else NULL,
    mean_se_shared = if (!is.null(object$sd.shared)) mean(object$sd.shared) else NULL,
    mean_se_dif2 = if (!is.null(object$sd.dif2)) mean(object$sd.dif2) else NULL
  )
  class(out) <- "summary.clgam"
  out
}

#' @param x a \code{"summary.clgam"} object
#' @param digits significant digits for \eqn{\tau^2} and \eqn{\kappa}
#' @rdname summary.clgam
#' @export
#' @method print summary.clgam
print.summary.clgam <- function(x, digits = 4, ...) {
  cat("CL-GAMM summary\n")
  if (!is.null(x$call)) {
    cat("Call: ")
    print(x$call)
  }
  cat("Type: ", x$type %||% "spatial",
      "    Family: ", x$family %||% "poisson",
      "    Method: ", x$method %||% "SOP",
      if (!is.null(x$re_term) && !identical(x$re_term, "none"))
        paste0("    re: ", x$re_term) else "",
      if (identical(x$type, "contrast") && !is.null(x$structure))
        paste0("    structure: ", x$structure) else "",
      "\n", sep = "")
  if (is.finite(x$n_coarse) || is.finite(x$n_fine)) {
    cat("n (coarse): ", x$n_coarse, "    n (fine): ", x$n_fine, "\n", sep = "")
  }

  cat("\nConvergence\n")
  conv <- x$converged
  conv_lab <- if (isTRUE(conv)) {
    "TRUE"
  } else if (identical(conv, FALSE)) {
    "FALSE"
  } else {
    "unknown (fit predates stored PIRLS diagnostics)"
  }
  cat("  Converged: ", conv_lab, "\n", sep = "")
  if (isTRUE(x$diverged)) {
    cat("  Variance-component update diverged; last finite tau^2 retained.\n")
  }
  maxit1 <- if (!is.null(x$maxit)) x$maxit[1] else NA
  cat("  Outer PIRLS: ", x$niter,
      if (is.finite(maxit1)) paste0(" / ", maxit1) else "", "\n", sep = "")
  if (!is.null(x$tol) && is.finite(x$tol)) {
    thr1 <- if (!is.null(x$thr)) x$thr[1] else NA
    cat("  Rel. change in eta: ", format(x$tol, scientific = TRUE, digits = 3),
        if (is.finite(thr1)) paste0("  (threshold ", format(thr1, scientific = TRUE, digits = 3), ")") else "",
        "\n", sep = "")
  }
  if (!is.null(x$elapsed) && is.finite(x$elapsed)) {
    cat("  Elapsed: ", round(x$elapsed, 3), " s\n", sep = "")
  }

  cat("\nDispersion\n")
  if (!is.null(x$phi) && is.finite(x$phi)) {
    cat("  Pearson phi: ", signif(x$phi, digits), sep = "")
    if (identical(x$family, "quasipoisson")) {
      cat("\n  Quasi-Poisson: SEs scaled by sqrt(phi) = ",
          signif(sqrt(x$phi), digits), "\n", sep = "")
    } else {
      cat("  (Poisson diagnostic; near 1 if mean-variance holds)\n")
    }
  } else {
    cat("  Pearson phi: not stored\n")
  }

  has_cov <- (!is.null(x$kappa) && is.finite(x$kappa)) ||
    (length(x$kappa_by) > 0L) ||
    !is.null(x$orth.smooth)
  if (has_cov) {
    cat("\nIdentifiability\n")
    if (!is.null(x$orth.smooth)) {
      cat("  orth.smooth: ", x$orth.smooth,
          if (isTRUE(x$orth.applied)) " (spatial bases residualized vs A_f)" else "",
          "\n", sep = "")
    }
    if (!is.null(x$kappa) && is.finite(x$kappa)) {
      cat("  Overall kappa (A_f): ", signif(x$kappa, digits), "\n", sep = "")
    }
  }

  cat("\nSmooth terms (SOP tau^2 and ED; kappa vs that covariate's A_k)\n")
  if (!is.null(x$terms) && nrow(x$terms)) {
    .clgam_print_term_table(x$terms, digits = digits)
  } else {
    cat("  (none stored)\n")
  }

  cat("\nModel fit\n")
  dev <- x$dev
  if (!is.null(dev)) {
    if (length(dev) > 1L) {
      cat("  Deviance (by group): ", paste(round(dev, 3), collapse = ", "),
          "  total: ", round(sum(dev), 3), "\n", sep = "")
    } else {
      cat("  Deviance: ", round(dev, 3), "\n", sep = "")
    }
  }
  if (!is.null(x$ed) && is.finite(x$ed)) {
    cat("  Model ED (trace, AIC): ", round(x$ed, 3), "\n", sep = "")
  }
  if (!is.null(x$aic) && is.finite(x$aic)) {
    cat("  AIC: ", round(x$aic, 3), sep = "")
    if (!is.null(x$bic) && is.finite(x$bic)) cat("    BIC: ", round(x$bic, 3), sep = "")
    cat("\n")
    if (identical(x$family, "quasipoisson")) {
      cat("  AIC/BIC use the Poisson log-likelihood (same fitted means).\n")
    }
    if (!is.null(x$aic_group)) {
      cat("  AIC by group: ", paste(round(x$aic_group, 3), collapse = ", "),
          "  (each includes the shared ED; do not sum)\n", sep = "")
    }
  } else {
    cat("  AIC/BIC: not computed (refit with elements=TRUE)\n")
  }
  if (identical(x$type, "contrast") && !is.null(x$mean_se_dif)) {
    cat("\nContrast inference (Bayesian SE, conditional on tau^2)\n")
    cat("  Mean se(eta1-eta2) [joint]: ", signif(x$mean_se_dif, digits), sep = "")
    if (!is.null(x$mean_se_dif2)) {
      cat("    indep. approx: ", signif(x$mean_se_dif2, digits), sep = "")
    }
    cat("\n")
    if (!is.null(x$mean_se_shared)) {
      cat("  Mean se((eta1+eta2)/2): ", signif(x$mean_se_shared, digits), "\n", sep = "")
    }
    cat("  Pointwise intervals: clgam_contrast_infer(fit); plot(fit, which=7)\n")
  }
  invisible(x)
}

#' Infer convergence for current and older fit objects.
#' @keywords internal
.clgam_infer_converged <- function(object) {
  if (!is.null(object$converged)) return(isTRUE(object$converged))
  if (isTRUE(object$diverged)) return(FALSE)
  NA
}

#' SOP term names: spatial penalties then smooth covariates.
#' @keywords internal
.clgam_edf_names <- function(fit) {
  n_ed <- length(fit$edf)
  if (identical(fit$type, "contrast") || inherits(fit, "clgam_contrast")) {
    nms <- if (identical(fit$structure, "contrast")) {
      c("spatial.shared.x1", "spatial.shared.x2",
        "spatial.contrast.x1", "spatial.contrast.x2")
    } else {
      c("spatial.g1.x1", "spatial.g1.x2", "spatial.g2.x1", "spatial.g2.x2")
    }
    return(nms[seq_len(min(n_ed, length(nms)))])
  }
  sm <- .clgam_mat_colnames(fit$nlcovfine, "s")
  nms <- c("spatial.x1", "spatial.x2", sm)
  if (identical(fit$re_term, "coarse") ||
      (!is.null(names(fit$var.comp)) && "re.coarse" %in% names(fit$var.comp))) {
    nms <- c(nms, "re.coarse")
  }
  if (length(nms) != n_ed) {
    extra <- max(n_ed - 2L, 0L)
    nms <- c("spatial.x1", "spatial.x2", paste0("s", seq_len(extra)))
    nms <- nms[seq_len(n_ed)]
  }
  nms
}

#' Term table: tau^2 / SOP ED / kappa by covariate (linear rows have no tau^2).
#' @keywords internal
.clgam_term_table <- function(object, kappa_by = NULL) {
  vc <- as.numeric(object$var.comp)
  edf <- as.numeric(object$edf)
  n <- max(length(vc), length(edf))
  if (!n && is.null(object$lcovfine)) {
    return(data.frame(
      term = character(0), tau2 = numeric(0), ed = numeric(0),
      kappa = numeric(0), stringsAsFactors = FALSE
    ))
  }
  nms <- names(object$var.comp)
  if (is.null(nms) || !length(nms) || !any(nzchar(nms))) {
    nms <- names(object$edf)
  }
  if (is.null(nms) || !length(nms) || !any(nzchar(nms))) {
    nms <- .clgam_edf_names(object)
  }
  if (length(nms) < n) nms <- c(nms, paste0("term", seq_len(n - length(nms))))
  nms <- nms[seq_len(n)]
  if (length(vc) < n) vc <- c(vc, rep(NA_real_, n - length(vc)))
  if (length(edf) < n) edf <- c(edf, rep(NA_real_, n - length(edf)))

  if (is.null(kappa_by)) {
    kappa_by <- tryCatch(.clgam_kappa_by(object), error = function(e) numeric(0))
  }
  kap <- rep(NA_real_, n)
  if (length(kappa_by)) {
    hit <- match(nms, names(kappa_by))
    ok <- !is.na(hit)
    kap[ok] <- unname(kappa_by[hit[ok]])
  }

  tab <- data.frame(
    term = nms, tau2 = vc[seq_len(n)], ed = edf[seq_len(n)], kappa = kap,
    stringsAsFactors = FALSE
  )

  lin_names <- .clgam_mat_colnames(object$lcovfine, "x")
  if (length(lin_names)) {
    for (nm in lin_names) {
      if (nm %in% tab$term) next
      tab <- rbind(tab, data.frame(
        term = nm, tau2 = NA_real_, ed = NA_real_,
        kappa = unname(kappa_by[nm]) %||% NA_real_,
        stringsAsFactors = FALSE
      ))
    }
  }
  rownames(tab) <- NULL
  tab
}

#' @keywords internal
.clgam_print_term_table <- function(d, digits = 4) {
  fmt_g <- function(x) {
    out <- rep("--", length(x))
    ok <- !is.na(x)
    if (any(ok)) out[ok] <- formatC(x[ok], digits = digits, format = "g")
    out
  }
  fmt_ed <- function(x) {
    out <- rep("--", length(x))
    ok <- !is.na(x)
    if (any(ok)) out[ok] <- sprintf("%.3f", x[ok])
    out
  }
  term <- as.character(d$term)
  tau <- fmt_g(d$tau2)
  ed <- fmt_ed(d$ed)
  kap <- fmt_g(d$kappa)
  w_term <- max(nchar(c("term", term)))
  w_tau <- max(nchar(c("tau2", tau)))
  w_ed <- max(nchar(c("ED", ed)))
  w_k <- max(nchar(c("kappa", kap)))
  row <- function(a, b, c, k) {
    sprintf("  %-*s  %*s  %*s  %*s", w_term, a, w_tau, b, w_ed, c, w_k, k)
  }
  cat(row("term", "tau2", "ED", "kappa"), "\n")
  for (i in seq_along(term)) {
    cat(row(term[i], tau[i], ed[i], kap[i]), "\n")
  }
}

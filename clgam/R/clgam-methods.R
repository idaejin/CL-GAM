#' Extractor and plotting methods for \code{clgam} fits
#'
#' @name clgam-methods
NULL

# ---- extractors --------------------------------------------------------------

#' @rdname clgam-methods
#' @param object a \code{clgam} fit
#' @param which for \code{coef}: \code{"all"} (default), \code{"fixed"}, or
#'   \code{"random"}
#' @param ... unused
#' @export
#' @method coef clgam
coef.clgam <- function(object, which = c("all", "fixed", "random"), ...) {
  which <- match.arg(which)
  bf <- object$b.fixed
  br <- object$b.random
  names(bf) <- paste0("beta", seq_along(bf))
  names(br) <- paste0("u", seq_along(br))
  switch(which,
    fixed = bf,
    random = br,
    all = c(bf, br)
  )
}

#' @rdname clgam-methods
#' @param type for \code{fitted}/\code{predict}: \code{"mu"} (coarse mean),
#'   \code{"eta"} (fine log-risk), \code{"gamma"} (fine intensity),
#'   \code{"response"} (alias of \code{mu}), or \code{"exp.eta"}
#' @export
#' @method fitted clgam
fitted.clgam <- function(object,
                         type = c("mu", "eta", "gamma", "response", "exp.eta"),
                         ...) {
  type <- match.arg(type)
  switch(type,
    mu = , response = object$mu,
    eta = object$eta,
    gamma = object$gamma,
    exp.eta = exp(object$eta)
  )
}

#' @rdname clgam-methods
#' @param type residual type: deviance, Pearson, or response (\code{y - mu})
#' @export
#' @method residuals clgam
residuals.clgam <- function(object,
                            type = c("deviance", "pearson", "response"),
                            ...) {
  type <- match.arg(type)
  y <- object$y
  mu <- object$mu
  if (is.null(y) || is.null(mu)) {
    stop("Fit does not store y/mu; refit with a current clgam version.", call. = FALSE)
  }
  switch(type,
    response = y - mu,
    pearson = (y - mu) / sqrt(pmax(mu, 1e-12)),
    deviance = {
      # signed deviance residual for Poisson
      r <- sign(y - mu) * sqrt(pmax(
        2 * (y * log(ifelse(y == 0, 1, y / pmax(mu, 1e-12))) - (y - mu)),
        0
      ))
      r
    }
  )
}

#' @rdname clgam-methods
#' @export
#' @method nobs clgam
nobs.clgam <- function(object, ...) {
  if (!is.null(object$y)) length(object$y) else length(object$mu)
}

#' @rdname clgam-methods
#' @export
#' @method deviance clgam
deviance.clgam <- function(object, ...) {
  d <- object$dev
  if (length(d) > 1L) sum(d) else d
}

#' @rdname clgam-methods
#' @export
#' @method logLik clgam
logLik.clgam <- function(object, ...) {
  # Poisson log-lik up to constant; ed used as df when available
  y <- object$y
  mu <- object$mu
  if (is.null(y) || is.null(mu)) {
    stop("Fit does not store y/mu; cannot compute logLik.", call. = FALSE)
  }
  ll <- sum(stats::dpois(y, lambda = pmax(mu, 1e-12), log = TRUE))
  df <- if (!is.null(object$ed)) object$ed else length(object$b.fixed) + length(object$b.random)
  structure(ll, df = df, nobs = length(y), class = "logLik")
}

#' @rdname clgam-methods
#' @param k penalty in AIC (default 2)
#' @export
#' @method AIC clgam
AIC.clgam <- function(object, ..., k = 2) {
  if (!is.null(object$aic) && isTRUE(all.equal(k, 2))) {
    a <- object$aic
    if (length(a) > 1L) {
      # A vector `aic` (e.g. from clgam_contrast()/pois_incat_SOP()) holds
      # one AIC per category, each already including the FULL shared `ed`
      # of the joint two-population fit. Summing them naively double-counts
      # `ed`; use the precomputed joint total when available (see
      # pois_incat_SOP.R), falling back to the old (double-counting)
      # behaviour only for older fit objects that lack it.
      return(if (!is.null(object$aic_total)) object$aic_total else sum(a))
    }
    return(a)
  }
  ll <- logLik(object)
  -2 * as.numeric(ll) + k * attr(ll, "df")
}

#' @rdname clgam-methods
#' @export
#' @method BIC clgam
BIC.clgam <- function(object, ...) {
  if (!is.null(object$bic)) {
    b <- object$bic
    if (length(b) > 1L) {
      # See AIC.clgam(): per-category `bic` entries each already include the
      # full shared `ed`; use the precomputed joint total to avoid
      # double-counting it (and mismatched per-category log(n) terms).
      return(if (!is.null(object$bic_total)) object$bic_total else sum(b))
    }
    return(b)
  }
  ll <- logLik(object)
  -2 * as.numeric(ll) + log(nobs(object)) * attr(ll, "df")
}

#' Predict from a CL-GAMM fit
#'
#' By default returns fitted values on the training fine/coarse grid.
#' Supplying \code{newdata} is not yet supported (would require rebuilding
#' P-spline bases at new locations).
#'
#' @param object a \code{clgam} fit
#' @param newdata currently must be \code{NULL}
#' @param type \code{"eta"}, \code{"mu"}, \code{"gamma"}, \code{"response"},
#'   or \code{"exp.eta"}
#' @param se.fit if \code{TRUE} and SEs were computed (\code{elements=TRUE}),
#'   return a list with \code{fit} and \code{se.fit}
#' @param ... unused
#' @export
#' @method predict clgam
predict.clgam <- function(object,
                          newdata = NULL,
                          type = c("eta", "mu", "gamma", "response", "exp.eta"),
                          se.fit = FALSE,
                          ...) {
  if (!is.null(newdata)) {
    stop("predict.clgam: newdata not implemented yet; use fitted grid.", call. = FALSE)
  }
  type <- match.arg(type)
  fit <- fitted(object, type = type)
  if (!isTRUE(se.fit)) return(fit)

  se <- switch(type,
    eta = object$sd.eta,
    exp.eta = object$sd.exp.eta,
    mu = , response = , gamma = NULL
  )
  if (is.null(se)) {
    warning("SEs not available for type='", type, "' (need elements=TRUE and type eta/exp.eta).",
            call. = FALSE)
  }
  list(fit = fit, se.fit = se)
}

# ---- plotting ----------------------------------------------------------------

#' Plot a CL-GAMM fit
#'
#' @param x a \code{clgam} fit
#' @param which which plot(s): \code{1} fine \eqn{\hat\eta}; \code{2} y vs
#'   \eqn{\hat\mu}; \code{3} residuals; \code{4} smooth; \code{5} contrast;
#'   \code{6} coarse map of \eqn{y-\hat\mu}
#' @param ask ask before each plot when length(which) > 1
#' @param sf_fine,sf_coarse optional \pkg{sf} layers (else \code{x$sf_fine} /
#'   \code{x$sf_coarse})
#' @param g_true optional true smooth vector aligned with fine cells
#' @param ... graphical parameters (\code{cex}, \code{pch}, \code{lwd}, ...)
#' @export
#' @method plot clgam
plot.clgam <- function(x,
                       which = 1:2,
                       ask = prod(graphics::par("mfcol")) < length(which) &&
                         grDevices::dev.interactive(),
                       sf_fine = NULL,
                       sf_coarse = NULL,
                       g_true = NULL,
                       ...) {
  show <- sort(unique(as.integer(which)))
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  dots <- list(...)
  if (is.null(sf_fine)) sf_fine <- x$sf_fine
  if (is.null(sf_coarse)) sf_coarse <- x$sf_coarse
  if (is.null(g_true)) g_true <- x$g_true

  if (1L %in% show) .clgam_plot_eta(x, dots, sf_fine = sf_fine)
  if (2L %in% show) .clgam_plot_obs_fit(x, dots, sf_coarse = sf_coarse)
  if (3L %in% show) .clgam_plot_resid(x, dots)
  if (4L %in% show) .clgam_plot_smooths(x, dots, g_true = g_true)
  if (5L %in% show) .clgam_plot_contrast(x, dots, sf_fine = sf_fine)
  if (6L %in% show) .clgam_plot_coarse_map(x, dots, sf_coarse = sf_coarse)
  invisible(x)
}

#' Map numeric values to colours
#' @keywords internal
.clgam_col_ramp <- function(v, palette = "Plasma", n = 100L) {
  v <- as.numeric(v)
  ok <- is.finite(v)
  cols <- rep(NA_character_, length(v))
  if (!any(ok)) return(cols)
  r <- range(v[ok])
  if (diff(r) < 1e-12) r <- r + c(-1, 1) * 1e-6
  pal <- grDevices::hcl.colors(n, palette)
  idx <- as.integer(cut(v[ok], breaks = seq(r[1], r[2], length.out = n + 1L),
                        include.lowest = TRUE))
  cols[ok] <- pal[idx]
  attr(cols, "range") <- r
  attr(cols, "palette") <- pal
  cols
}

#' Draw a vertical colour bar in the right margin
#' @keywords internal
.clgam_colorbar <- function(zlim, palette, title = NULL) {
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  din <- graphics::par("din")
  # thin strip just outside plot
  x0 <- usr[2] + 0.02 * diff(usr[1:2])
  x1 <- usr[2] + 0.06 * diff(usr[1:2])
  ys <- seq(usr[3], usr[4], length.out = length(palette) + 1L)
  for (i in seq_along(palette)) {
    graphics::rect(x0, ys[i], x1, ys[i + 1L], col = palette[i], border = NA)
  }
  graphics::rect(x0, usr[3], x1, usr[4], border = "grey30")
  labs <- pretty(zlim, n = 4)
  labs <- labs[labs >= zlim[1] & labs <= zlim[2]]
  ylab <- usr[3] + (labs - zlim[1]) / diff(zlim) * diff(usr[3:4])
  graphics::text(x1 + 0.01 * diff(usr[1:2]), ylab, labels = signif(labs, 3),
                 adj = 0, cex = 0.7, xpd = TRUE)
  if (!is.null(title)) {
    graphics::mtext(title, side = 4, line = 0.2, cex = 0.75)
  }
  invisible(NULL)
}

#' @keywords internal
.clgam_plot_eta <- function(x, dots, sf_fine = NULL) {
  if (is.null(x$eta)) {
    warning("Cannot plot eta: missing eta on fit object.", call. = FALSE)
    return(invisible(NULL))
  }
  eta <- x$eta
  cols <- .clgam_col_ramp(eta, palette = dots$palette %||% "Plasma")
  main <- dots$main %||% expression(hat(eta) ~ "(fine)")

  if (!is.null(sf_fine) && requireNamespace("sf", quietly = TRUE) &&
      inherits(sf_fine, "sf") && nrow(sf_fine) == length(eta)) {
    op <- graphics::par(mar = c(4, 4, 2.5, 3.5))
    on.exit(graphics::par(op), add = TRUE)
    plot(sf::st_geometry(sf_fine), col = cols, border = "grey30", lwd = 0.3,
         main = main, axes = TRUE, reset = FALSE)
    if (!is.null(x$sf_coarse) && inherits(x$sf_coarse, "sf")) {
      plot(sf::st_geometry(x$sf_coarse), add = TRUE, border = "black", lwd = 1.1)
    }
    .clgam_colorbar(attr(cols, "range"), attr(cols, "palette"))
  } else if (!is.null(x$x1) && !is.null(x$x2)) {
    op <- graphics::par(mar = c(4, 4, 2.5, 3.5))
    on.exit(graphics::par(op), add = TRUE)
    cex <- dots$cex %||% 0.85
    pch <- dots$pch %||% 16
    graphics::plot(x$x1, x$x2, col = cols, pch = pch, cex = cex,
                   xlab = "x1", ylab = "x2", main = main, asp = 1)
    .clgam_colorbar(attr(cols, "range"), attr(cols, "palette"))
  } else {
    warning("Need sf_fine or x1/x2 to plot eta.", call. = FALSE)
  }
}

#' @keywords internal
.clgam_plot_obs_fit <- function(x, dots, sf_coarse = NULL) {
  if (is.null(x$y) || is.null(x$mu)) {
    warning("Cannot plot obs vs fit: missing y/mu.", call. = FALSE)
    return(invisible(NULL))
  }
  pch <- dots$pch %||% 16
  r <- stats::cor(x$mu, x$y)
  graphics::plot(
    x$mu, x$y, pch = pch,
    xlab = expression(hat(mu)), ylab = "y",
    main = sprintf("Coarse y vs fit  (r=%.2f)", r)
  )
  graphics::abline(0, 1, lty = 2, col = "grey40")
}

#' @keywords internal
.clgam_plot_coarse_map <- function(x, dots, sf_coarse = NULL, field = c("resid", "y", "mu")) {
  field <- match.arg(field)
  if (is.null(sf_coarse) || !requireNamespace("sf", quietly = TRUE) ||
      !inherits(sf_coarse, "sf")) {
    message("Coarse map needs sf_coarse.")
    return(invisible(NULL))
  }
  v <- switch(field,
    resid = x$y - x$mu,
    y = x$y,
    mu = x$mu
  )
  if (length(v) != nrow(sf_coarse)) {
    warning("sf_coarse rows != length(y).", call. = FALSE)
    return(invisible(NULL))
  }
  pal <- if (field == "resid") "RdBu" else "YlOrRd"
  cols <- .clgam_col_ramp(v, palette = pal)
  op <- graphics::par(mar = c(4, 4, 2.5, 3.5))
  on.exit(graphics::par(op), add = TRUE)
  main <- switch(field,
    resid = expression(y - hat(mu)),
    y = "y (coarse)",
    mu = expression(hat(mu))
  )
  plot(sf::st_geometry(sf_coarse), col = cols, border = "grey20", lwd = 0.7,
       main = main, axes = TRUE, reset = FALSE)
  .clgam_colorbar(attr(cols, "range"), attr(cols, "palette"))
}

#' @keywords internal
.clgam_plot_resid <- function(x, dots) {
  r <- tryCatch(residuals(x, type = "deviance"), error = function(e) NULL)
  if (is.null(r) || is.null(x$mu)) {
    warning("Cannot plot residuals.", call. = FALSE)
    return(invisible(NULL))
  }
  pch <- dots$pch %||% 16
  graphics::plot(
    x$mu, r, pch = pch,
    xlab = expression(hat(mu)), ylab = "Dev. residual",
    main = "Deviance residuals"
  )
  graphics::abline(h = 0, lty = 2, col = "grey40")
}

#' @keywords internal
.clgam_plot_smooths <- function(x, dots, g_true = NULL) {
  if (is.null(x$nleffects)) {
    message("No nonlinear effects stored on this fit.")
    return(invisible(NULL))
  }
  nk <- ncol(x$nleffects)
  for (k in seq_len(nk)) {
    xx <- if (!is.null(x$nlcovfine)) {
      x$nlcovfine[, min(k, ncol(x$nlcovfine))]
    } else {
      seq_len(nrow(x$nleffects))
    }
    ord <- order(xx)
    xx <- xx[ord]
    yy <- x$nleffects[ord, k]
    # center for display comparability with g_true
    yy <- yy - mean(yy)

    ylim <- range(yy, finite = TRUE)
    if (!is.null(x$sdnleffects)) {
      se <- x$sdnleffects[ord, k]
      ylim <- range(c(ylim, yy - 1.96 * se, yy + 1.96 * se), finite = TRUE)
    }
    g_plot <- NULL
    if (!is.null(g_true) && length(g_true) == length(ord)) {
      g_plot <- g_true[ord] - mean(g_true)
      ylim <- range(c(ylim, g_plot), finite = TRUE)
    }

    xlab <- {
      cn <- colnames(x$nlcovfine)
      if (!is.null(cn) && length(cn) >= k && nzchar(cn[k])) cn[k] else paste0("z[", k, "]")
    }
    graphics::plot(
      xx, yy, type = "n",
      xlab = if (!is.null(x$nlcovfine)) xlab else "index",
      ylab = "effect (centred)",
      main = paste("Nonlinear effect", k),
      ylim = ylim + c(-1, 1) * 0.05 * max(diff(ylim), 1e-6)
    )
    if (!is.null(x$sdnleffects)) {
      se <- x$sdnleffects[ord, k]
      graphics::polygon(
        c(xx, rev(xx)),
        c(yy - 1.96 * se, rev(yy + 1.96 * se)),
        col = grDevices::adjustcolor("steelblue", 0.25),
        border = NA
      )
    }
    graphics::lines(xx, yy, lwd = dots$lwd %||% 2, col = "steelblue")
    if (!is.null(g_plot)) {
      graphics::lines(xx, g_plot, lwd = 2, lty = 2, col = "grey20")
      graphics::legend("topright",
                       legend = c(expression(hat(g)(z)), expression(g[true](z))),
                       col = c("steelblue", "grey20"), lwd = 2, lty = c(1, 2),
                       bty = "n", cex = 0.85)
    }
    graphics::abline(h = 0, lty = 3, col = "grey70")
  }
}

#' @keywords internal
.clgam_plot_contrast <- function(x, dots, sf_fine = NULL) {
  if (!inherits(x, "clgam_contrast") &&
      !identical(x$type, "contrast") &&
      !identical(x$family, "contrast")) {
    message("Contrast plot only for two-group fits (clgam_contrast).")
    return(invisible(NULL))
  }
  if (is.null(x$eta)) return(invisible(NULL))
  n2 <- length(x$eta) %/% 2L
  dif <- x$eta[seq_len(n2)] - x$eta[(n2 + 1L):(2L * n2)]
  cols <- .clgam_col_ramp(dif, palette = "RdYlBu")

  if (!is.null(sf_fine) && requireNamespace("sf", quietly = TRUE) &&
      inherits(sf_fine, "sf") && nrow(sf_fine) >= n2) {
    geom <- sf::st_geometry(sf_fine)[seq_len(n2)]
    op <- graphics::par(mar = c(4, 4, 2.5, 3.5))
    on.exit(graphics::par(op), add = TRUE)
    plot(geom, col = cols, border = "grey30", lwd = 0.3,
         main = expression(hat(eta)[1] - hat(eta)[2]), axes = TRUE, reset = FALSE)
    .clgam_colorbar(attr(cols, "range"), attr(cols, "palette"))
  } else if (!is.null(x$x1)) {
    graphics::plot(
      x$x1[seq_len(n2)], x$x2[seq_len(n2)],
      col = cols, pch = dots$pch %||% 16, cex = dots$cex %||% 0.7,
      xlab = "x1", ylab = "x2", asp = 1,
      main = expression(hat(eta)[1] - hat(eta)[2])
    )
  }
}

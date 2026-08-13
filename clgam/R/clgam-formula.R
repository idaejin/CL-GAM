#' Univariate or spatial P-spline term (formula interface)
#'
#' Used on the right-hand side of a \code{\link{clgam}} formula.
#' \code{s(x1, x2)} is the anisotropic spatial field; \code{s(z)} is a
#' univariate P-spline. A bare name \code{z} (no \code{s()}) is a
#' \emph{linear} covariate, equivalent to the \code{linear=} argument.
#'
#' P-spline settings follow the package notation:
#' \code{s(x1, x2, ndx = c(15, 15), bdeg = 3, pord = 2)}.
#'
#' \code{\link{clgam}} walks the formula language object and looks for the
#' symbol \code{s}; it does not call \code{mgcv::s}. If both packages are
#' attached, \code{s(...)} in a \code{clgam} formula still works. Use
#' \code{?clgam::s} (or attach \pkg{clgam} last) for this help page.
#'
#' @param ... one name (\code{s(z)}) or two (\code{s(x1, x2)})
#' @param ndx number of internal B-spline intervals (passed as \code{ndx} /
#'   \code{ndxnl} to \code{\link{pois_SOP}}). Length 1 for a univariate
#'   smooth; length 2 (or recycled) for space. Alias: \code{k}.
#' @param bdeg B-spline degree (default \code{3})
#' @param pord difference penalty order (default \code{2})
#' @param k alias for \code{ndx} (mgcv-style). Ignored if \code{ndx} is set.
#' @param level \code{"fine"} (default) or \code{"coarse"} (Case B/C
#'   aggregated-scale smooth)
#' @return An object of class \code{"clgam.smooth.spec"} (not used directly)
#' @seealso \code{\link{clgam}}, \code{\link{simulate_ata_scenarios}}
#' @export
s <- function(..., ndx = NULL, bdeg = 3, pord = 2, k = NULL,
              level = "fine") {
  level <- match.arg(level, c("fine", "coarse"))
  vars <- as.list(substitute(list(...)))[-1L]
  varn <- vapply(vars, function(v) {
    if (is.character(v)) v else paste(deparse(v), collapse = "")
  }, character(1))
  ndx <- .clgam_resolve_ndx(ndx, k)
  structure(
    list(term = varn, ndx = ndx, k = ndx, bdeg = bdeg, pord = pord,
         level = level),
    class = "clgam.smooth.spec"
  )
}

#' Resolve \code{ndx} vs alias \code{k}
#' @noRd
.clgam_resolve_ndx <- function(ndx, k) {
  if (!is.null(ndx) && !is.null(k)) {
    if (!isTRUE(all.equal(as.numeric(ndx), as.numeric(k)))) {
      warning("Both ndx= and k= in s(); using ndx=.", call. = FALSE)
    }
    return(ndx)
  }
  if (!is.null(ndx)) return(ndx)
  k
}

#' Recycle a length-1 P-spline setting to length 2 (spatial margins)
#' @noRd
.clgam_recycle2 <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.numeric(x)
  if (length(x) == 1L) c(x, x) else x
}

#' Parse a CL-GAM formula \code{y ~ s(x1, x2) + s(z)} / \code{+ z}
#' @keywords internal
.clgam_parse_formula <- function(formula, data = NULL, env = parent.frame()) {
  if (!inherits(formula, "formula") || length(formula) < 3L) {
    stop("Expected a two-sided formula, e.g. y ~ s(x1, x2) + s(z).",
         call. = FALSE)
  }
  env_data <- .clgam_eval_env(data, env)
  terms <- .clgam_collect_terms(formula[[3L]], env = env_data)
  specs <- terms$smooth
  if (!length(specs)) {
    stop("Formula must contain at least s(x1, x2).", call. = FALSE)
  }

  spatial <- NULL
  uni <- list()
  for (sp in specs) {
    if (length(sp$term) == 2L) {
      if (!is.null(spatial)) {
        stop("Only one spatial term s(x1, x2) is allowed.", call. = FALSE)
      }
      spatial <- sp
    } else if (length(sp$term) == 1L) {
      uni[[length(uni) + 1L]] <- sp
    } else {
      stop("s() takes one covariate or two spatial coordinates.", call. = FALSE)
    }
  }
  if (is.null(spatial)) {
    stop("Formula needs a spatial term s(x1, x2).", call. = FALSE)
  }

  y <- eval(formula[[2L]], envir = env_data)
  x1 <- .clgam_eval_term(spatial$term_expr[[1L]], env_data)
  x2 <- .clgam_eval_term(spatial$term_expr[[2L]], env_data)
  knots <- .clgam_recycle2(spatial$ndx %||% spatial$k)
  bdeg <- .clgam_recycle2(spatial$bdeg %||% 3)
  pord <- .clgam_recycle2(spatial$pord %||% 2)

  smooth_vars <- NULL
  nl_level <- NULL
  knots_nl <- NULL
  bdegnl <- NULL
  pordnl <- NULL
  if (length(uni)) {
    smooth_vars <- vector("list", length(uni))
    nl_level <- character(length(uni))
    knots_nl <- rep(NA_real_, length(uni))
    bdegnl <- rep(3, length(uni))
    pordnl <- rep(2, length(uni))
    cn <- character(length(uni))
    for (i in seq_along(uni)) {
      ui <- uni[[i]]
      smooth_vars[[i]] <- .clgam_eval_term(ui$term_expr[[1L]], env_data)
      nl_level[i] <- ui$level %||% "fine"
      cn[i] <- ui$term
      ndx_i <- ui$ndx %||% ui$k
      if (!is.null(ndx_i)) knots_nl[i] <- as.numeric(ndx_i)[1L]
      if (!is.null(ui$bdeg)) bdegnl[i] <- as.numeric(ui$bdeg)[1L]
      if (!is.null(ui$pord)) pordnl[i] <- as.numeric(ui$pord)[1L]
    }
    names(smooth_vars) <- cn
    if (all(is.na(knots_nl))) knots_nl <- NULL
  }

  linear_vars <- NULL
  if (length(terms$linear)) {
    linear_vars <- vector("list", length(terms$linear))
    cn <- character(length(terms$linear))
    for (i in seq_along(terms$linear)) {
      ui <- terms$linear[[i]]
      linear_vars[[i]] <- .clgam_eval_term(ui$term_expr, env_data)
      cn[i] <- ui$term
    }
    names(linear_vars) <- cn
  }

  list(
    y = y, x1 = x1, x2 = x2,
    smooth_vars = smooth_vars,
    linear_vars = linear_vars,
    smooth_level = nl_level,
    knots = knots,
    knots_nl = knots_nl,
    bdeg = bdeg,
    pord = pord,
    bdegnl = bdegnl,
    pordnl = pordnl,
    formula = formula
  )
}

#' Evaluate a formula term (symbol, call, or name string)
#' @keywords internal
.clgam_eval_term <- function(expr, env) {
  if (is.character(expr) && length(expr) == 1L) {
    expr <- as.name(expr)
  }
  eval(expr, envir = env)
}

#' Bind formula covariates; expand length-\code{n_coarse} vectors via partition \code{C}
#' @keywords internal
.clgam_bind_smooth <- function(cols, C) {
  if (is.null(cols) || !length(cols)) return(NULL)
  n_fine <- ncol(C)
  n_coarse <- nrow(C)
  out <- vector("list", length(cols))
  nms <- names(cols)
  if (is.null(nms) || !length(nms)) nms <- paste0("s", seq_along(cols))
  for (i in seq_along(cols)) {
    v <- as.numeric(cols[[i]])
    nm <- nms[i]
    if (length(v) == n_fine) {
      out[[i]] <- v
    } else if (length(v) == n_coarse) {
      if (!isTRUE(.is_partition_C(C))) {
        stop(
          "Covariate '", nm, "' has length n_coarse but C is not a 0-1 ",
          "partition. Expand it to the fine support first.",
          call. = FALSE
        )
      }
      out[[i]] <- v[.partition_groups(C)]
    } else {
      stop(
        "Covariate '", nm, "' has length ", length(v),
        "; expected ", n_fine, " (fine) or ", n_coarse, " (coarse).",
        call. = FALSE
      )
    }
  }
  M <- do.call(cbind, out)
  colnames(M) <- nms
  M
}

#' Fill C / exposure from a \code{simulate_ata()} list when omitted
#' @keywords internal
.clgam_pull_design <- function(data, C, exposure) {
  if (is.null(data) || !is.list(data)) {
    return(list(C = C, exposure = exposure))
  }
  if (is.null(C) && !is.null(data$C)) C <- data$C
  if (is.null(exposure)) {
    if (!is.null(data$efine)) {
      exposure <- data$efine
    } else if (!is.null(data$exposure)) {
      exposure <- data$exposure
    } else if (!is.null(data$ef)) {
      exposure <- data$ef
    }
  }
  list(C = C, exposure = exposure)
}

#' Evaluation environment: list/data.frame first, then calling env
#' @keywords internal
.clgam_eval_env <- function(data, env) {
  if (is.null(data)) return(env)
  if (is.environment(data)) return(data)
  if (is.list(data)) {
    return(list2env(as.list(data), parent = env, hash = TRUE))
  }
  env
}

.clgam_empty_terms <- function() list(smooth = list(), linear = list())

.clgam_merge_terms <- function(a, b) {
  list(smooth = c(a$smooth, b$smooth), linear = c(a$linear, b$linear))
}

#' Walk \code{+} / \code{s()} / linear names on a formula RHS
#' @keywords internal
.clgam_collect_terms <- function(expr, env) {
  if (is.null(expr)) return(.clgam_empty_terms())
  if (identical(expr, 1) || identical(expr, 1L)) return(.clgam_empty_terms())
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    if (!nzchar(nm) || identical(nm, ".")) {
      stop("Unsupported formula term '", nm, "'.", call. = FALSE)
    }
    return(list(
      smooth = list(),
      linear = list(list(term = nm, term_expr = expr))
    ))
  }
  if (!is.call(expr)) {
    stop("Cannot parse formula term: ", paste(deparse(expr), collapse = " "),
         call. = FALSE)
  }
  head <- expr[[1L]]
  if (identical(head, as.name("("))) {
    return(.clgam_collect_terms(expr[[2L]], env = env))
  }
  if (identical(head, as.name("+"))) {
    return(.clgam_merge_terms(
      .clgam_collect_terms(expr[[2L]], env = env),
      .clgam_collect_terms(expr[[3L]], env = env)
    ))
  }
  if (identical(head, as.name("s")) ||
      (is.call(head) && identical(head[[1L]], as.name("::")) &&
       identical(head[[2L]], as.name("clgam")) &&
       identical(head[[3L]], as.name("s")))) {
    return(list(
      smooth = list(.clgam_parse_s_call(expr, env = env)),
      linear = list()
    ))
  }
  if (identical(head, as.name("I"))) {
    lab <- paste(deparse(expr), collapse = "")
    return(list(
      smooth = list(),
      linear = list(list(term = lab, term_expr = expr))
    ))
  }
  stop(
    "Unsupported formula term '", paste(deparse(expr), collapse = " "),
    "'. Use s(x1, x2, ndx, bdeg, pord), s(z) / s(z, level = \"coarse\"), ",
    "or a bare name for a linear covariate (optionally inside I()).",
    call. = FALSE
  )
}

#' @keywords internal
.clgam_collect_s <- function(expr, env) {
  .clgam_collect_terms(expr, env)$smooth
}

#' @keywords internal
.clgam_parse_s_call <- function(expr, env) {
  args <- as.list(expr)[-1L]
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  pos <- args[nms == ""]
  named <- args[nms != ""]
  term <- vapply(pos, function(a) {
    if (is.character(a) && length(a) == 1L) a
    else paste(deparse(a), collapse = "")
  }, character(1))
  ndx <- NULL
  k <- NULL
  bdeg <- 3
  pord <- 2
  level <- "fine"
  if (length(named)) {
    extra <- setdiff(names(named), c("k", "ndx", "bdeg", "pord", "level"))
    if (length(extra)) {
      warning(
        "Ignoring s() arguments not used by clgam: ",
        paste(extra, collapse = ", "),
        call. = FALSE
      )
    }
    if (!is.null(named$ndx)) ndx <- eval(named$ndx, envir = env)
    if (!is.null(named$k)) k <- eval(named$k, envir = env)
    if (!is.null(named$bdeg)) bdeg <- eval(named$bdeg, envir = env)
    if (!is.null(named$pord)) pord <- eval(named$pord, envir = env)
    if (!is.null(named$level)) {
      level <- as.character(eval(named$level, envir = env))[1L]
    }
  }
  ndx <- .clgam_resolve_ndx(ndx, k)
  level <- match.arg(level, c("fine", "coarse"))
  list(
    term = term, term_expr = pos,
    ndx = ndx, k = ndx, bdeg = bdeg, pord = pord, level = level
  )
}

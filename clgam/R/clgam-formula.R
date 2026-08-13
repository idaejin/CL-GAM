#' Univariate or spatial P-spline term (formula interface)
#'
#' Used on the right-hand side of a \code{\link{clgam}} formula.
#' \code{s(x1, x2)} is the anisotropic spatial field; \code{s(z)} is a
#' univariate smooth. \code{\link{clgam}} parses the \code{s()} \emph{call}
#' and does not evaluate \pkg{mgcv}'s \code{s}. This helper exists so the
#' syntax is documented and so \code{s} is found if a formula is ever
#' evaluated.
#'
#' @param ... one name (\code{s(z)}) or two (\code{s(x1, x2)})
#' @param k knot count: length 1 for a univariate smooth, length 2 (or
#'   recycled) for space
#' @param level \code{"fine"} (default) or \code{"coarse"} (Case B/C
#'   aggregated-scale smooth)
#' @return An object of class \code{"clgam.smooth.spec"} (not used directly)
#' @seealso \code{\link{clgam}}, \code{\link{simulate_ata_scenarios}}
#' @export
s <- function(..., k = NULL, level = "fine") {
  level <- match.arg(level, c("fine", "coarse"))
  vars <- as.list(substitute(list(...)))[-1L]
  varn <- vapply(vars, function(v) {
    if (is.character(v)) v else paste(deparse(v), collapse = "")
  }, character(1))
  structure(
    list(term = varn, k = k, level = level),
    class = "clgam.smooth.spec"
  )
}

#' Parse a CL-GAM formula \code{y ~ s(x1, x2) + s(z)}
#' @keywords internal
.clgam_parse_formula <- function(formula, data = NULL, env = parent.frame()) {
  if (!inherits(formula, "formula") || length(formula) < 3L) {
    stop("Expected a two-sided formula, e.g. y ~ s(x1, x2) + s(z).",
         call. = FALSE)
  }
  env_data <- .clgam_eval_env(data, env)
  specs <- .clgam_collect_s(formula[[3L]], env = env_data)
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
  knots <- spatial$k
  if (!is.null(knots) && length(knots) == 1L) knots <- c(knots, knots)

  smooth_vars <- NULL
  nl_level <- NULL
  knots_nl <- NULL
  if (length(uni)) {
    smooth_vars <- vector("list", length(uni))
    nl_level <- character(length(uni))
    knots_nl <- rep(NA_real_, length(uni))
    cn <- character(length(uni))
    for (i in seq_along(uni)) {
      ui <- uni[[i]]
      smooth_vars[[i]] <- .clgam_eval_term(ui$term_expr[[1L]], env_data)
      nl_level[i] <- ui$level %||% "fine"
      cn[i] <- ui$term
      if (!is.null(ui$k)) knots_nl[i] <- as.numeric(ui$k)[1L]
    }
    names(smooth_vars) <- cn
    if (all(is.na(knots_nl))) knots_nl <- NULL
  }

  list(
    y = y, x1 = x1, x2 = x2,
    smooth_vars = smooth_vars,
    smooth_level = nl_level,
    knots = knots,
    knots_nl = knots_nl,
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

#' Bind formula smooths; expand length-\code{n_coarse} vectors via partition \code{C}
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

#' Walk \code{+} / \code{s()} on a formula RHS
#' @keywords internal
.clgam_collect_s <- function(expr, env) {
  if (is.null(expr)) return(list())
  if (is.symbol(expr)) {
    stop(
      "Linear terms are not parsed from the formula yet; use s(z) or the ",
      "linear= argument. Found: ", as.character(expr),
      call. = FALSE
    )
  }
  if (!is.call(expr)) {
    stop("Cannot parse formula term: ", paste(deparse(expr), collapse = " "),
         call. = FALSE)
  }
  head <- expr[[1L]]
  if (identical(head, as.name("("))) {
    return(.clgam_collect_s(expr[[2L]], env = env))
  }
  if (identical(head, as.name("+"))) {
    return(c(
      .clgam_collect_s(expr[[2L]], env = env),
      .clgam_collect_s(expr[[3L]], env = env)
    ))
  }
  if (identical(head, as.name("s")) ||
      (is.call(head) && identical(head[[1L]], as.name("::")) &&
       identical(head[[2L]], as.name("clgam")) &&
       identical(head[[3L]], as.name("s")))) {
    return(list(.clgam_parse_s_call(expr, env = env)))
  }
  stop(
    "Unsupported formula term '", paste(deparse(expr), collapse = " "),
    "'. Use s(x1, x2) and s(z) / s(z, level = \"coarse\").",
    call. = FALSE
  )
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
  k <- NULL
  level <- "fine"
  if (length(named)) {
    extra <- setdiff(names(named), c("k", "level"))
    if (length(extra)) {
      warning(
        "Ignoring s() arguments not used by clgam: ",
        paste(extra, collapse = ", "),
        call. = FALSE
      )
    }
    if (!is.null(named$k)) k <- eval(named$k, envir = env)
    if (!is.null(named$level)) {
      level <- as.character(eval(named$level, envir = env))[1L]
    }
  }
  level <- match.arg(level, c("fine", "coarse"))
  list(term = term, term_expr = pos, k = k, level = level)
}

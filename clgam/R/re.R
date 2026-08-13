#' Resolve the optional coarse-level random effect
#'
#' @param re \code{NULL}/\code{"none"} (default) or \code{"coarse"}
#' @return \code{"none"} or \code{"coarse"}
#' @keywords internal
.clgam_resolve_re <- function(re = NULL) {
  if (is.null(re) || isFALSE(re) || identical(re, "none")) {
    return("none")
  }
  if (isTRUE(re) || identical(re, "coarse")) {
    return("coarse")
  }
  stop(
    "re must be NULL/'none' or 'coarse' (iid random effect per coarse unit).",
    call. = FALSE
  )
}

#' Incidence matrix mapping fine units to coarse cells (partition C)
#' @keywords internal
.clgam_Z_re <- function(groups, n_coarse) {
  groups <- as.integer(groups)
  m <- length(groups)
  n_coarse <- as.integer(n_coarse)
  Z <- matrix(0, m, n_coarse)
  Z[cbind(seq_len(m), groups)] <- 1
  Z
}

#' Log-determinant of a symmetric matrix
#'
#' Numerically preferable to \code{log(det(m))}.
#'
#' @param m a symmetric numeric matrix
#' @return a scalar, the log of the determinant
#' @keywords internal
logdetsym <- function(m){
  sum(log(eigen(m, symmetric = TRUE, only.values = TRUE)$values))
}

#' Trace of the product of two matrices
#'
#' This function efficiently computes the trace of the product of two matrices \code{A} and \code{B}.
#'
#' @param A TODO
#' @param B TODO
#'
#' @return TODO
#' @keywords internal
trprod <- function(A, B){
  # Equivalent to sum(diag(A %*% B)), i.e. the trace of A %*% B
  tr <- as.numeric(crossprod(c(A), c(t(B))))
}

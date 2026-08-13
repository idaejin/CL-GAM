#' Trace of the product of two matrices
#'
#' This function efficiently computes the trace of the product of two matrices \code{A} and \code{B}.
#'
#' @param A,B numeric matrices of compatible dimensions
#'
#' @return Scalar \eqn{\mathrm{tr}(A B)}, computed as
#'   \code{sum(A * t(B))} without forming \code{A \%*\% B}.
#' @keywords internal
trprod <- function(A, B){
  # Equivalent to sum(diag(A %*% B)), i.e. the trace of A %*% B
  tr <- as.numeric(crossprod(c(A), c(t(B))))
}

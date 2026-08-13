#' Row-wise or Box Kronecker product of two matrices
#'
#' @param X1,X2 numeric matrices with the same number of rows
#'
#' @return Row-wise (Box) Kronecker product: each row of the result is
#'   the Kronecker product of the corresponding rows of \code{X1} and
#'   \code{X2}. Dimensions \eqn{n \times (c_1 c_2)}.
#' @export

rten2 <- function(X1, X2) {
  vec1 <- matrix(1, 1, ncol(X1))
  vec2 <- matrix(1, 1, ncol(X2))
  kronecker(X1, vec2)*kronecker(vec1, X2)
}
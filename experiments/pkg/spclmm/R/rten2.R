#' Row-wise or Box Kronecker product of two matrices
#'
#' @param X1 TODO
#' @param X2 TODO
#'
#' @return TODO
#' @export

rten2 <- function(X1, X2) {
  vec1 <- matrix(1, 1, ncol(X1))
  vec2 <- matrix(1, 1, ncol(X2))
  kronecker(X1, vec2)*kronecker(vec1, X2)
}
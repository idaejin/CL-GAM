#' Inverse of a two-by-two nonsingular partitioned matrix
#'
#' @param A11 TODO
#' @param A12 TODO
#' @param A21 TODO
#' @param A22 TODO
#'
#' @return TODO
#' @export

inv_bblock2 <- function(A11, A12, A21, A22){

  # Assumptions:
  # A11 and A22 are symmetric matrices
  # A21 = t(A12)

  M1 <- chol2inv(chol(A11))
  M2 <- crossprod(A12, M1)
  M3 <- tcrossprod(M2, A21)

  S22 <- chol2inv(chol(A22 - M3))
  S21 <- - S22 %*% M2
  S11 <- M1 - crossprod(M2, S21)

  S <- bblock2(S11, t(S21), S21, S22)

  output <- list(S11 = S11, S12 = t(S21), S21 = S21, S22 = S22, S = S)
  return(output)
}

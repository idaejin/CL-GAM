#' Inverse of a two-by-two nonsingular partitioned matrix
#'
#' @param A11,A12,A21,A22 blocks of the partitioned matrix
#' @return list with \code{S11}, \code{S12}, \code{S21}, \code{S22}, \code{S}
#' @keywords internal
inv_bblock2 <- function(A11, A12, A21, A22) {
  .safe_inv <- function(A) {
    A <- as.matrix(A)
    tryCatch(
      chol2inv(chol(A)),
      error = function(e) {
        tryCatch(
          solve(A),
          error = function(e2) MASS::ginv(A)
        )
      }
    )
  }

  M1 <- .safe_inv(A11)
  M2 <- crossprod(A12, M1)
  M3 <- tcrossprod(M2, A21)

  S22 <- .safe_inv(A22 - M3)
  S21 <- -S22 %*% M2
  S11 <- M1 - crossprod(M2, S21)

  S <- bblock2(S11, t(S21), S21, S22)
  list(S11 = S11, S12 = t(S21), S21 = S21, S22 = S22, S = S)
}

#' Construct elements of the mixed model representation
#'
#' @inheritParams bbasis
#' @param pord TODO
#' @param decom TODO
#'
#' @return TODO
#' @seealso \code{\link{bbasis}}
#' @export

mm_basis <- function(x, xl, xr, ndx, bdeg, pord, decom = 1){
  Bb <- spclmm::bbasis(x, xl, xr, ndx, bdeg)
  knots <- Bb$knots
  B <- Bb$B
  m <- ncol(B)
  D <- diff(diag(m), differences = pord)
  P.svd <- svd(crossprod(D))      # Equivalent to svd(t(D)%*%D)
  Us <- (P.svd$u)[, 1:(m - pord)] # Non-null eigenvectors part
  d <- (P.svd$d)[1:(m - pord)]    # Non-null eigenvalues
  Z <- B %*% Us                   # Random model matrix
  if (decom == 1) {
    Un <- (P.svd$u)[, -(1:(m - pord))] # Null eigenvectors part
    X <- B %*% Un                 # Fixed model matrix (alt. 1)
  } else if (decom == 2) {
    X <- NULL
    for (i in 0:(pord - 1)) {
      X <- cbind(X, x^i)          # Fixed model matrix (alt. 2)
    }
  }

  output <- list(X = X, Z = Z, d = d, B = B, m = m, D = D, U = P.svd$u, knots = knots)
  return(output)
}

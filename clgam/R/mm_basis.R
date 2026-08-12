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
  Bb <- bbasis(x, xl, xr, ndx, bdeg)
  knots <- Bb$knots
  B <- Bb$B
  m <- ncol(B)
  D <- diff(diag(m), differences = pord)
  # crossprod(D) is symmetric PSD by construction: eigen(symmetric=TRUE)
  # (LAPACK dsyevr) is faster than a general svd() (dgesdd) for this case and
  # gives the same eigenvalues/eigenvectors (U=V, singular values = |eigenvalues|,
  # here already >= 0). Called once per spatial dimension / smooth covariate
  # per fit, so the saving compounds across Monte Carlo replicates.
  P.eig <- eigen(crossprod(D), symmetric = TRUE)
  # eigen() sorts by decreasing eigenvalue, same convention as svd()$d.
  Us <- (P.eig$vectors)[, 1:(m - pord)] # Non-null eigenvectors part
  d <- (P.eig$values)[1:(m - pord)]     # Non-null eigenvalues
  Z <- B %*% Us                   # Random model matrix
  if (decom == 1) {
    Un <- (P.eig$vectors)[, -(1:(m - pord))] # Null eigenvectors part
    X <- B %*% Un                 # Fixed model matrix (alt. 1)
  } else if (decom == 2) {
    X <- NULL
    for (i in 0:(pord - 1)) {
      X <- cbind(X, x^i)          # Fixed model matrix (alt. 2)
    }
  }

  output <- list(X = X, Z = Z, d = d, B = B, m = m, D = D, U = P.eig$vectors, knots = knots)
  return(output)
}

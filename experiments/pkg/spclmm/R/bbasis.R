#' Construct B-spline basis
#'
#' It constructs a B-spline basis of equally-spaced knots from a (numerical) covariate within the function \code{mm_basis}.
#'
#' @param x A numerical vector.
#' @param xl A numerical value that indicates a left boundary for \code{x}.
#' @param xr A numerical value that indicates a right boundary for \code{x}.
#' @param ndx A numerical value that indicates the number of internal intervals (or the number of internal knots minus one).
#' @param bdeg A numerical value that indicates the degree of the B-splines.
#'
#' @return A list that contains the B-spline basis (\code{B}), and the internal knots (\code{knots}).
#'
#' @seealso \code{\link{tpower}}, \code{\link{mm_basis}}
#'
#' @export
#'
#' @examples
#' x <- seq(0, 1, length = 200)
#' xls <- min(x) - 0.001
#' xrs <- max(x) + 0.001
#' ndxs <- 10
#' ## A quadratic B-spline basis
#' B2 <- bbasis(x, xl = xls, xr = xrs, ndx = ndxs, bdeg = 2)$B
#' matplot(x, B2, t = "l", lty = 1, lwd = 2, col = rainbow(ncol(B2)))
#' ## A cubic B-spline basis
#' B3 <- bbasis(x, xl = xls, xr = xrs, ndx = ndxs, bdeg = 3)$B
#' matplot(x, B3, t = "l", lty = 1, lwd = 2, col = rainbow(ncol(B3)))

bbasis <- function(x, xl, xr, ndx, bdeg){
  dx <- (xr - xl)/ndx
  knots <- seq(xl - bdeg*dx, xr + bdeg*dx, by = dx)
  P <- outer(x, knots, spclmm::tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1)*dx^bdeg)
  B <- (-1)^(bdeg + 1)*P %*% t(D)

  output <- list(B = B, knots = knots)
  return(output)
}

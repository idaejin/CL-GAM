#' Truncated p-th power function
#'
#' It constructs a p-th degree truncated power function along a (numerical) vector within the function \code{bspline}. This is an internal function of the package \code{spclmm} used for the construction of a B-spline basis.
#'
#' @param x A numerical vector.
#' @param t A numerical value of truncation point.
#' @param p A numerical value that indicates the degree of the power.
#'
#' @return A numerical vector whose length is equal to \code{length(x)}.
#'
#' @seealso \code{\link{bbasis}}
#'
#' @export
#'
#' @examples
#' x <- seq(0, 1, length = 100)
#' f1 <- tpower(x = x, t = 0.2, p = 2)
#' f2 <- tpower(x = x, t = 0.4, p = 2)
#' f3 <- tpower(x = x, t = 0.6, p = 2)
#' f4 <- tpower(x = x, t = 0.8, p = 2)
#'
#' ## Basis formed by four truncated power functions of degree 2 with equally-spaced knots
#' plot(x, f1, t = "l", col = 2, lwd = 2, main = "Truncated power function basis of degree 2", xlab = "", ylab = "")
#' lines(x, f2, col = 3, lwd = 2)
#' lines(x, f3, col = 4, lwd = 2)
#' lines(x, f4, col = 5, lwd = 2)
#'
#' ## Quadratic B-spline formed as f1 - 3f2 + 3f3 - f4
#' plot(x, f1 - 3*f2 + 3*f3 - f4, t = "l", lwd = 2, main = "A quadratic B-spline constructed from f1, f2, f3, and f4")

tpower <- function(x, t, p){
  (x - t)^p*(x > t)
}

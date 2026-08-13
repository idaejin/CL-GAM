#' Truncated p-th power function
#'
#' Internal helper for B-spline construction (\code{\link{bbasis}}).
#'
#' @param x A numerical vector.
#' @param t A numerical value of truncation point.
#' @param p A numerical value that indicates the degree of the power.
#'
#' @return A numerical vector whose length is equal to \code{length(x)}.
#'
#' @seealso \code{\link{bbasis}}
#' @keywords internal
#' @examples
#' \dontrun{
#' x <- seq(0, 1, length = 100)
#' f1 <- tpower(x = x, t = 0.2, p = 2)
#' f2 <- tpower(x = x, t = 0.4, p = 2)
#' f3 <- tpower(x = x, t = 0.6, p = 2)
#' f4 <- tpower(x = x, t = 0.8, p = 2)
#' plot(x, f1, type = "l", col = 2, lwd = 2, xlab = "", ylab = "")
#' lines(x, f2, col = 3, lwd = 2)
#' lines(x, f3, col = 4, lwd = 2)
#' lines(x, f4, col = 5, lwd = 2)
#' }
tpower <- function(x, t, p){
  (x - t)^p*(x > t)
}

#' logdetsym function
#'
#' @param TODO
#'
#' @return TODO
#' @keywords internal
logdetsym <- function(m){
  # It computes the log of the det of a symmetric matrix
  # Note: Numerical superior to log(det(m))
  erg <- sum(log(eigen(m, symmetric = TRUE, only.values = TRUE)$values))
}

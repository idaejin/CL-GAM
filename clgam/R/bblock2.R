#' Construct matrix from two-by-two blocks
#'
#' This function recovers a partitioned matrix from two-by-two blocks.
#'
#' @param A11,A12,A21,A22 blocks of the partitioned matrix
#'
#' @return The assembled matrix \code{rbind(cbind(A11, A12), cbind(A21, A22))}.
#' @keywords internal
bblock2 <- function(A11, A12, A21, A22){
  block <- rbind(cbind(A11, A12), cbind(A21, A22))
  unname(block, force = FALSE)
}
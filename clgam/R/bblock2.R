#' Construct matrix from two-by-two blocks
#'
#' This function recovers a partitioned matrix from two-by-two blocks.
#'
#' @param A11 TODO
#' @param A12 TODO
#' @param A21 TODO
#' @param A22 TODO
#'
#' @return TODO
#' @keywords internal
bblock2 <- function(A11, A12, A21, A22){
  block <- rbind(cbind(A11, A12), cbind(A21, A22))
  unname(block, force = FALSE)
}
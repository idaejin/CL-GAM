#' Construct matrix from three-by-three blocks
#'
#' This function recovers a partitioned matrix from three-by-three blocks.
#'
#' @param A11,A12,A13,A21,A22,A23,A31,A32,A33 blocks of the partitioned matrix
#'
#' @return The assembled three-by-three block matrix.
#' @keywords internal
bblock3 <- function(A11, A12, A13,
                    A21, A22, A23,
                    A31, A32, A33){

  block <- rbind(cbind(A11, A12, A13),
                 cbind(A21, A22, A23),
                 cbind(A31, A32, A33))

  unname(block, force = FALSE)

}
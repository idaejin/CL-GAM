#' Construct matrix from three-by-three blocks
#'
#' This function recovers a partitioned matrix from three-by-three blocks.
#'
#' @param A11 TODO
#' @param A12 TODO
#' @param A13 TODO
#' @param A21 TODO
#' @param A22 TODO
#' @param A23 TODO
#' @param A31 TODO
#' @param A32 TODO
#' @param A33 TODO
#'
#' @return TODO
#' @keywords internal
bblock3 <- function(A11, A12, A13,
                    A21, A22, A23,
                    A31, A32, A33){

  block <- rbind(cbind(A11, A12, A13),
                 cbind(A21, A22, A23),
                 cbind(A31, A32, A33))

  unname(block, force = FALSE)

}
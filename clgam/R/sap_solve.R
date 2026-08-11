#' Prefer Rcpp SOP kernels when compiled; fall back to pure R.
#' @keywords internal
.sap_solve_schur <- function(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache = NULL) {
  if (isTRUE(getOption("clgam.use_rcpp", TRUE)) &&
      exists("sap_solve_schur_cpp", mode = "function")) {
    A11inv <- if (!is.null(cache$A11inv)) cache$A11inv else matrix(0, 0, 0)
    out <- tryCatch(
      sap_solve_schur_cpp(XtX, ZtX, ZtZ, ZtXtZ, u, G, A11inv),
      error = function(e) NULL
    )
    if (!is.null(out)) {
      return(list(
        b.fixed = as.numeric(out$b.fixed),
        b.random = as.numeric(out$b.random),
        dZtNZ = as.numeric(out$dZtNZ),
        cache = list(A11inv = out$A11inv)
      ))
    }
  }
  .sap_solve_schur_R(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache)
}

#' Pure-R Schur SOP solve (fallback)
#' @keywords internal
.sap_solve_schur_R <- function(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache = NULL) {
  p <- ncol(XtX)
  q <- length(G)
  u1 <- u[seq_len(p)]
  u2 <- u[(p + 1L):(p + q)]

  A21 <- ZtX
  A12 <- t(G * ZtX)
  A22 <- t(G * ZtZ)
  diag(A22) <- diag(A22) + 1

  if (is.null(cache$A11inv)) {
    cache <- list(A11inv = tryCatch(
      solve(XtX),
      error = function(e) MASS::ginv(XtX)
    ))
  }
  A11inv <- cache$A11inv

  A11inv_u1 <- as.numeric(A11inv %*% u1)
  A11inv_A12 <- A11inv %*% A12
  S <- A22 - A21 %*% A11inv_A12
  rhs2 <- u2 - as.numeric(A21 %*% A11inv_u1)

  Sinv <- try(solve(S), silent = TRUE)
  if (inherits(Sinv, "try-error")) {
    Sinv <- MASS::ginv(S)
  }
  b2 <- as.numeric(Sinv %*% rhs2)
  b1 <- as.numeric(A11inv %*% (u1 - as.numeric(A12 %*% b2)))

  H22 <- Sinv
  H21 <- -Sinv %*% (A21 %*% A11inv)
  H_bottom <- cbind(H21, H22)

  dZtNZ <- colSums(t(H_bottom) * ZtXtZ)

  list(
    b.fixed = b1,
    b.random = G * b2,
    dZtNZ = dZtNZ,
    cache = cache
  )
}

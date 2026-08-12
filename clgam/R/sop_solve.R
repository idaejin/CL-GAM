#' Prefer Rcpp SOP kernels when compiled; fall back to pure R.
#'
#' Both this and \code{.sop_solve_schur_R} cache, in addition to
#' \code{A11inv = solve(XtX)}, the G-FREE quantities
#' \code{N = ZtZ - ZtX \%*\% A11inv \%*\% t(ZtX)} and
#' \code{rhs2 = u2 - ZtX \%*\% A11inv \%*\% u1}: within one outer PIRLS
#' iteration only \code{G} (from the variance components) changes across
#' inner SOP iterations, and \code{S = A22 - A21 A11inv A12} algebraically
#' equals \code{N \%*\% diag(G) + I} (verified numerically), so \code{N} and
#' \code{rhs2} need computing only once per outer iteration rather than
#' being rebuilt (at O(q^2 p) each) on every inner iteration.
#' @keywords internal
.sop_solve_schur <- function(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache = NULL) {
  if (isTRUE(getOption("clgam.use_rcpp", TRUE)) &&
      exists("sop_solve_schur_cpp", envir = environment(), inherits = FALSE)) {
    A11inv <- if (!is.null(cache$A11inv)) cache$A11inv else matrix(0, 0, 0)
    N <- if (!is.null(cache$N)) cache$N else matrix(0, 0, 0)
    rhs2 <- if (!is.null(cache$rhs2)) cache$rhs2 else numeric(0)
    out <- tryCatch(
      sop_solve_schur_cpp(XtX, ZtX, ZtZ, ZtXtZ, u, G, A11inv, N, rhs2),
      error = function(e) NULL
    )
    if (!is.null(out)) {
      return(list(
        b.fixed = as.numeric(out$b.fixed),
        b.random = as.numeric(out$b.random),
        dZtNZ = as.numeric(out$dZtNZ),
        cache = list(A11inv = out$A11inv, N = out$N, rhs2 = as.numeric(out$rhs2))
      ))
    }
  }
  .sop_solve_schur_R(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache)
}

#' Pure-R Schur SOP solve (fallback)
#' @keywords internal
.sop_solve_schur_R <- function(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache = NULL) {
  p <- ncol(XtX)
  q <- length(G)
  u1 <- u[seq_len(p)]
  u2 <- u[(p + 1L):(p + q)]

  A21 <- ZtX
  XtZ <- t(ZtX)

  have_cache <- !is.null(cache$A11inv) && !is.null(cache$N) && !is.null(cache$rhs2)
  if (!have_cache) {
    A11inv <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
    N <- ZtZ - A21 %*% (A11inv %*% XtZ)
    rhs2 <- as.numeric(u2 - A21 %*% (A11inv %*% u1))
    cache <- list(A11inv = A11inv, N = N, rhs2 = rhs2)
  }
  A11inv <- cache$A11inv
  N <- cache$N
  rhs2 <- cache$rhs2

  # S = N %*% diag(G) + I (column rescale of the cached, G-free N).
  S <- sweep(N, 2L, G, `*`)
  diag(S) <- diag(S) + 1
  A12 <- sweep(XtZ, 2L, G, `*`)

  # S is generally NOT symmetric (only N is); a plain solve() (LU-based)
  # handles that correctly, unlike a Cholesky/sympd-only method.
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

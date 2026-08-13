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
#'
#' Set \code{options(clgam.sop.backend = "kron_hybrid")} to invert \code{S}
#' by an exact B3-vs-rest block factorization (same estimator). Default
#' \code{"dense"} keeps the compiled or dense-\code{solve} path.
#' @param meta optional \code{.sop_kron_meta()} list; used only when the
#'   backend is \code{"kron_hybrid"}.
#' @keywords internal
.sop_solve_schur <- function(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache = NULL,
                             meta = NULL) {
  backend <- getOption("clgam.sop.backend", "dense")
  use_kron <- identical(backend, "kron_hybrid") &&
    .sop_kron_applicable(meta, length(G))
  # Kronecker-Schur block inverse is R-only (needs meta). Dense keeps Rcpp.
  if (!use_kron &&
      isTRUE(getOption("clgam.use_rcpp", TRUE)) &&
      exists("sop_solve_schur_cpp", mode = "function")) {
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
  .sop_solve_schur_R(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache, meta = meta)
}

#' Pure-R Schur SOP solve (fallback / Kronecker-Schur hybrid)
#' @keywords internal
.sop_solve_schur_R <- function(XtX, ZtX, ZtZ, ZtXtZ, u, G, cache = NULL,
                               meta = NULL) {
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

  backend <- getOption("clgam.sop.backend", "dense")
  if (identical(backend, "kron_hybrid") && .sop_kron_applicable(meta, q)) {
    Sinv <- .sop_Sinv_kron(N, G, meta)
  } else {
    Sinv <- .sop_Sinv_dense(N, G)
  }
  A12 <- sweep(XtZ, 2L, G, `*`)

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

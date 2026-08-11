#' Sparse backends for composition and Kronecker penalties
#'
#' CL-GAM / \code{pois_SOP} already densifies the SOP system in coefficient
#' space; sparse algebra helps for (1) the composition matrix \code{C} and
#' (2) Kronecker P-spline precisions (LMMsolver-style), not for wrapping the
#' dense \code{V+D} block.
#'
#' @name sparse_backends
NULL

#' Default sparse backend: \code{"Matrix"}, \code{"spam"}, or \code{"dense"}
#' @keywords internal
.clgam_sparse_backend <- function(backend = NULL) {
  if (is.null(backend) || identical(backend, "auto")) {
    backend <- getOption("clgam.sparse.backend", "Matrix")
  }
  match.arg(backend, c("Matrix", "spam", "dense", "auto"))
}

#' Coerce Matrix / matrix to spam without relying on missing S4 methods
#' @keywords internal
.as_spam <- function(C) {
  if (inherits(C, "spam")) return(C)
  spam::as.spam(as.matrix(C))
}

#' Coerce spam to dgCMatrix
#' @keywords internal
.spam_to_dgC <- function(C) {
  Matrix::Matrix(as.matrix(C), sparse = TRUE)
}

#' Coerce composition matrix to a sparse (or dense) backend
#'
#' @param C composition matrix
#' @param backend \code{"Matrix"} (dgCMatrix), \code{"spam"}, \code{"dense"},
#'   or \code{"auto"} (option \code{clgam.sparse.backend}, default Matrix)
#' @return matrix-like object suitable for \code{\%*\%}
#' @export
as_comp_C <- function(C, backend = "auto") {
  backend <- .clgam_sparse_backend(backend)
  if (identical(backend, "auto")) backend <- "Matrix"

  if (identical(backend, "dense")) {
    if (inherits(C, "Matrix") || inherits(C, "spam")) {
      return(as.matrix(C))
    }
    return(C)
  }

  if (identical(backend, "Matrix") && inherits(C, "dgCMatrix")) return(C)
  if (identical(backend, "spam") && inherits(C, "spam")) return(C)

  if (inherits(C, "Matrix")) {
    if (identical(backend, "Matrix")) {
      return(as(as(C, "generalMatrix"), "CsparseMatrix"))
    }
    return(.as_spam(C))
  }
  if (inherits(C, "spam")) {
    if (identical(backend, "spam")) return(C)
    return(.spam_to_dgC(C))
  }

  if (!is.matrix(C)) C <- as.matrix(C)
  nc <- ncol(C)
  nz <- sum(C != 0)
  if (nc < 500L || nz >= 0.25 * length(C)) {
    return(C)
  }

  if (identical(backend, "spam")) {
    return(spam::as.spam(C))
  }
  Matrix::Matrix(C, sparse = TRUE)
}

# Keep internal alias used by pois_*
.as_comp_C <- function(C, backend = "auto") as_comp_C(C, backend = backend)

#' Detect 0-1 partition composition (each fine cell → one coarse)
#' @keywords internal
.is_partition_C <- function(C) {
  if (inherits(C, "Matrix")) {
    if (!all(C@x %in% c(0, 1))) return(FALSE)
    return(all(Matrix::colSums(C) == 1))
  }
  if (inherits(C, "spam")) {
    if (!all(C@entries %in% c(0, 1))) return(FALSE)
    return(all(spam::colSums.spam(C) == 1))
  }
  if (!is.matrix(C)) return(FALSE)
  if (!all(C %in% c(0, 1))) return(FALSE)
  all(colSums(C) == 1)
}

#' Fine→coarse group id for partition C (length = ncol(C))
#' @keywords internal
.partition_groups <- function(C) {
  if (inherits(C, "dgCMatrix")) {
    # one nonzero per column ⇒ @i is already the (0-based) row id per column
    if (length(C@x) != ncol(C)) {
      stop("partition C expected one nonzero per column")
    }
    return(as.integer(C@i) + 1L)
  }
  if (inherits(C, "Matrix")) {
    return(.partition_groups(as(C, "dgCMatrix")))
  }
  if (inherits(C, "spam")) {
    return(.partition_groups(.spam_to_dgC(C)))
  }
  max.col(t(C), ties.method = "first")
}

#' C %*% (gamma * A) with Matrix, spam, or partition rowsum
#' @keywords internal
.comp_mul <- function(C, gamma, A, groups = NULL) {
  if (!is.null(groups)) {
    # rowsum aggregates fine rows → coarse in group order 1..n_coarse
    n_coarse <- nrow(C)
    out <- rowsum(gamma * A, group = factor(groups, levels = seq_len(n_coarse)))
    return(unname(as.matrix(out)))
  }
  out <- C %*% (gamma * A)
  if (inherits(out, "Matrix") || inherits(out, "spam")) {
    out <- as.matrix(out)
  }
  out
}

#' Anisotropic 2D P-spline precision (Kronecker), SparseMatrix
#'
#' \deqn{P = \lambda_1 (I_{m_2} \otimes D_1'D_1) + \lambda_2 (D_2'D_2 \otimes I_{m_1})}
#'
#' @param m1,m2 basis sizes along each axis
#' @param lambda1,lambda2 precision weights (e.g. \code{1/tau^2})
#' @param pord difference order
#' @param backend \code{"Matrix"} (default) or \code{"spam"}
#' @return sparse precision of size \code{(m1*m2) x (m1*m2)}
#' @export
sparse_P2D <- function(m1, m2, lambda1 = 1, lambda2 = 1, pord = 2L,
                       backend = c("Matrix", "spam")) {
  backend <- match.arg(backend)
  if (identical(backend, "Matrix")) {
    D1 <- Matrix::Matrix(diff(diag(m1), differences = pord), sparse = TRUE)
    D2 <- Matrix::Matrix(diff(diag(m2), differences = pord), sparse = TRUE)
    P1 <- Matrix::crossprod(D1)
    P2 <- Matrix::crossprod(D2)
    P <- lambda1 * kronecker(Matrix::Diagonal(m2), P1) +
      lambda2 * kronecker(P2, Matrix::Diagonal(m1))
    return(P)
  }
  P <- sparse_P2D(m1, m2, lambda1, lambda2, pord, backend = "Matrix")
  .as_spam(P)
}

#' Sparse Cholesky solve for Kronecker precision + diagonal weights (demo/API)
#'
#' @param P sparse precision from \code{sparse_P2D}
#' @param w positive diagonal weights (length \code{nrow(P)})
#' @param b right-hand side
#' @return solution vector
#' @keywords internal
sparse_chol_solve <- function(P, w, b) {
  if (inherits(P, "spam")) {
    A <- P + spam::diag.spam(w, nrow = length(w), ncol = length(w))
    return(as.numeric(spam::solve.spam(A, b)))
  }
  A <- P + Matrix::Diagonal(x = w)
  ch <- Matrix::Cholesky(A, perm = TRUE, LDL = FALSE, super = TRUE)
  as.numeric(Matrix::solve(ch, b))
}

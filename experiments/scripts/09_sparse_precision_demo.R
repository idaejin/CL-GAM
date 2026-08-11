# Sparse 2D P-spline precision (Kronecker) — demo, not full CLMM
#
# Why this matters: Ayma's mm_basis uses an SVD null-space reparameterization
# that densifies the random design Z. LMMsolver / Boer (2023) keep a sparse
# precision on B-spline coefficients instead.
#
# This script times forming/factoring a sparse anisotropic penalty
#   P = λ1 (I ⊗ D'D) + λ2 (D'D ⊗ I)
# vs a dense (p+q) solve of the same nominal size — the regime of pois_SAP.

suppressPackageStartupMessages({
  library(Matrix)
})

.sparse_diff_penalty <- function(m, pord = 2L) {
  # D (m-pord) x m first-difference order pord; return D'D sparse
  D <- diff(Diagonal(m), differences = pord)
  crossprod(D)
}

.kronecker_aniso_precision <- function(m1, m2, lambda1, lambda2, pord = 2L) {
  P1 <- .sparse_diff_penalty(m1, pord)
  P2 <- .sparse_diff_penalty(m2, pord)
  I1 <- Diagonal(m1)
  I2 <- Diagonal(m2)
  # θ arranged as vec of m1 x m2 array (column-major): (I2 ⊗ P1) and (P2 ⊗ I1)
  lambda1 * kronecker(I2, P1) + lambda2 * kronecker(P2, I1)
}

m1 <- m2 <- 23L # ~ ndx=20 + bdeg
q <- m1 * m2
cat(sprintf("Basis size m1*m2 = %d\n", q))

set.seed(1)
lambda1 <- 1 / 0.15
lambda2 <- 1 / 1.9

t_sp <- system.time({
  P <- .kronecker_aniso_precision(m1, m2, lambda1, lambda2)
  # Fake sparse normal equations: B'WB + P ≈ diag(w) + P for illustration
  w <- runif(q, 0.5, 2)
  A <- P + Diagonal(x = w)
  ch <- Cholesky(A, perm = TRUE, LDL = FALSE, super = TRUE)
  x <- solve(ch, rnorm(q))
})
cat(sprintf("Sparse Matrix::Cholesky solve: %.3fs  nnz(A)=%d  density=%.4f\n",
            t_sp["elapsed"], nnzero(A), nnzero(A) / q^2))

t_dn <- system.time({
  Ad <- as.matrix(A)
  xd <- solve(Ad, rnorm(q))
})
cat(sprintf("Dense solve same matrix:       %.3fs\n", t_dn["elapsed"]))

# Larger grid like LMMsolver examples (nseg=41)
m1b <- m2b <- 44L
qb <- m1b * m2b
t_big <- system.time({
  Pb <- .kronecker_aniso_precision(m1b, m2b, lambda1, lambda2)
  Ab <- Pb + Diagonal(x = runif(qb, 0.5, 2))
  chb <- Cholesky(Ab, perm = TRUE, LDL = FALSE, super = TRUE)
  xb <- solve(chb, rnorm(qb))
})
cat(sprintf("Sparse Cholesky m=%d (nseg~41): %.3fs  nnz=%d  density=%.5f\n",
            qb, t_big["elapsed"], nnzero(Ab), nnzero(Ab) / qb^2))

cat("\nTakeaway: SparseMatrix wins when we KEEP Kronecker precision on\n")
cat("B-spline coeffs (LMMsolver style). Densified SOP (V+D) is already dense;\n")
cat("wrapping it in Matrix::sparseMatrix does not help.\n")

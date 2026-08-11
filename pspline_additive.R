# ======================================================================
# P-spline additive model with identifiability (sum-to-zero) constraints
#
# Model:      y_i = beta_0 + f_1(x_i1) + ... + f_d(x_id) + eps_i
#             f_j(x) = B_j(x) %*% beta_j          (B-spline basis)
#
# Penalty:    lambda_j * || D_2 beta_j ||^2       (Eilers & Marx)
#
# Constraint: 1' B_j beta_j = 0  <=>  sum_i fhat_j(x_ij) = 0
#             absorbed into the penalty as the rank-1 term
#
#                 kappa * B_j' 1 1' B_j           (p x p, as in the notes)
#
# Solve:      ( B'B + blockdiag(lambda_j D'D + kappa B_j'11'B_j) ) beta = B'y
# ======================================================================

library(splines)

## ---------------------------------------------------------------------
## Building blocks
## ---------------------------------------------------------------------

# B-spline design matrix on the range of x (equally spaced inner knots)
bspline_basis <- function(x, n_basis = 20, degree = 3) {
  xl <- min(x); xr <- max(x)
  n_inner <- n_basis - degree - 1
  inner <- seq(xl, xr, length.out = n_inner + 2)[-c(1, n_inner + 2)]
  knots <- c(rep(xl, degree + 1), inner, rep(xr, degree + 1))
  splineDesign(knots, x, ord = degree + 1, outer.ok = TRUE)
}

# D'D for q-th order difference matrix D
diff_penalty <- function(p, order = 2) {
  D <- diff(diag(p), differences = order)
  crossprod(D)                      # t(D) %*% D
}

## ---------------------------------------------------------------------
## Fit: penalized LS with centering constraints in the penalty
## ---------------------------------------------------------------------
fit_pspline_additive <- function(y, X, n_basis = 20, degree = 3,
                                 diff_order = 2, lambdas = NULL,
                                 kappa = 1e8) {
  n <- length(y); d <- ncol(X)
  if (is.null(lambdas)) lambdas <- rep(10, d)
  ones <- rep(1, n)

  # Per-covariate bases and full design [1 | B_1 | ... | B_d]
  bases <- lapply(seq_len(d), function(j) bspline_basis(X[, j], n_basis, degree))
  Bfull <- cbind(ones, do.call(cbind, bases))

  # Block-diagonal penalty: 0 for intercept, then per smooth
  #   lambda_j * D'D  +  kappa * (B_j' 1)(1' B_j)   <- the trick on the paper
  ptot <- 1 + d * n_basis
  P <- matrix(0, ptot, ptot)
  DtD <- diff_penalty(n_basis, diff_order)
  for (j in seq_len(d)) {
    idx <- (1 + (j - 1) * n_basis + 1):(1 + j * n_basis)
    Bt1 <- crossprod(bases[[j]], ones)            # p x 1
    P[idx, idx] <- lambdas[j] * DtD + kappa * tcrossprod(Bt1)   # p x p
  }

  # Normal equations: (B'B + P) beta = B'y
  beta <- solve(crossprod(Bfull) + P, crossprod(Bfull, y))

  beta0  <- beta[1]
  blocks <- lapply(seq_len(d), function(j)
    beta[(1 + (j - 1) * n_basis + 1):(1 + j * n_basis)])
  fhat <- lapply(seq_len(d), function(j) drop(bases[[j]] %*% blocks[[j]]))
  mu <- beta0 + Reduce(`+`, fhat)

  list(beta0 = beta0, blocks = blocks, bases = bases,
       fhat = fhat, mu = mu, n_basis = n_basis, degree = degree)
}

## ---------------------------------------------------------------------
## Demo: two smooths
## ---------------------------------------------------------------------
set.seed(1)
n  <- 400
x1 <- runif(n); x2 <- runif(n)

f1_true <- sin(2 * pi * x1)                    # ~centered already
f2_true <- exp(2 * x2) - (exp(2) - 1) / 2      # deliberately not centered
beta0_true <- 3
y <- beta0_true + f1_true + f2_true + rnorm(n, 0, 0.3)

fit <- fit_pspline_additive(y, cbind(x1, x2), lambdas = c(5, 5))

## --- check the identifiability constraints  1' B_j beta_j = 0 ---
for (j in 1:2)
  cat(sprintf("sum_i fhat_%d(x_i) = % .3e   (should be ~0)\n",
              j, sum(fit$fhat[[j]])))
cat(sprintf("intercept beta0_hat = %.3f\n", fit$beta0))
cat(sprintf("mean(y)             = %.3f   (beta0 absorbs the mean)\n", mean(y)))
cat(sprintf("RMSE of fit         = %.3f\n", sqrt(mean((y - fit$mu)^2))))

## --- plot: estimated vs true (centered) component functions ---
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
truths <- list(f1_true - mean(f1_true), f2_true - mean(f2_true))
xs     <- list(x1, x2)
for (j in 1:2) {
  o <- order(xs[[j]])
  plot(xs[[j]][o], truths[[j]][o], type = "l", lty = 2, lwd = 1.5,
       xlab = bquote(x[.(j)]), ylab = bquote(hat(f)[.(j)]),
       main = sprintf("Smooth %d:  sum fhat = %.1e", j, sum(fit$fhat[[j]])))
  lines(xs[[j]][o], fit$fhat[[j]][o], col = "orange", lwd = 2)
  abline(h = 0, col = "grey")
  legend("topright", c("true (centered)", "P-spline estimate"),
         lty = c(2, 1), lwd = c(1.5, 2), col = c("black", "orange"), bty = "n")
}
par(op)

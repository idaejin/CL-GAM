# ======================================================================
# P-spline additive model with sum-to-zero constraints + confidence bands
#
# Model:      y_i = beta_0 + f_1(x_i1) + ... + f_d(x_id) + eps_i
# Penalty:    lambda_j * ||D_2 beta_j||^2
# Constraint: 1' B_j beta_j = 0, absorbed into the penalty as
#             kappa * B_j' 1 1' B_j              (rank-1, p x p)
#
# Inference:  beta_hat = (B'B + P)^{-1} B'y  =: A^{-1} B'y
#   Bayesian covariance (Wahba/Silverman, used by mgcv):
#       V_bayes = sigma^2 * A^{-1}
#   Frequentist "sandwich" covariance (ignores smoothing bias):
#       V_freq  = sigma^2 * A^{-1} B'B A^{-1}
#   sigma^2_hat = RSS / (n - edf),  edf = tr( B A^{-1} B' ) = tr(A^{-1} B'B)
#
#   Pointwise band for component j at new x:
#       fhat_j(x) +/- z * sqrt( b_j(x)' V[jj] b_j(x) )
# ======================================================================

library(splines)

## ---------------------------------------------------------------------
## Building blocks
## ---------------------------------------------------------------------
make_knots <- function(x, n_basis = 20, degree = 3) {
  xl <- min(x); xr <- max(x)
  n_inner <- n_basis - degree - 1
  inner <- seq(xl, xr, length.out = n_inner + 2)[-c(1, n_inner + 2)]
  c(rep(xl, degree + 1), inner, rep(xr, degree + 1))
}

bspline_basis <- function(x, knots, degree = 3)
  splineDesign(knots, x, ord = degree + 1, outer.ok = TRUE)

diff_penalty <- function(p, order = 2) {
  D <- diff(diag(p), differences = order)
  crossprod(D)
}

## ---------------------------------------------------------------------
## Fit + covariance
## ---------------------------------------------------------------------
fit_pspline_additive <- function(y, X, n_basis = 20, degree = 3,
                                 diff_order = 2, lambdas = NULL,
                                 kappa = 1e8, cov_type = c("bayes", "freq")) {
  cov_type <- match.arg(cov_type)
  n <- length(y); d <- ncol(X)
  if (is.null(lambdas)) lambdas <- rep(10, d)
  ones <- rep(1, n)

  knots <- lapply(seq_len(d), function(j) make_knots(X[, j], n_basis, degree))
  bases <- lapply(seq_len(d), function(j) bspline_basis(X[, j], knots[[j]], degree))
  Bfull <- cbind(ones, do.call(cbind, bases))

  # penalty blocks: lambda_j D'D + kappa (B_j'1)(1'B_j)
  ptot <- 1 + d * n_basis
  P <- matrix(0, ptot, ptot)
  DtD <- diff_penalty(n_basis, diff_order)
  for (j in seq_len(d)) {
    idx <- (1 + (j - 1) * n_basis + 1):(1 + j * n_basis)
    Bt1 <- crossprod(bases[[j]], ones)
    P[idx, idx] <- lambdas[j] * DtD + kappa * tcrossprod(Bt1)
  }

  BtB  <- crossprod(Bfull)
  A    <- BtB + P
  Ainv <- solve(A)
  beta <- drop(Ainv %*% crossprod(Bfull, y))

  # fitted values, effective df, error variance
  mu   <- drop(Bfull %*% beta)
  edf  <- sum(diag(Ainv %*% BtB))                 # tr(hat matrix)
  rss  <- sum((y - mu)^2)
  sigma2 <- rss / (n - edf)

  # coefficient covariance
  V <- if (cov_type == "bayes") sigma2 * Ainv
       else                     sigma2 * Ainv %*% BtB %*% Ainv

  blocks_idx <- lapply(seq_len(d), function(j)
    (1 + (j - 1) * n_basis + 1):(1 + j * n_basis))
  fhat <- lapply(seq_len(d), function(j)
    drop(bases[[j]] %*% beta[blocks_idx[[j]]]))

  list(beta = beta, beta0 = beta[1], V = V, sigma2 = sigma2, edf = edf,
       knots = knots, degree = degree, blocks_idx = blocks_idx,
       fhat = fhat, mu = mu, d = d)
}

## Component-wise prediction with pointwise CI (and Bonferroni-ish option)
predict_component <- function(fit, j, xgrid, level = 0.95) {
  Bg  <- bspline_basis(xgrid, fit$knots[[j]], fit$degree)
  idx <- fit$blocks_idx[[j]]
  f   <- drop(Bg %*% fit$beta[idx])
  se  <- sqrt(pmax(0, rowSums((Bg %*% fit$V[idx, idx]) * Bg)))  # diag(Bg V Bg')
  z   <- qnorm(1 - (1 - level) / 2)
  data.frame(x = xgrid, fit = f, se = se,
             lwr = f - z * se, upr = f + z * se)
}

## ---------------------------------------------------------------------
## Demo
## ---------------------------------------------------------------------
set.seed(1)
n  <- 400
x1 <- runif(n); x2 <- runif(n)
f1_true <- sin(2 * pi * x1)
f2_true <- exp(2 * x2) - (exp(2) - 1) / 2
y <- 3 + f1_true + f2_true + rnorm(n, 0, 0.3)

fit <- fit_pspline_additive(y, cbind(x1, x2), lambdas = c(5, 5))

for (j in 1:2)
  cat(sprintf("sum_i fhat_%d(x_i) = % .3e\n", j, sum(fit$fhat[[j]])))
cat(sprintf("beta0_hat = %.3f | mean(y) = %.3f\n", fit$beta0, mean(y)))
cat(sprintf("edf = %.2f | sigma_hat = %.3f (true 0.3)\n",
            fit$edf, sqrt(fit$sigma2)))

## empirical coverage check of the 95%% bands at the observed points
truths <- list(f1_true - mean(f1_true), f2_true - mean(f2_true))
xs <- list(x1, x2)
for (j in 1:2) {
  pr <- predict_component(fit, j, xs[[j]])
  cov_j <- mean(truths[[j]] >= pr$lwr & truths[[j]] <= pr$upr)
  cat(sprintf("coverage of 95%% band, smooth %d: %.1f%%\n", j, 100 * cov_j))
}

## plot with shaded bands
png("pspline_ci.png", width = 1100, height = 420, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (j in 1:2) {
  g  <- seq(min(xs[[j]]), max(xs[[j]]), length.out = 200)
  pr <- predict_component(fit, j, g)
  plot(NA, xlim = range(g), ylim = range(pr$lwr, pr$upr, truths[[j]]),
       xlab = bquote(x[.(j)]), ylab = bquote(hat(f)[.(j)]),
       main = sprintf("Smooth %d with 95%% pointwise CI", j))
  polygon(c(g, rev(g)), c(pr$lwr, rev(pr$upr)),
          col = adjustcolor("orange", 0.30), border = NA)
  lines(g, pr$fit, col = "orange3", lwd = 2)
  o <- order(xs[[j]])
  lines(xs[[j]][o], truths[[j]][o], lty = 2)
  abline(h = 0, col = "grey")
  legend("topright", c("estimate", "95% CI", "true (centered)"),
         lwd = c(2, 8, 1), lty = c(1, 1, 2),
         col = c("orange3", adjustcolor("orange", .3), "black"), bty = "n")
}
par(op); dev.off()

# ======================================================================
# PS-ANOVA: smooth-ANOVA decomposition with tensor-product P-splines
# (Lee & Durban style), constraints imposed by the penalty trick.
#
#   y = beta0 + f1(x1) + f2(x2) + f12(x1,x2) + eps
#
# Bases:
#   main effects : B1 (n x p1), B2 (n x p2)
#   interaction  : B12 = B1 box B2  (row-wise Kronecker, n x p1p2)
#
# Roughness penalties:
#   f1 : lambda1 * D'D            f2 : lambda2 * D'D
#   f12: lam12_1 * (D1'D1 x I) + lam12_2 * (I x D2'D2)   (anisotropic)
#
# Identifiability (generalizing 1'B beta = 0 from the notes):
#   f1, f2   : orthogonal to the constant  ->  kappa * B_j'11'B_j
#   f12      : orthogonal to constant AND both main-effect spaces.
#              Let Z = [1 | B1 | B2]  (n x (1+p1+p2)). Constraint
#              Z' B12 beta12 = 0, absorbed as   kappa * B12' Z Z' B12.
#              (Same trick: replace the vector 1 by the matrix Z.)
# ======================================================================

suppressMessages({library(splines); library(mgcv)})

make_knots <- function(x, p, deg = 3) {
  xl <- min(x); xr <- max(x)
  inner <- seq(xl, xr, length.out = p - deg + 1)[-c(1, p - deg + 1)]
  c(rep(xl, deg + 1), inner, rep(xr, deg + 1))
}
bbase <- function(x, knots, deg = 3) splineDesign(knots, x, ord = deg + 1, outer.ok = TRUE)
DtDm  <- function(p, q = 2) crossprod(diff(diag(p), differences = q))
rowkron <- function(A, B) {                     # row-wise Kronecker (box product)
  A[, rep(seq_len(ncol(A)), each = ncol(B))] *
  B[, rep(seq_len(ncol(B)), times = ncol(A))]
}

## ---- simulate data with a genuine interaction ---------------------------
set.seed(2)
n  <- 1500
x1 <- runif(n); x2 <- runif(n)
f1_t  <- sin(2 * pi * x1)
f2_t  <- 0.5 * (exp(2 * x2) - (exp(2) - 1) / 2)
f12_t <- 1.5 * sin(2 * pi * x1) * (x2 - 0.5)
y <- 3 + f1_t + f2_t + f12_t + rnorm(n, 0, 0.3)

## ---- bases ---------------------------------------------------------------
p1 <- 12; p2 <- 12; deg <- 3
k1 <- make_knots(x1, p1, deg); k2 <- make_knots(x2, p2, deg)
B1  <- bbase(x1, k1, deg)
B2  <- bbase(x2, k2, deg)
B12 <- rowkron(B1, B2)                          # n x (p1*p2)

## ---- penalties -------------------------------------------------------------
kappa <- 1e8
P1 <- DtDm(p1); P2 <- DtDm(p2)
C1 <- tcrossprod(crossprod(B1, rep(1, n)))      # B1'11'B1
C2 <- tcrossprod(crossprod(B2, rep(1, n)))      # B2'11'B2

S12a <- kronecker(P1, diag(p2))                 # smooth in x1 direction
S12b <- kronecker(diag(p1), P2)                 # smooth in x2 direction
Z    <- cbind(1, B1, B2)                        # lower-order model space
K    <- crossprod(B12, Z)                       # (p1p2) x (1+p1+p2)
C12  <- tcrossprod(K)                           # B12' Z Z' B12

## ---- fit via gam + paraPen (REML selects the 4 lambdas) --------------------
fit <- gam(y ~ B1 + B2 + B12,
           paraPen = list(
             B1  = list(P1, C1,        sp = c(-1, kappa)),
             B2  = list(P2, C2,        sp = c(-1, kappa)),
             B12 = list(S12a, S12b, C12, sp = c(-1, -1, kappa))),
           method = "REML")

cat(sprintf("REML lambdas: f1 %.2f | f2 %.2f | f12 (%.2f, %.2f)\n",
            fit$sp[1], fit$sp[2], fit$sp[3], fit$sp[4]))
cat(sprintf("edf = %.1f | sigma_hat = %.3f (true 0.3)\n",
            sum(fit$edf), sqrt(fit$sig2)))

## ---- verify all identifiability constraints --------------------------------
b   <- coef(fit)
# coefficient blocks follow the formula order: intercept, B1, B2, B12
i1  <- 1 + seq_len(p1)
i2  <- 1 + p1 + seq_len(p2)
i12 <- 1 + p1 + p2 + seq_len(p1 * p2)
fh1 <- drop(B1 %*% b[i1]); fh2 <- drop(B2 %*% b[i2]); fh12 <- drop(B12 %*% b[i12])

cat(sprintf("sum f1hat            = % .2e\n", sum(fh1)))
cat(sprintf("sum f2hat            = % .2e\n", sum(fh2)))
cat(sprintf("max |Z' f12hat|      = % .2e   (orthogonal to 1, B1, B2)\n",
            max(abs(crossprod(Z, fh12)))))
cat(sprintf("intercept = %.3f | mean(y) = %.3f\n", b[1], mean(y)))

## ---- accuracy vs truth and vs mgcv's built-in ti() decomposition -----------
mu <- fitted(fit)
cat(sprintf("RMSE(mu, true surface) = %.3f\n",
            sqrt(mean((mu - (3 + f1_t + f2_t + f12_t))^2))))
ref <- gam(y ~ ti(x1, bs = "ps", k = p1) + ti(x2, bs = "ps", k = p2) +
               ti(x1, x2, bs = c("ps", "ps"), k = c(p1, p2)), method = "REML")
cat(sprintf("cor with ti() fit = %.5f | RMS diff = %.3f\n",
            cor(mu, fitted(ref)), sqrt(mean((mu - fitted(ref))^2))))

## ---- plot: main effects with CI + interaction surface ----------------------
se_band <- function(Bg, idx, level = .95) {
  f  <- drop(Bg %*% b[idx])
  se <- sqrt(pmax(0, rowSums((Bg %*% fit$Vp[idx, idx]) * Bg)))
  z  <- qnorm(1 - (1 - level)/2); cbind(f, f - z*se, f + z*se)
}
png("pspline_psanova.png", width = 1400, height = 430, res = 110)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
g <- seq(0.01, 0.99, length.out = 60)
for (j in 1:2) {
  Bg <- bbase(g, if (j == 1) k1 else k2, deg)
  pr <- se_band(Bg, if (j == 1) i1 else i2)
  tr <- if (j == 1) sin(2*pi*g) - mean(f1_t) else
        0.5*(exp(2*g) - (exp(2)-1)/2) - mean(f2_t)
  plot(NA, xlim = 0:1, ylim = range(pr, tr), xlab = bquote(x[.(j)]),
       ylab = bquote(hat(f)[.(j)]), main = sprintf("Main effect %d", j))
  polygon(c(g, rev(g)), c(pr[,2], rev(pr[,3])),
          col = adjustcolor("orange", .3), border = NA)
  lines(g, pr[,1], col = "orange3", lwd = 2); lines(g, tr, lty = 2)
  abline(h = 0, col = "grey")
}
G   <- expand.grid(x1 = g, x2 = g)
B12g <- rowkron(bbase(G$x1, k1, deg), bbase(G$x2, k2, deg))
F12  <- matrix(B12g %*% b[i12], length(g), length(g))
image(g, g, F12, col = hcl.colors(40, "RdBu", rev = TRUE),
      xlab = expression(x[1]), ylab = expression(x[2]),
      main = expression(hat(f)[12](x[1], x[2])))
contour(g, g, F12, add = TRUE, lwd = .5)
par(op); dev.off()

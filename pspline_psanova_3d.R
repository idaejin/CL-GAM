# ======================================================================
# PS-ANOVA in THREE dimensions, constraints by the penalty trick.
#
#  y = b0 + f1 + f2 + f3 + f12 + f13 + f23 + f123 + eps
#
# Bases:  B_j marginal;  B_jk = B_j box B_k;  B_123 = B_1 box B_2 box B_3
# Coefs:  beta_j (mains), theta_jk (2-way), theta_123 (3-way)
#
# Roughness: each term gets one difference penalty PER DIRECTION
#   f_j   : lam * D'D
#   f_jk  : lam_a (D'D x I) + lam_b (I x D'D)
#   f_123 : lam_a (D'D x I x I) + lam_b (I x D'D x I) + lam_c (I x I x D'D)
#   -> 3 + 6 + 3 = 12 REML-estimated smoothing parameters.
#
# Identifiability, RECURSIVE: each term orthogonal to ALL terms below it.
#   C_j   = 1                              -> kappa * B_j'11'B_j
#   C_jk  = [1 | B_j | B_k]                -> kappa * B_jk' C C' B_jk
#   C_123 = [1 | B1|B2|B3 | B12|B13|B23]   -> kappa * B123' C C' B123
# ======================================================================
suppressMessages({library(splines); library(mgcv)})

make_knots <- function(x, p, deg = 3) {
  xl <- min(x); xr <- max(x)
  inner <- seq(xl, xr, length.out = p - deg + 1)[-c(1, p - deg + 1)]
  c(rep(xl, deg + 1), inner, rep(xr, deg + 1))
}
bbase   <- function(x, k, deg = 3) splineDesign(k, x, ord = deg+1, outer.ok = TRUE)
DtDm    <- function(p, q = 2) crossprod(diff(diag(p), differences = q))
rowkron <- function(A, B) A[, rep(seq_len(ncol(A)), each  = ncol(B))] *
                          B[, rep(seq_len(ncol(B)), times = ncol(A))]
Kcon    <- function(Bterm, C) tcrossprod(crossprod(Bterm, C))  # B'CC'B

## ---- data with genuine 2- and 3-way interactions ---------------------------
set.seed(4)
n  <- 1000
x1 <- runif(n); x2 <- runif(n); x3 <- runif(n)
f1_t   <- sin(2*pi*x1)
f2_t   <- 0.5*(exp(2*x2) - (exp(2)-1)/2)
f3_t   <- 2*(x3 - 0.5)^2 - 1/6
f12_t  <- 1.5*sin(2*pi*x1)*(x2 - 0.5)
f13_t  <- 0
f23_t  <- (x2 - 0.5)*(x3 - 0.5)*2
f123_t <- 2*sin(2*pi*x1)*(x2 - 0.5)*(x3 - 0.5)*4
y <- 3 + f1_t + f2_t + f3_t + f12_t + f23_t + f123_t + rnorm(n, 0, 0.3)

## ---- bases ------------------------------------------------------------------
p <- 6; deg <- 3
k1 <- make_knots(x1,p,deg); k2 <- make_knots(x2,p,deg); k3 <- make_knots(x3,p,deg)
B1 <- bbase(x1,k1); B2 <- bbase(x2,k2); B3 <- bbase(x3,k3)
B12  <- rowkron(B1,B2); B13 <- rowkron(B1,B3); B23 <- rowkron(B2,B3)
B123 <- rowkron(B12, B3)                     # n x p^3

## ---- roughness penalties ------------------------------------------------------
P <- DtDm(p); I <- diag(p)
S12a <- kronecker(P,I); S12b <- kronecker(I,P)      # also serves 13, 23
S3a  <- kronecker(P, kronecker(I,I))
S3b  <- kronecker(I, kronecker(P,I))
S3c  <- kronecker(I, kronecker(I,P))

## ---- constraint penalties (recursive C's) -------------------------------------
kappa <- 1e8
one <- rep(1, n)
K1 <- Kcon(B1, one); K2 <- Kcon(B2, one); K3 <- Kcon(B3, one)
C12  <- cbind(1, B1, B2);  K12  <- Kcon(B12,  C12)
C13  <- cbind(1, B1, B3);  K13  <- Kcon(B13,  C13)
C23  <- cbind(1, B2, B3);  K23  <- Kcon(B23,  C23)
C123 <- cbind(1, B1, B2, B3, B12, B13, B23)
K123 <- Kcon(B123, C123)

## ---- fit -----------------------------------------------------------------------
t0 <- proc.time()
fit <- gam(y ~ B1 + B2 + B3 + B12 + B13 + B23 + B123,
  paraPen = list(
    B1   = list(P, K1, sp = c(-1, kappa)),
    B2   = list(P, K2, sp = c(-1, kappa)),
    B3   = list(P, K3, sp = c(-1, kappa)),
    B12  = list(S12a, S12b, K12, sp = c(-1, -1, kappa)),
    B13  = list(S12a, S12b, K13, sp = c(-1, -1, kappa)),
    B23  = list(S12a, S12b, K23, sp = c(-1, -1, kappa)),
    B123 = list(S3a, S3b, S3c, K123, sp = c(-1, -1, -1, kappa))),
  method = "REML")
cat(sprintf("fit time %.1f s | edf = %.1f | sigma_hat = %.3f (true 0.3)\n",
            (proc.time()-t0)[3], sum(fit$edf), sqrt(fit$sig2)))

## ---- coefficient blocks (formula order) ------------------------------------------
b <- coef(fit); off <- 1
blk <- function(m) { i <- off + seq_len(m); off <<- off + m; i }
i1 <- blk(p); i2 <- blk(p); i3 <- blk(p)
i12 <- blk(p^2); i13 <- blk(p^2); i23 <- blk(p^2); i123 <- blk(p^3)

## ---- verify the full constraint hierarchy -----------------------------------------
chk <- function(Bt, idx, C) max(abs(crossprod(C, Bt %*% b[idx])))
cat(sprintf("mains   |1'f|      : %.1e  %.1e  %.1e\n",
            chk(B1,i1,one), chk(B2,i2,one), chk(B3,i3,one)))
cat(sprintf("2-way   |C'f|      : %.1e  %.1e  %.1e\n",
            chk(B12,i12,C12), chk(B13,i13,C13), chk(B23,i23,C23)))
cat(sprintf("3-way   |C123'f|   : %.1e\n", chk(B123,i123,C123)))
cat(sprintf("intercept = %.3f | mean(y) = %.3f\n", b[1], mean(y)))

## ---- accuracy ------------------------------------------------------------------------
tru <- 3 + f1_t + f2_t + f3_t + f12_t + f23_t + f123_t
cat(sprintf("RMSE(mu, truth) = %.3f\n", sqrt(mean((fitted(fit) - tru)^2))))
## component-level: is the (true-zero) f13 estimated as ~0? REML lambda should be huge
cat(sprintf("||f13hat|| / ||f12hat|| = %.3f  (f13 true = 0)\n",
            sqrt(mean((B13 %*% b[i13])^2)) / sqrt(mean((B12 %*% b[i12])^2))))
cat(sprintf("REML lambdas (12): %s\n", paste(sprintf("%.2g", fit$sp), collapse=" ")))

## ---- plot: two slices of the estimated 3-way interaction ------------------------------
png("pspline_psanova_3d.png", width = 1100, height = 430, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
g <- seq(0.02, 0.98, length.out = 50)
G <- expand.grid(x1 = g, x2 = g)
for (x3v in c(0.15, 0.85)) {
  Bs <- rowkron(rowkron(bbase(G$x1,k1), bbase(G$x2,k2)),
                bbase(rep(x3v, nrow(G)), k3))
  F3 <- matrix(Bs %*% b[i123], length(g), length(g))
  image(g, g, F3, col = hcl.colors(40, "RdBu", rev = TRUE),
        xlab = expression(x[1]), ylab = expression(x[2]),
        main = bquote(hat(f)[123] ~ "at" ~ x[3] == .(x3v)))
  contour(g, g, F3, add = TRUE, lwd = .5)
}
par(op); dev.off()

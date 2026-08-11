#!/usr/bin/env Rscript
# Benchmark composition backends: dense / Matrix / spam (+ partition rowsum).
# From experiments/: Rscript scripts/13_sparse_backend_bench.R
root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))

dat <- clgam_load_madrid(FALSE)
cat("spclmm", as.character(packageVersion("spclmm")), "\n")

# Kronecker P demo
P_M <- sparse_P2D(23, 23, 1 / 0.15, 1 / 1.9, backend = "Matrix")
P_S <- sparse_P2D(23, 23, 1 / 0.15, 1 / 1.9, backend = "spam")
cat(sprintf("sparse_P2D Matrix nnz=%d  spam nnz=%d\n",
            Matrix::nnzero(P_M), length(P_S@entries)))

for (be in c("dense", "Matrix", "spam")) {
  t1 <- system.time({
    m1 <- pois_SAP(
      y = dat$ym, x1 = dat$xxc[, 1], x2 = dat$xxc[, 2],
      efine = dat$ec, C = dat$C_m, ndx = c(20, 20),
      elements = TRUE, sparse.backend = be
    )
  })
  cat(sprintf(
    "backend=%-6s  AIC=%.4f  tau=(%.6f,%.6f)  time=%.2fs  Cclass=%s\n",
    be, m1$aic, m1$var.comp[1], m1$var.comp[2], t1["elapsed"],
    paste(class(m1$matlist$C), collapse = "/")
  ))
}

stopifnot(abs(m1$aic - 396.9185) < 1e-3)
cat("SPARSE_BACKEND_OK\n")

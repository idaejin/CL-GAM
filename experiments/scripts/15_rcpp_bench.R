#!/usr/bin/env Rscript
# Compare pure-R vs RcppArmadillo SOP kernels (same AIC / τ²).
root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))

dat <- clgam_load_madrid(FALSE)
cat("spclmm", as.character(packageVersion("spclmm")), "\n")
cat("Rcpp kernels loaded:", exists("sap_solve_schur_cpp", mode = "function"), "\n")

fit_pair <- function(use_rcpp) {
  options(spclmm.use_rcpp = use_rcpp)
  label <- if (use_rcpp) "Rcpp" else "pure-R"
  t0 <- system.time({
    m0 <- pois_SAP(
      y = dat$ym, x1 = dat$xxm[, 1], x2 = dat$xxm[, 2],
      efine = dat$em, C = diag(length(dat$ym)), ndx = c(20, 20), elements = TRUE
    )
  })
  t1 <- system.time({
    m1 <- pois_SAP(
      y = dat$ym, x1 = dat$xxc[, 1], x2 = dat$xxc[, 2],
      efine = dat$ec, C = dat$C_m, ndx = c(20, 20), elements = TRUE
    )
  })
  list(
    label = label,
    m0_aic = m0$aic, m1_aic = m1$aic,
    m0_tau = m0$var.comp, m1_tau = m1$var.comp,
    m0_sec = unname(t0["elapsed"]), m1_sec = unname(t1["elapsed"])
  )
}

# Pure R first (no compiled path), then Rcpp
r <- fit_pair(FALSE)
c <- fit_pair(TRUE)

cat(sprintf(
  "%s model0: AIC=%.4f  time=%.2fs\n%s model0: AIC=%.4f  time=%.2fs\n",
  r$label, r$m0_aic, r$m0_sec, c$label, c$m0_aic, c$m0_sec
))
cat(sprintf(
  "%s model1: AIC=%.4f  time=%.2fs\n%s model1: AIC=%.4f  time=%.2fs\n",
  r$label, r$m1_aic, r$m1_sec, c$label, c$m1_aic, c$m1_sec
))

stopifnot(abs(r$m0_aic - 374.5841) < 1e-3, abs(c$m0_aic - 374.5841) < 1e-3)
stopifnot(abs(r$m1_aic - 396.9185) < 1e-3, abs(c$m1_aic - 396.9185) < 1e-3)
stopifnot(max(abs(r$m0_tau - c$m0_tau)) < 1e-6)
stopifnot(max(abs(r$m1_tau - c$m1_tau)) < 1e-6)

out <- list(pure_R = r, Rcpp = c, version = as.character(packageVersion("spclmm")))
saveRDS(out, file.path(CLGAM_OUTPUT, "rcpp_benchmark.rds"))
cat(sprintf(
  "speedup model0=%.2fx  model1=%.2fx\n",
  r$m0_sec / c$m0_sec, r$m1_sec / c$m1_sec
))
cat("RCPP BENCH OK\n")

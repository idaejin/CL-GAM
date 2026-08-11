#!/usr/bin/env Rscript
# Benchmark + numerical check: sped-up spclmm vs paper reference numbers.
root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))

dat <- clgam_load_madrid(FALSE)
cat("spclmm", as.character(packageVersion("spclmm")), "\n")

t0 <- system.time({
  m0 <- pois_SAP(
    y = dat$ym, x1 = dat$xxm[, 1], x2 = dat$xxm[, 2],
    efine = dat$em, C = diag(length(dat$ym)), ndx = c(20, 20), elements = TRUE
  )
})
cat(sprintf("model0: AIC=%.4f (ref 374.5841)  tau=%s  time=%.2fs\n",
            m0$aic, paste(round(m0$var.comp, 6), collapse = ","), t0["elapsed"]))

t1 <- system.time({
  m1 <- pois_SAP(
    y = dat$ym, x1 = dat$xxc[, 1], x2 = dat$xxc[, 2],
    efine = dat$ec, C = dat$C_m, ndx = c(20, 20), elements = TRUE
  )
})
cat(sprintf("model1: AIC=%.4f (ref 396.9185)  tau=%s  time=%.2fs\n",
            m1$aic, paste(round(m1$var.comp, 6), collapse = ","), t1["elapsed"]))

stopifnot(abs(m0$aic - 374.5841) < 1e-3)
stopifnot(abs(m1$aic - 396.9185) < 1e-3)

saveRDS(
  list(model0_sec = unname(t0["elapsed"]), model1_sec = unname(t1["elapsed"]),
       m0_aic = m0$aic, m1_aic = m1$aic, version = as.character(packageVersion("spclmm"))),
  file.path(CLGAM_OUTPUT, "speed_benchmark.rds")
)
cat("SPEEDCHECK OK\n")

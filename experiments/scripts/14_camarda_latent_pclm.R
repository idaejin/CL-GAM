#!/usr/bin/env Rscript
# Camarda & Durbán (arXiv:2412.04956) latent-response PCLM vs pois_SAP (Madrid ATA).
# From experiments/: Rscript scripts/14_camarda_latent_pclm.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))

dat <- clgam_load_madrid(FALSE)
cat("spclmm", as.character(packageVersion("spclmm")), "\n")
stopifnot("pois_PCLM_latent" %in% getNamespaceExports("spclmm"))

message("=== pois_SAP (SOP / paper) ===")
t_sap <- system.time({
  sap <- pois_SAP(
    y = dat$ym, x1 = dat$xxc[, 1], x2 = dat$xxc[, 2],
    efine = dat$ec, C = dat$C_m, ndx = c(20, 20), elements = TRUE
  )
})

# Map SAP variance components to Camarda λ roughly as 1/τ² (order-of-magnitude)
lam <- 1 / pmax(as.numeric(sap$var.comp), 1e-8)
message(sprintf("SAP tau2 -> lambda trial: %s", paste(round(lam, 4), collapse = ", ")))

message("=== pois_PCLM_latent (Camarda–Durbán) ===")
t_cd <- system.time({
  cd <- pois_PCLM_latent(
    y = dat$ym, x1 = dat$xxc[, 1], x2 = dat$xxc[, 2],
    efine = dat$ec, C = dat$C_m, ndx = c(20, 20),
    lambda = lam, thr = 1e-6, maxit = 80L, trace = FALSE
  )
})

# Also try a coarser λ grid if correlation is poor
if (cor(sap$eta, cd$eta) < 0.9) {
  message("Low corr — retry with lambda = c(0.1, 0.1)")
  t_cd <- system.time({
    cd <- pois_PCLM_latent(
      y = dat$ym, x1 = dat$xxc[, 1], x2 = dat$xxc[, 2],
      efine = dat$ec, C = dat$C_m, ndx = c(20, 20),
      lambda = c(0.1, 0.1), thr = 1e-6, maxit = 80L
    )
  })
}

cat(sprintf(
  "SAP:  AIC=%.3f  time=%.2fs  ed=%.2f\n",
  sap$aic, t_sap["elapsed"], sap$ed
))
cat(sprintf(
  "CD:   dev=%.3f  time=%.2fs  niter=%d  lambda=(%.4g,%.4g)\n",
  cd$dev, t_cd["elapsed"], cd$niter, cd$lambda[1], cd$lambda[2]
))
cat(sprintf(
  "corr(eta)=%.4f  RMSE=%.4f  speedup=%.2fx\n",
  cor(sap$eta, cd$eta),
  sqrt(mean((sap$eta - cd$eta)^2)),
  as.numeric(t_sap["elapsed"] / t_cd["elapsed"])
))

saveRDS(
  list(
    sap_aic = sap$aic,
    sap_sec = unname(t_sap["elapsed"]),
    cd_dev = cd$dev,
    cd_sec = unname(t_cd["elapsed"]),
    cd_lambda = cd$lambda,
    corr_eta = cor(sap$eta, cd$eta),
    rmse_eta = sqrt(mean((sap$eta - cd$eta)^2)),
    ref = "Camarda & Durban arXiv:2412.04956"
  ),
  file.path(CLGAM_OUTPUT, "camarda_latent_pclm.rds")
)
message("Wrote ", file.path(CLGAM_OUTPUT, "camarda_latent_pclm.rds"))

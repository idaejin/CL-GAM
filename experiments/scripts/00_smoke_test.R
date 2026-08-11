#!/usr/bin/env Rscript
# Smoke test: paths, tidy data dims, maps, optional spclmm install.
# Run from experiments/:  Rscript scripts/00_smoke_test.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/02_load_madrid.R"))

dat <- clgam_load_madrid(load_maps = TRUE)

cat("Municipalities (y):", length(dat$ym), "\n")
cat("Census tracts (fine):", length(dat$yc), "\n")
cat("C_m dim:", paste(dim(dat$C_m), collapse = " x "), "\n")
cat("Covariates ncol:", ncol(dat$covariates), "\n")
cat("sum(C_m %*% yc) vs sum(ym):", sum(as.numeric(dat$C_m %*% dat$yc)), "vs", sum(dat$ym), "\n")
stopifnot(abs(sum(as.numeric(dat$C_m %*% dat$yc)) - sum(dat$ym)) < 1e-6)

cat("map mun polygons:", length(dat$map_mun$sp), "\n")
cat("map ct polygons:", length(dat$map_ct$sp), "\n")

# Package load / install (may need write access to R library)
tryCatch(
  {
    source(file.path(root, "R/01_load_spclmm.R"))
    cat("spclmm version:", as.character(packageVersion("spclmm")), "\n")
    cat("pois_SAP exported:", "pois_SAP" %in% getNamespaceExports("spclmm"), "\n")
  },
  error = function(e) {
    message("spclmm install/load skipped or failed: ", conditionMessage(e))
    message("Install manually: install.packages(CLGAM_SPCLMM_SRC, repos=NULL, type='source')")
  }
)

cat("SMOKE OK\n")

#!/usr/bin/env Rscript
# Adapt disaggregation (Nandi et al., JSS 2023) to Madrid MEDEA.
#
# Bridge: mun counts (ym) + CT covariates/expected rasterized to a grid.
# Not identical to CL-GAMM Case A (areal P-splines); fair modern competitor.
#
# From experiments/:
#   Sys.setenv(CLGAM_FAST = "1")   # coarse grid, fewer iters (smoke)
#   Rscript scripts/06_disaggregation_madrid.R
#
# Requires: disaggregation, terra, sf (install into R_libs if needed).

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/02_load_madrid.R"))
source(file.path(root, "R/04_disaggregation.R"))

# Prefer local experiment library (disaggregation installed here)
if (dir.exists(CLGAM_RLIBS)) {
  .libPaths(c(CLGAM_RLIBS, .libPaths()))
}

if (!requireNamespace("disaggregation", quietly = TRUE)) {
  stop(
    "Install disaggregation first, e.g.\n",
    "  install.packages('disaggregation', lib = '", CLGAM_RLIBS, "')"
  )
}
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Package 'terra' is required")
}

# Fine enough that almost all CTs intersect a cell; FAST only shrinks mesh/iters.
# 100 m ≈ 23 empty CT; 200 m ≈ 810 empty (then claimed to nearest NA cell).
res_m <- if (CLGAM_FAST) 200 else 100
iters <- if (CLGAM_FAST) 80 else 300
# iid polygon effect can be unstable on smoke meshes; keep on for full runs
use_iid <- !CLGAM_FAST

message("=== disaggregation Madrid bridge ===")
message("res_m=", res_m, "  iterations=", iters, "  FAST=", CLGAM_FAST)

dat <- clgam_load_madrid(load_maps = TRUE)

t0 <- system.time({
  out <- clgam_disag_fit_madrid(
    dat,
    cov_names = c("unemployed", "ageing"),
    res_m = res_m,
    iterations = iters,
    field = TRUE,
    iid = use_iid,
    fast = CLGAM_FAST,
    silent = TRUE
  )
})

message(sprintf("fit elapsed: %.1fs", unname(t0["elapsed"])))

# Predictions (pixel incidence / intensity × aggregation → counts)
message("Predicting pixel surface...")
preds <- tryCatch(
  disaggregation::predict_model(out$fit, newdata = NULL),
  error = function(e) {
    message("predict_model failed: ", conditionMessage(e))
    NULL
  }
)

# Aggregate predicted mean intensity*agg back to municipalities for mass check
mass <- NULL
if (!is.null(preds)) {
  # preds structure varies by version; try common slots
  mean_r <- NULL
  if (inherits(preds, "SpatRaster")) {
    mean_r <- preds
  } else if (is.list(preds)) {
    if (!is.null(preds$mean_predictions)) mean_r <- preds$mean_predictions
    else if (!is.null(preds$prediction)) mean_r <- preds$prediction
    else if (!is.null(preds$mean)) mean_r <- preds$mean
  }
  if (!is.null(mean_r) && inherits(mean_r, "SpatRaster")) {
    # For poisson/log, package predict usually returns rate or field; use
    # polygon prediction from fitted object if available
    mass <- list(note = "see fit$sdreport / plot; mass check via extract if rate*agg")
  }
}

# Quick covariate raster sanity + mun sf
summary_df <- data.frame(
  fast = CLGAM_FAST,
  res_m = res_m,
  iterations = iters,
  field = TRUE,
  iid = use_iid,
  n_mun = nrow(out$mun_sf),
  n_cov = length(out$cov_names),
  n_cells_cov = terra::ncell(out$rasters$covariates),
  n_ct_zero_cells = sum(out$rasters$n_cells_per_ct == 0),
  n_empty_painted = out$rasters$n_empty_painted,
  elapsed_sec = unname(t0["elapsed"]),
  stringsAsFactors = FALSE
)

# Parameter table if present
par_df <- tryCatch({
  s <- summary(out$fit$sd_out)
  as.data.frame(s)
}, error = function(e) NULL)

saveRDS(
  list(
    note = paste(
      "disaggregation on Madrid via CT→raster bridge.",
      "Covariates piecewise-constant on census tracts;",
      "aggregation_raster = expected per pixel within CT.",
      "Join mun response by GEOCODIGO (tidy order ≠ shapefile order)."
    ),
    summary = summary_df,
    parameters = par_df,
    fit = out$fit,
    data = out$data,
    mun_sf = out$mun_sf,
    cov_names = out$cov_names,
    res_m = res_m,
    preds = preds,
    elapsed = t0
  ),
  file.path(CLGAM_OUTPUT, "disaggregation_madrid_fit.rds")
)
utils::write.csv(
  summary_df,
  file.path(CLGAM_OUTPUT, "disaggregation_madrid_summary.csv"),
  row.names = FALSE
)

# Map: municipal observed SMR vs placeholder — save covariate rasters PDF
pdf(file.path(CLGAM_OUTPUT, "disaggregation_madrid_rasters.pdf"),
    width = 9, height = 4)
old <- par(mfrow = c(1, 3))
terra::plot(out$rasters$covariates[[1]], main = "unemployed (CT→raster)")
terra::plot(out$rasters$covariates[[2]], main = "ageing (CT→raster)")
terra::plot(out$rasters$aggregation, main = "expected per cell")
par(old)
dev.off()

message("Wrote ", file.path(CLGAM_OUTPUT, "disaggregation_madrid_fit.rds"))
message("Wrote ", file.path(CLGAM_OUTPUT, "disaggregation_madrid_summary.csv"))
message("DISAGG Madrid OK")

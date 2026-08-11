# ST-PCLM composition adapted to nested ATA (+ optional raster fine support).
#
# Lee et al. (2022) PLoS ONE: C_st = C_s ⊗ C_t with C_s typically mun→grid (ATP).
# Nested ATA: C_s = mun→CT (0/1), as in Madrid CL-GAMM / spclmm::pois_SAP.
# With raster: C_s = mun→grid = C_mun_ct %*% C_ct_grid (two-step nested).
#
# See experiments/improvements/ST_PCLM_ATA.md

#' Spatio-temporal composition C_st = kronecker(C_t, C_s) (column-major space-fast).
#' Matches Lee et al. (2022) Eq. (4) when vec is stacked by time blocks of spatial units.
clgam_C_st <- function(C_s, C_t) {
  C_s <- as.matrix(C_s)
  C_t <- as.matrix(C_t)
  kronecker(C_t, C_s)
}

#' Identity / aggregation time composition: each coarse interval owns `reps` fine steps.
#' Example: 12 months → 52/53 weeks needs a custom C_t; this helper is for equal splits.
clgam_C_t_block <- function(n_coarse, reps) {
  stopifnot(length(reps) == 1L || length(reps) == n_coarse)
  if (length(reps) == 1L) reps <- rep(as.integer(reps), n_coarse)
  m <- sum(reps)
  C <- matrix(0, n_coarse, m)
  j <- 1L
  for (i in seq_len(n_coarse)) {
    idx <- j:(j + reps[i] - 1L)
    C[i, idx] <- 1
    j <- j + reps[i]
  }
  C
}

#' CT → grid incidence from a terra raster of CT indices (cell value = CT row id).
#' Returns list(C_ct_grid, cell_id, ct_id, xy) for non-NA cells only.
clgam_C_ct_grid_from_idx <- function(idx_r) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required")
  }
  vals <- as.integer(terra::values(idx_r, mat = FALSE))
  cells <- which(!is.na(vals))
  ct_id <- vals[cells]
  n_ct <- max(ct_id, na.rm = TRUE)
  m <- length(cells)
  # sparse-friendly dense build (Madrid grids can be large — keep only active cells)
  C <- matrix(0, n_ct, m)
  C[cbind(ct_id, seq_len(m))] <- 1
  xy <- terra::xyFromCell(idx_r, cells)
  list(
    C_ct_grid = C,
    cell = cells,
    ct_id = ct_id,
    xy = xy,
    n_ct = n_ct,
    n_cell = m
  )
}

#' Mun → grid via nested product C_mun_ct %*% C_ct_grid.
clgam_C_mun_grid <- function(C_mun_ct, C_ct_grid) {
  C_mun_ct <- as.matrix(C_mun_ct)
  C_ct_grid <- as.matrix(C_ct_grid)
  stopifnot(ncol(C_mun_ct) == nrow(C_ct_grid))
  C_mun_ct %*% C_ct_grid
}

#' Build nested raster composition bundle from Madrid dat + CT rasters helper.
#'
#' Uses the same CT→raster painting as disaggregation (`clgam_ct_to_rasters`),
#' then forms:
#'   C_ct_grid, C_mun_grid = C_m %*% C_ct_grid
clgam_stpclm_ata_raster_C <- function(
  dat,
  cov_names = c("unemployed", "ageing"),
  res_m = 200
) {
  if (is.null(dat$map_ct) || is.null(dat$map_mun)) {
    stop("Need maps: clgam_load_madrid(load_maps = TRUE)")
  }
  if (!exists("clgam_ct_sf_attrs", mode = "function")) {
    stop("Source R/04_disaggregation.R first")
  }
  ct_sf <- clgam_ct_sf_attrs(dat, cov_names = cov_names)
  ras <- clgam_ct_to_rasters(ct_sf, cov_names = cov_names, res_m = res_m)
  if (is.null(ras$ct_index)) {
    stop("clgam_ct_to_rasters must return ct_index raster")
  }

  cg <- clgam_C_ct_grid_from_idx(ras$ct_index)
  C_mun_ct <- as.matrix(dat$C_m)
  stopifnot(ncol(C_mun_ct) == cg$n_ct)
  C_mun_grid <- clgam_C_mun_grid(C_mun_ct, cg$C_ct_grid)

  e_cell <- as.numeric(terra::values(ras$aggregation, mat = FALSE)[cg$cell])
  if (anyNA(e_cell)) {
    # fallback: split CT expected equally across its cells
    n_per <- as.numeric(table(factor(cg$ct_id, levels = seq_len(cg$n_ct))))
    e_cell <- as.numeric(dat$ec[cg$ct_id] / pmax(n_per[cg$ct_id], 1))
  }

  list(
    C_mun_ct = C_mun_ct,
    C_ct_grid = cg$C_ct_grid,
    C_mun_grid = C_mun_grid,
    xy_grid = cg$xy,
    cell = cg$cell,
    ct_id = cg$ct_id,
    e_cell = e_cell,
    e_ct = as.numeric(dat$ec),
    rasters = ras,
    res_m = res_m,
    note = paste(
      "ST-PCLM ATA nested + raster: fine support = grid cells;",
      "C_mun_grid = C_mun_ct %*% C_ct_grid (nested).",
      "Spatial ATA (CL-GAMM) uses C_mun_ct; ATP-via-nested uses C_mun_grid."
    )
  )
}

#' Aggregate a grid vector to CT via C_ct_grid (sum of cells in each CT).
clgam_grid_to_ct <- function(v_grid, C_ct_grid) {
  as.numeric(as.matrix(C_ct_grid) %*% as.numeric(v_grid))
}

#' Fit spatial PCLM on nested ATA (mun→CT) — time-degenerate ST-PCLM.
#' Thin wrapper around pois_SAP for provenance clarity.
clgam_stpclm_ata_fit <- function(
  y_mun,
  e_ct,
  C_mun_ct,
  lon_ct,
  lat_ct,
  nlcovfine = NULL,
  lcovfine = NULL,
  ndx = c(20L, 20L),
  ndxnl = NULL,
  trace = TRUE
) {
  if (!exists("pois_SAP", mode = "function")) {
    stop("pois_SAP not found; source R/01_load_spclmm.R")
  }
  args <- list(
    y = y_mun,
    x1 = lon_ct,
    x2 = lat_ct,
    efine = e_ct,
    C = C_mun_ct,
    ndx = ndx,
    elements = TRUE,
    trace = trace
  )
  if (!is.null(nlcovfine)) {
    args$nlcovfine <- nlcovfine
    if (is.null(ndxnl)) {
      args$ndxnl <- rep(12L, ncol(as.matrix(nlcovfine)))
    } else {
      args$ndxnl <- ndxnl
    }
  }
  if (!is.null(lcovfine)) args$lcovfine <- lcovfine
  fit <- do.call(pois_SAP, args)
  fit$composition <- "ATA_nested_mun_ct"
  fit$st_pclm_note <- paste(
    "Spatial section of ST-PCLM (Lee et al. 2022) with C_s = nested mun→CT.",
    "No C_t (Madrid CVD tidy has no time). Estimator: pois_SAP (SOP)."
  )
  fit
}

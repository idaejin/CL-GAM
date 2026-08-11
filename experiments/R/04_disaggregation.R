# Bridge Madrid MEDEA (areal CT covariates) → disaggregation::prepare_data inputs.
#
# disaggregation wants: sf polygons (coarse response) + SpatRaster covariates
# (+ optional aggregation raster). We rasterize CT polygons so that cov/expected
# are piecewise-constant on tracts — an honest but approximate competitor to
# Case A CL-GAMM (which keeps areal support).

#' Municipality code (3-digit) for each tidy mun row, via composition matrix.
clgam_mun_codes_from_C <- function(C, ct_geocodigo) {
  C <- as.matrix(C)
  ct_geocodigo <- as.character(ct_geocodigo)
  stopifnot(ncol(C) == length(ct_geocodigo))
  vapply(seq_len(nrow(C)), function(i) {
    j <- which(C[i, ] > 0)[1]
    substr(ct_geocodigo[j], 1, 3)
  }, character(1))
}

#' Build municipality sf with response aligned to tidy `ym` (join by GEOCODIGO).
clgam_mun_sf_for_disag <- function(dat, id_var = "area_id", response_var = "response") {
  if (is.null(dat$map_mun)) {
    stop("dat$map_mun missing; call clgam_load_madrid(load_maps = TRUE)")
  }
  mun_codes <- clgam_mun_codes_from_C(dat$C_m, dat$map_ct$sf$GEOCODIGO)
  stopifnot(length(unique(mun_codes)) == length(mun_codes))

  shp <- dat$map_mun$sf
  shp$GEOCODIGO <- as.character(shp$GEOCODIGO)
  idx <- match(mun_codes, shp$GEOCODIGO)
  if (anyNA(idx)) {
    stop("Could not match mun codes to map_mun$GEOCODIGO: ",
         paste(mun_codes[is.na(idx)], collapse = ", "))
  }

  out <- shp[idx, ]
  out[[id_var]] <- mun_codes
  out[[response_var]] <- as.numeric(dat$ym)
  out$expected_mun <- as.numeric(dat$em)
  # tiny buffer fixes occasional invalid geometries for terra extract
  out <- sf::st_buffer(out, dist = 0)
  out
}

#' Census-tract sf with covariates + expected (same row order as tidy CT).
clgam_ct_sf_attrs <- function(dat, cov_names = c("unemployed", "ageing")) {
  if (is.null(dat$map_ct)) {
    stop("dat$map_ct missing; call clgam_load_madrid(load_maps = TRUE)")
  }
  ct <- dat$map_ct$sf
  stopifnot(nrow(ct) == length(dat$yc), nrow(ct) == nrow(dat$covariates))

  cov <- as.data.frame(dat$covariates)[, cov_names, drop = FALSE]
  for (nm in cov_names) ct[[nm]] <- cov[[nm]]
  ct$expected <- as.numeric(dat$ec)
  ct$GEOCODIGO <- as.character(ct$GEOCODIGO)
  ct <- sf::st_buffer(ct, dist = 0)
  ct
}

#' Raster template covering Madrid maps (resolution in CRS units, metres).
clgam_madrid_raster_template <- function(sf_obj, res_m) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required")
  }
  bb <- sf::st_bbox(sf_obj)
  terra::rast(
    xmin = bb["xmin"], xmax = bb["xmax"],
    ymin = bb["ymin"], ymax = bb["ymax"],
    resolution = res_m,
    crs = sf::st_crs(sf_obj)$wkt
  )
}

#' Build covariate + aggregation SpatRasters from CT polygons.
#'
#' Primary: `terra::rasterize(..., touches = TRUE)` so each cell inherits its
#' intersecting census tract. Tiny tracts with no cell get the nearest empty
#' cell at their centroid (last-resort paint).
#'
#' Aggregation layer = expected count **per pixel** within each CT
#' (`ec / n_pixels_in_ct`), so municipal sums ≈ `em`.
clgam_ct_to_rasters <- function(
  ct_sf,
  cov_names = c("unemployed", "ageing"),
  res_m = 100,
  template = NULL
) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required")
  }
  if (is.null(template)) {
    template <- clgam_madrid_raster_template(ct_sf, res_m)
  }

  # touches=TRUE so small CTs still claim at least one cell when possible
  vct <- terra::vect(ct_sf)
  cov_list <- lapply(cov_names, function(nm) {
    r <- terra::rasterize(vct, template, field = nm, touches = TRUE)
    names(r) <- nm
    r
  })
  cov_stack <- terra::rast(cov_list)

  ct_sf$.ct_i <- seq_len(nrow(ct_sf))
  vct <- terra::vect(ct_sf)
  idx_r <- terra::rasterize(vct, template, field = ".ct_i", touches = TRUE)
  cell_ct <- terra::freq(idx_r)
  n_cells <- integer(nrow(ct_sf))
  if (!is.null(cell_ct) && nrow(cell_ct)) {
    ok <- !is.na(cell_ct$value)
    n_cells[cell_ct$value[ok]] <- as.integer(cell_ct$count[ok])
  }

  # CTs still empty: claim nearest currently-NA cell to centroid
  miss <- which(n_cells == 0L)
  if (length(miss)) {
    message("Raster res leaves ", length(miss),
            " CTs with 0 cells; claiming nearest empty cell.")
    cents <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(ct_sf[miss, ])))
    na_cells <- which(is.na(terra::values(idx_r, mat = FALSE)))
    if (!length(na_cells)) {
      stop("No empty raster cells left to assign tiny census tracts")
    }
    xy_na <- terra::xyFromCell(idx_r, na_cells)
    for (k in seq_along(miss)) {
      d2 <- (xy_na[, 1] - cents[k, 1])^2 + (xy_na[, 2] - cents[k, 2])^2
      j <- which.min(d2)
      cell <- na_cells[j]
      idx_r[cell] <- miss[k]
      n_cells[miss[k]] <- 1L
      # remove claimed cell from pool
      na_cells <- na_cells[-j]
      xy_na <- xy_na[-j, , drop = FALSE]
      if (!length(na_cells) && k < length(miss)) {
        stop("Ran out of empty cells while painting tiny CTs")
      }
    }
    # refresh covariate values on newly claimed cells from CT attrs
    for (nm in cov_names) {
      vals <- terra::values(cov_stack[[nm]], mat = FALSE)
      for (k in seq_along(miss)) {
        # find cell now holding this CT
        cells_k <- which(terra::values(idx_r, mat = FALSE) == miss[k])
        vals[cells_k] <- ct_sf[[nm]][miss[k]]
      }
      terra::values(cov_stack[[nm]]) <- vals
    }
  }

  e_per <- ct_sf$expected / pmax(n_cells, 1L)
  # build aggregation from idx_r (avoids second polygon rasterize fights)
  agg <- terra::rast(idx_r)
  terra::values(agg) <- NA_real_
  idx_vals <- terra::values(idx_r, mat = FALSE)
  agg_vals <- rep(NA_real_, length(idx_vals))
  ok <- !is.na(idx_vals)
  agg_vals[ok] <- e_per[idx_vals[ok]]
  terra::values(agg) <- agg_vals
  names(agg) <- "expected_per_cell"

  list(
    covariates = cov_stack,
    aggregation = agg,
    ct_index = idx_r,
    n_cells_per_ct = n_cells,
    template = template,
    res_m = res_m,
    n_empty_painted = length(miss)
  )
}

#' Default mesh args in projected metres (Madrid UTM).
clgam_disag_mesh_args <- function(fast = FALSE) {
  if (fast) {
    list(
      max.edge = c(5000, 20000),
      offset = c(5000, 20000),
      cutoff = 2000
    )
  } else {
    list(
      max.edge = c(2000, 8000),
      offset = c(3000, 15000),
      cutoff = 800
    )
  }
}

#' Fit disaggregation Poisson model on Madrid bridge objects.
clgam_disag_fit_madrid <- function(
  dat,
  cov_names = c("unemployed", "ageing"),
  res_m = 1000,
  iterations = 200,
  field = TRUE,
  iid = TRUE,
  fast = FALSE,
  silent = TRUE
) {
  if (!requireNamespace("disaggregation", quietly = TRUE)) {
    stop(
      "Package 'disaggregation' not installed. ",
      "install.packages('disaggregation') or use experiments/R_libs."
    )
  }

  mun_sf <- clgam_mun_sf_for_disag(dat)
  ct_sf <- clgam_ct_sf_attrs(dat, cov_names = cov_names)
  ras <- clgam_ct_to_rasters(ct_sf, cov_names = cov_names, res_m = res_m)

  # Scale covariates like typical disaggregation workflows
  cov_scaled <- terra::scale(ras$covariates)

  data_for_model <- disaggregation::prepare_data(
    polygon_shapefile = mun_sf,
    covariate_rasters = cov_scaled,
    aggregation_raster = ras$aggregation,
    id_var = "area_id",
    response_var = "response",
    mesh_args = clgam_disag_mesh_args(fast = fast),
    na_action = TRUE
  )

  # disag_model checks inherits(iterations, "numeric"); integer fails that test
  fit <- disaggregation::disag_model(
    data_for_model,
    iterations = as.numeric(iterations),
    family = "poisson",
    link = "log",
    field = field,
    iid = iid,
    silent = silent
  )

  list(
    fit = fit,
    data = data_for_model,
    mun_sf = mun_sf,
    rasters = ras,
    cov_names = cov_names,
    res_m = res_m,
    iterations = iterations
  )
}

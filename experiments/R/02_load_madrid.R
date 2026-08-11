# Load Madrid MEDEA tidy counts, covariates, composition matrices, and maps.
# Prefer sf::st_read; convert to sp::SpatialPolygonsDataFrame for spclmm::choromap.

.clgam_source_paths <- function() {
  here <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(here) && nzchar(here)) {
    source(file.path(dirname(here), "00_paths.R"), local = FALSE)
  } else if (file.exists("R/00_paths.R")) {
    source("R/00_paths.R")
  } else if (file.exists("../R/00_paths.R")) {
    source("../R/00_paths.R")
  } else {
    stop("Cannot find R/00_paths.R; source it first or setwd to experiments/")
  }
}

if (!exists("CLGAM_TIDY")) .clgam_source_paths()

.clgam_read_shape <- function(shp_path, id_var = "GEOCODIGO") {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required to read cartography. install.packages('sf')")
  }
  sf_obj <- sf::st_read(shp_path, quiet = TRUE)
  if (!id_var %in% names(sf_obj)) {
    warning("ID column ", id_var, " not found in ", basename(shp_path))
  }
  # choromap / legacy plot expect sp
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package 'sp' is required for Spatial* conversion used by choromap.")
  }
  sp_obj <- as(sf_obj, "Spatial")
  list(sf = sf_obj, sp = sp_obj)
}

#' Load all objects needed for paper Codes2-style fits
#' @return named list
clgam_load_madrid <- function(load_maps = TRUE) {
  mun_f <- utils::read.table(
    file.path(CLGAM_TIDY, "femaleGG7_mun.txt"),
    header = TRUE
  )
  mun_m <- utils::read.table(
    file.path(CLGAM_TIDY, "maleGG7_mun.txt"),
    header = TRUE
  )
  ct_f <- utils::read.table(
    file.path(CLGAM_TIDY, "femaleGG7_ct.txt"),
    header = TRUE
  )
  ct_m <- utils::read.table(
    file.path(CLGAM_TIDY, "maleGG7_ct.txt"),
    header = TRUE
  )

  ym <- mun_f$observed + mun_m$observed
  em <- mun_f$expected + mun_m$expected
  yc <- ct_f$observed + ct_m$observed
  ec <- ct_f$expected + ct_m$expected

  xxm <- cbind(mun_f$lon, mun_f$lat)
  xxc <- cbind(ct_f$lon, ct_f$lat)

  # Composition municipality → census tract (179 × 3906)
  env_C <- new.env(parent = emptyenv())
  load(file.path(CLGAM_TIDY, "spC_mun-ct.Rdata"), envir = env_C)
  C_m <- env_C$C.m
  stopifnot(is.matrix(C_m) || inherits(C_m, "Matrix"))
  stopifnot(nrow(C_m) == length(ym), ncol(C_m) == length(yc))

  # Covariates at CT (paper Codes2)
  evc <- utils::read.csv2(
    file.path(CLGAM_TIDY, "covariates2001.csv"),
    header = TRUE,
    sep = ","
  )
  # columns as in PAPER - Codes2.R
  expl <- cbind(
    as.numeric(as.character(evc$permanuales)),
    as.numeric(as.character(evc$perdesempleados)),
    as.numeric(as.character(evc$perasalariados)),
    as.numeric(as.character(evc$perpolucion)),
    as.numeric(as.character(evc$perenvejecimiento)),
    as.numeric(as.character(evc$perinstrinsuf))
  )
  colnames(expl) <- c(
    "manual", "unemployed", "salaried",
    "pollution", "ageing", "low_edu"
  )
  expl_cen <- scale(expl, center = TRUE, scale = FALSE)

  out <- list(
    mun_female = mun_f,
    mun_male = mun_m,
    ct_female = ct_f,
    ct_male = ct_m,
    ym = ym,
    em = em,
    yc = yc,
    ec = ec,
    xxm = xxm,
    xxc = xxc,
    C_m = C_m,
    covariates = expl,
    covariates_cen = expl_cen,
    logSMR_mun = log((ym + 1) / em)
  )

  if (load_maps) {
    mun_shp <- file.path(CLGAM_CARTO, "Municipios 2001", "200001099.shp")
    ct_shp <- file.path(CLGAM_CARTO, "Seccionados 2001", "200001108.shp")
    out$map_mun <- .clgam_read_shape(mun_shp)
    out$map_ct <- .clgam_read_shape(ct_shp)
  }

  out
}

# Pennsylvania lung cancer (SpatialEpi::pennLC) nested ATA scaffold.
#
# Fine units = 67 counties. Coarse units = k-means regions of county centroids
# (administrative multi-county regions are not in the package). Smoking is
# county-level (Case A); a coarse-mean covariate is expanded for Case B/C.

clgam_ensure_penn_packages <- function() {
  lib <- file.path(CLGAM_EXP, "R_libs")
  if (!dir.exists(lib)) dir.create(lib, recursive = TRUE)
  .libPaths(c(lib, .libPaths()))
  if (!requireNamespace("SpatialEpi", quietly = TRUE)) {
    stop("Need SpatialEpi in experiments/R_libs")
  }
  invisible(TRUE)
}

#' Load pennLC nested ATA bundle
#'
#' @param n_coarse Number of k-means coarse regions (default 20).
#' @param seed RNG seed for k-means.
clgam_load_pennLC <- function(n_coarse = 20L, seed = 1L) {
  clgam_ensure_penn_packages()
  suppressPackageStartupMessages(library(SpatialEpi))

  data(pennLC, package = "SpatialEpi", envir = environment())

  # Strata order required by expected()
  d <- pennLC$data
  d <- d[order(d$county, d$race, d$gender, d$age), ]
  e_strata <- expected(
    population = d$population,
    cases = d$cases,
    n.strata = 16L
  )
  # County totals in factor-level order (do not as.numeric() before name-index)
  y_tab <- tapply(d$cases, d$county, sum)
  pop_tab <- tapply(d$population, d$county, sum)
  counties <- names(y_tab)
  stopifnot(length(e_strata) == length(counties))
  y_county <- as.numeric(y_tab)
  pop_county <- as.numeric(pop_tab)
  e_county <- as.numeric(e_strata)
  names(y_county) <- names(pop_county) <- names(e_county) <- counties
  stopifnot(!anyNA(y_county), !anyNA(e_county), all(e_county > 0))

  geo <- pennLC$geo
  geo$county <- as.character(geo$county)
  geo <- geo[match(counties, geo$county), , drop = FALSE]
  stopifnot(identical(geo$county, counties))

  sm <- pennLC$smoking
  sm$county <- as.character(sm$county)
  sm <- sm[match(counties, sm$county), , drop = FALSE]
  smoking <- as.numeric(sm$smoking)
  stopifnot(!anyNA(smoking))

  x1 <- as.numeric(geo$x)
  x2 <- as.numeric(geo$y)
  n_fine <- length(counties)
  n_coarse <- as.integer(n_coarse)
  stopifnot(n_coarse >= 2L, n_coarse < n_fine)

  set.seed(as.integer(seed))
  km <- stats::kmeans(cbind(x1, x2), centers = n_coarse, nstart = 25L)
  group <- factor(km$cluster)
  n_group <- nlevels(group)

  C <- matrix(0, nrow = n_group, ncol = n_fine)
  for (j in seq_len(n_fine)) C[as.integer(group[j]), j] <- 1
  rownames(C) <- levels(group)
  colnames(C) <- counties

  y <- as.numeric(tapply(y_county, group, sum))
  names(y) <- levels(group)
  e_g <- as.numeric(tapply(e_county, group, sum))
  names(e_g) <- levels(group)

  # Recompute fine expected so C %*% e_fine = e_g using overall rate from coarse
  # Prefer stratified E from SpatialEpi; check mass
  # Use county E directly (age-sex-race standardised); y is aggregated from same data
  rate_check <- max(abs(as.numeric(C %*% e_county) - e_g))

  cov_fine <- cbind(
    smoking = as.numeric(scale(smoking, scale = FALSE))
  )
  rownames(cov_fine) <- counties

  # Aggregated covariate: region-mean latitude (climate / northness proxy), expanded
  lat_by_g <- as.numeric(tapply(x2, group, mean))
  names(lat_by_g) <- levels(group)
  lat_exp <- lat_by_g[as.integer(group)]
  cov_agg_exp <- cbind(
    north_agg = as.numeric(scale(lat_exp, scale = FALSE))
  )
  rownames(cov_agg_exp) <- counties

  list(
    y = y,
    e_fine = e_county,
    C = C,
    x1 = x1,
    x2 = x2,
    county = counties,
    group = as.character(group),
    group_levels = levels(group),
    y_fine = y_county,
    pop_fine = pop_county,
    smoking = smoking,
    covariates_fine = cov_fine,
    covariates_agg_expanded = cov_agg_exp,
    n_coarse = n_group,
    n_fine = n_fine,
    coarse_scheme = sprintf("kmeans_%d_seed%d", n_coarse, as.integer(seed)),
    e_mass_check = rate_check,
    source = "SpatialEpi::pennLC (PA lung cancer 2002 + smoking)",
    note = paste(
      "Fine = 67 counties; coarse = k-means regions.",
      "Case A: smoking at county; Case B: region-mean latitude expanded;",
      "Case C: both. Expected counts via SpatialEpi::expected (16 strata)."
    )
  )
}

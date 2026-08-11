# NY leukemia county→tract nested ATA (SpatialEpi::NYleukemia + DClusterm::NY8).
#
# Builds composition C (8 counties × 281 tracts) and fine/aggregated covariates
# for CL-GAMM Case A / B / C illustrations with pois_SAP.
#
# Nested equivalence (Factual): for partition C, a county-level additive effect
# h(x_a) outside the link equals expanding x_a as constant within each county and
# placing it in nlcovfine / lcovfine inside η before C. See manuscript Case B/C.

clgam_ensure_ny_packages <- function() {
  lib <- file.path(CLGAM_EXP, "R_libs")
  if (!dir.exists(lib)) dir.create(lib, recursive = TRUE)
  .libPaths(c(lib, .libPaths()))
  need <- c("SpatialEpi", "DClusterm", "sp")
  miss <- need[!vapply(need, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss)) {
    stop(
      "Missing packages for NY example: ", paste(miss, collapse = ", "),
      ". Install into experiments/R_libs (repos = cloud.r-project.org)."
    )
  }
  invisible(TRUE)
}

#' Load NY leukemia nested ATA bundle
#'
#' @param n_coarse Number of coarse units. \code{NULL} or \code{"county"} uses the
#'   8 true counties (often too few for P-spline ATA). An integer (default 40)
#'   builds nested aggregates by k-means on tract coordinates (reproducible seed).
#' @param seed RNG seed for k-means coarse units.
#' @return list with y, e_tract, C, coords, covariates, ...
clgam_load_ny_leukemia <- function(n_coarse = 40L, seed = 1L) {
  clgam_ensure_ny_packages()
  suppressPackageStartupMessages({
    library(SpatialEpi)
    library(DClusterm)
    library(sp)
  })

  data(NYleukemia, package = "SpatialEpi", envir = environment())
  data(NY8, package = "DClusterm", envir = environment())

  ny8 <- as.data.frame(NY8)
  ny8$AREAKEY <- as.character(ny8$AREAKEY)

  epi <- NYleukemia$data
  epi$censustract.FIPS <- as.character(epi$censustract.FIPS)
  geo <- NYleukemia$geo
  geo$censustract.FIPS <- as.character(geo$censustract.FIPS)

  # Align on FIPS / AREAKEY
  stopifnot(all(epi$censustract.FIPS %in% ny8$AREAKEY))
  ord <- match(epi$censustract.FIPS, ny8$AREAKEY)
  ny8 <- ny8[ord, , drop = FALSE]
  stopifnot(identical(epi$censustract.FIPS, ny8$AREAKEY))

  geo <- geo[match(epi$censustract.FIPS, geo$censustract.FIPS), , drop = FALSE]

  fips <- epi$censustract.FIPS
  n_tract <- length(fips)
  x1 <- as.numeric(geo$x)
  x2 <- as.numeric(geo$y)

  cases_tr <- as.numeric(ny8$Cases)
  if (anyNA(cases_tr)) cases_tr <- as.numeric(epi$cases)
  pop <- as.numeric(ny8$POP8)
  if (anyNA(pop)) pop <- as.numeric(epi$population)
  stopifnot(!anyNA(pop), all(pop > 0), !anyNA(cases_tr))

  county <- substr(fips, 1L, 5L)

  if (is.null(n_coarse) || identical(n_coarse, "county")) {
    group <- factor(county)
    coarse_scheme <- "county_FIPS"
  } else {
    n_coarse <- as.integer(n_coarse)
    stopifnot(n_coarse >= 2L, n_coarse < n_tract)
    set.seed(as.integer(seed))
    km <- stats::kmeans(cbind(x1, x2), centers = n_coarse, nstart = 20L)
    group <- factor(km$cluster)
    coarse_scheme <- sprintf("kmeans_%d_seed%d", n_coarse, as.integer(seed))
  }

  n_group <- nlevels(group)
  C <- matrix(0, nrow = n_group, ncol = n_tract)
  for (j in seq_len(n_tract)) {
    C[as.integer(group[j]), j] <- 1
  }
  rownames(C) <- levels(group)
  colnames(C) <- fips

  y <- as.numeric(tapply(cases_tr, group, function(z) round(sum(z))))
  names(y) <- levels(group)
  pop_g <- as.numeric(tapply(pop, group, sum))
  names(pop_g) <- levels(group)
  stopifnot(!anyNA(y), !anyNA(pop_g), sum(pop_g) > 0)

  rate <- sum(y) / sum(pop_g)
  e_tract <- as.numeric(pop * rate)
  stopifnot(!anyNA(e_tract), all(e_tract > 0))

  age65 <- as.numeric(ny8$PCTAGE65P)
  ownhome <- as.numeric(ny8$PCTOWNHOME)
  pexp <- as.numeric(ny8$PEXPOSURE)
  stopifnot(!anyNA(age65), !anyNA(ownhome), !anyNA(pexp))

  cov_fine <- cbind(
    ageing = as.numeric(scale(age65, scale = FALSE)),
    ownhome = as.numeric(scale(ownhome, scale = FALSE))
  )
  rownames(cov_fine) <- fips

  # Aggregated TCE: mean within each coarse unit, expanded to tracts
  pexp_by_g <- as.numeric(tapply(pexp, group, mean))
  names(pexp_by_g) <- levels(group)
  pexp_exp <- pexp_by_g[as.integer(group)]
  cov_agg_exp <- cbind(
    tce_agg = as.numeric(scale(pexp_exp, scale = FALSE))
  )
  rownames(cov_agg_exp) <- fips

  cov_agg <- data.frame(
    group = levels(group),
    y = y,
    pop = pop_g,
    tce = pexp_by_g,
    ageing = as.numeric(tapply(age65, group, mean)),
    stringsAsFactors = FALSE
  )

  list(
    y = y,
    e_tract = e_tract,
    C = C,
    x1 = x1,
    x2 = x2,
    fips = fips,
    group = as.character(group),
    group_levels = levels(group),
    county = county,
    pop_tract = pop,
    cases_tract = cases_tr,
    covariates_fine = cov_fine,
    covariates_agg_expanded = cov_agg_exp,
    covariates_agg = cov_agg,
    n_coarse = n_group,
    n_tract = n_tract,
    coarse_scheme = coarse_scheme,
    source = c(
      "SpatialEpi::NYleukemia",
      "DClusterm::NY8 (PCTAGE65P, PCTOWNHOME, PEXPOSURE)"
    ),
    note = paste(
      "Case B/C aggregated effects enter pois_SAP via coarse→tract expansion",
      "(constant within coarse unit); equivalent to h(x_a) for nested partition C.",
      "Default coarse units: k-means (8 true counties are too few for P-spline ATA)."
    )
  )
}

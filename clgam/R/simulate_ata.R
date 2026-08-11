#' Unit-square envelope as an \pkg{sf} polygon
#' @keywords internal
.clgam_unit_square <- function() {
  sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1, ymax = 1)))
}

#' Voronoi partition of \code{envelope} into exactly \code{n} polygons
#'
#' Retries until the tessellation has \code{n} cells that cover \code{envelope}
#' (no dropped slivers). Returns an \code{sfc} of length \code{n}.
#' @keywords internal
.clgam_voronoi <- function(n, envelope, min_area = 1e-8, max_try = 40L) {
  n <- as.integer(n)
  stopifnot(n >= 1L)
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))

  env <- sf::st_make_valid(envelope)
  env_area <- as.numeric(sf::st_area(env))
  box <- sf::st_as_sfc(sf::st_bbox(env))

  for (it in seq_len(max_try)) {
    seeds <- sf::st_sample(env, size = n, type = "random")
    tries <- 0L
    while (length(seeds) < n && tries < 25L) {
      seeds <- c(seeds, sf::st_sample(env, size = n - length(seeds), type = "random"))
      tries <- tries + 1L
    }
    if (length(seeds) < n) next
    seeds <- seeds[seq_len(n)]

    vor <- sf::st_voronoi(sf::st_combine(seeds), envelope = box)
    polys <- sf::st_collection_extract(sf::st_sfc(vor), "POLYGON")
    polys <- sf::st_intersection(polys, env)
    polys <- polys[sf::st_geometry_type(polys) %in% c("POLYGON", "MULTIPOLYGON")]
    # keep multipolygons as single cells (do NOT explode) so count stays n
    polys <- sf::st_cast(polys, "MULTIPOLYGON", warn = FALSE)
    areas <- as.numeric(sf::st_area(polys))
    polys <- polys[is.finite(areas) & areas > min_area]
    if (length(polys) != n) next

    uni <- tryCatch(sf::st_make_valid(sf::st_union(polys)), error = function(e) NULL)
    if (is.null(uni) || length(uni) == 0L) next
    # empty sym_difference => perfect cover (area length 0)
    sdif <- suppressWarnings(sf::st_sym_difference(uni, env))
    gap_area <- if (is.null(sdif) || length(sdif) == 0L) {
      0
    } else {
      aa <- as.numeric(sf::st_area(sf::st_make_valid(sdif)))
      if (length(aa) == 0L || !is.finite(aa[1L])) 0 else aa[1L]
    }
    if ((gap_area / env_area) < 1e-3) {
      return(sf::st_make_valid(polys))
    }
  }
  stop("Failed to build a ", n, "-cell Voronoi partition of the envelope.",
       call. = FALSE)
}

#' Simulate ATA counts on nested areal polygons (Voronoi)
#'
#' Generative story (fine \eqn{\rightarrow} coarse):
#' \enumerate{
#'   \item Build a \strong{smooth} latent field \(\eta(s)\) evaluated at fine
#'     centroids (piecewise constant on fine polygons): spatial surface plus
#'     optional nonlinear \eqn{g(z)}.
#'   \item Fine intensity \(\gamma = e_f \odot \exp(\eta)\), with exposures
#'     \eqn{e_f} proportional to fine area.
#'   \item Simulate fine counts \(y_f \sim \mathrm{Poisson}(\gamma)\).
#'   \item Aggregate the \emph{raw} coarse map \eqn{y = C y_f}.
#' }
#' For independent Poisson fine counts, \(y\sim\mathrm{Poisson}(C\gamma)\)
#' marginally, matching the composite-link mean in \code{clgam()}.
#'
#' Geometry: coarse Voronoi of the unit square; each coarse cell is itself
#' partitioned into exactly \code{n_fine_per} fine Voronoi subpolygons
#' (a nested partition; \code{C} is 0-1 with one 1 per column).
#'
#' Requires \pkg{sf}. No MEDEA / proprietary data.
#'
#' @param n_coarse number of coarse polygons
#' @param n_fine_per number of fine subpolygons \strong{per} coarse polygon
#' @param seed RNG seed
#' @param include_covariate logical; include a nonlinear covariate effect
#' @param covariate_level \code{"fine"} (Case A: \eqn{g(z_f)}), \code{"coarse"}
#'   (Case B: \eqn{h(z_a)} constant within each coarse unit), or \code{"both"}
#'   (Case C: fine \eqn{g} and coarse \eqn{h}). Cases B/C use fine expansion of
#'   \eqn{z_a} so \(\mu=\exp(\log(C\gamma)+h(z_a))\) when \eqn{h} is within-unit
#'   constant.
#' @param spatial_amp amplitude of the spatial field on fine centroids
#' @param nl_amp amplitude of nonlinear effect(s); length 1 or 2 for Case C
#'   as \code{c(g, h)}
#' @param nl_fun function of \(z\in[0,1]\), or length-2 \code{list} for Case C.
#'   Default is a centred quadratic \eqn{(z-1/2)^2}.
#'
#' @return A list with \code{y}, \code{y_fine}, \code{C}, \code{x1}/\code{x2},
#'   \code{efine}, \code{gamma}, \code{sf_coarse}/\code{sf_fine}, and truth.
#' @export
simulate_ata <- function(n_coarse = 12L,
                         n_fine_per = 8L,
                         seed = 1L,
                         include_covariate = FALSE,
                         covariate_level = c("fine", "coarse", "both"),
                         spatial_amp = 0.45,
                         nl_amp = 0.85,
                         nl_fun = function(z) (z - 0.5)^2) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("simulate_ata() requires the 'sf' package. install.packages(\"sf\")",
         call. = FALSE)
  }
  covariate_level <- match.arg(covariate_level)
  set.seed(as.integer(seed))
  n_coarse <- as.integer(n_coarse)
  n_fine_per <- as.integer(n_fine_per)
  stopifnot(n_coarse >= 2L, n_fine_per >= 2L)

  envelope <- .clgam_unit_square()
  coarse_geom <- .clgam_voronoi(n_coarse, envelope)
  stopifnot(length(coarse_geom) == n_coarse)

  fine_list <- vector("list", n_coarse)
  for (i in seq_len(n_coarse)) {
    env_i <- coarse_geom[i]
    fine_list[[i]] <- .clgam_voronoi(n_fine_per, env_i, min_area = 1e-10)
  }

  fine_geom <- do.call(c, unname(fine_list))
  n_fine <- length(fine_geom)
  stopifnot(n_fine == n_coarse * n_fine_per)
  coarse_id <- rep(seq_len(n_coarse), each = n_fine_per)

  C <- matrix(0, n_coarse, n_fine)
  C[cbind(coarse_id, seq_len(n_fine))] <- 1
  stopifnot(all(colSums(C) == 1), all(rowSums(C) == n_fine_per))

  cents <- sf::st_coordinates(sf::st_centroid(fine_geom))
  x1 <- as.numeric(cents[, 1])
  x2 <- as.numeric(cents[, 2])

  eta_spatial <- spatial_amp * sin(2 * pi * x1) * cos(2 * pi * x2) -
    0.15 * (x1 - 0.5)
  eta <- eta_spatial
  nlcovfine <- NULL
  z <- NULL
  z_f <- NULL
  z_a <- NULL
  g_true <- NULL
  h_true <- NULL
  h_a <- NULL
  case <- "spatial"

  .as_nl_funs <- function(nl_fun, n) {
    if (is.list(nl_fun)) {
      stopifnot(length(nl_fun) >= n)
      return(nl_fun[seq_len(n)])
    }
    rep(list(nl_fun), n)
  }
  .as_amps <- function(nl_amp, n) {
    a <- rep_len(as.numeric(nl_amp), n)
    stopifnot(length(a) == n, all(is.finite(a)))
    a
  }

  if (isTRUE(include_covariate)) {
    if (identical(covariate_level, "fine")) {
      case <- "A"
      funs <- .as_nl_funs(nl_fun, 1L)
      amps <- .as_amps(nl_amp, 1L)
      z_f <- stats::runif(n_fine, 0, 1)
      z <- z_f
      g_true <- as.numeric(amps[1] * funs[[1]](z_f))
      g_true <- g_true - mean(g_true)
      nlcovfine <- cbind(z_f = z_f)
      eta <- eta + g_true
    } else if (identical(covariate_level, "coarse")) {
      case <- "B"
      funs <- .as_nl_funs(nl_fun, 1L)
      amps <- .as_amps(nl_amp, 1L)
      z_a <- stats::runif(n_coarse, 0, 1)
      h_a <- as.numeric(amps[1] * funs[[1]](z_a))
      h_a <- h_a - mean(h_a)
      z <- z_a[coarse_id]
      h_true <- h_a[coarse_id]
      g_true <- h_true
      nlcovfine <- cbind(z_a = z)
      eta <- eta + h_true
    } else {
      # Case C: g(z_f) + h(z_a)
      case <- "C"
      funs <- .as_nl_funs(nl_fun, 2L)
      amps <- .as_amps(nl_amp, 2L)
      z_f <- stats::runif(n_fine, 0, 1)
      z_a <- stats::runif(n_coarse, 0, 1)
      g_true <- as.numeric(amps[1] * funs[[1]](z_f))
      g_true <- g_true - mean(g_true)
      h_a <- as.numeric(amps[2] * funs[[2]](z_a))
      h_a <- h_a - mean(h_a)
      h_true <- h_a[coarse_id]
      z <- z_f
      nlcovfine <- cbind(z_f = z_f, z_a = z_a[coarse_id])
      eta <- eta + g_true + h_true
    }
  }

  area <- as.numeric(sf::st_area(fine_geom))
  efine <- area / mean(area)

  gamma <- efine * exp(eta)
  y_fine <- stats::rpois(n_fine, lambda = pmax(gamma, 1e-12))
  y <- as.numeric(C %*% y_fine)
  mu_coarse <- as.numeric(C %*% gamma)

  stopifnot(identical(y, as.numeric(C %*% y_fine)))
  stopifnot(all.equal(sum(y), sum(y_fine)))
  stopifnot(all.equal(gamma, efine * exp(eta)))
  if (identical(case, "B")) {
    gamma0 <- efine * exp(eta_spatial)
    mu0 <- as.numeric(C %*% gamma0)
    stopifnot(all.equal(mu_coarse, mu0 * exp(h_a), tolerance = 1e-8))
  }
  if (identical(case, "C")) {
    gamma_fg <- efine * exp(eta_spatial + g_true)
    mu_fg <- as.numeric(C %*% gamma_fg)
    stopifnot(all.equal(mu_coarse, mu_fg * exp(h_a), tolerance = 1e-8))
  }

  sf_coarse <- sf::st_sf(
    id = seq_len(n_coarse),
    y = y,
    mu = mu_coarse,
    geometry = coarse_geom
  )
  if (!is.null(z_a)) {
    sf_coarse$z_a <- z_a
    sf_coarse$h_true <- h_a
  }
  sf_fine <- sf::st_sf(
    id = seq_len(n_fine),
    coarse_id = coarse_id,
    eta_true = eta,
    gamma = gamma,
    y_fine = y_fine,
    efine = efine,
    geometry = fine_geom
  )
  if (!is.null(z)) {
    sf_fine$z <- z
    if (!is.null(g_true)) sf_fine$g_true <- g_true
  }
  if (!is.null(z_f)) sf_fine$z_f <- z_f
  if (!is.null(h_true)) sf_fine$h_true <- h_true

  out <- list(
    y = y, y_fine = y_fine, C = C, x1 = x1, x2 = x2, efine = efine,
    gamma = gamma, eta_true = eta, eta_spatial_true = eta_spatial,
    mu_coarse = mu_coarse, coarse_id = coarse_id,
    n_coarse = n_coarse, n_fine = n_fine,
    sf_coarse = sf_coarse, sf_fine = sf_fine,
    geometry = "voronoi",
    case = case,
    covariate_level = if (isTRUE(include_covariate)) covariate_level else NA_character_
  )
  if (!is.null(nlcovfine)) {
    out$nlcovfine <- nlcovfine
    out$z <- z
    if (!is.null(z_f)) out$z_f <- z_f
    if (!is.null(g_true)) out$g_true <- g_true
  }
  if (!is.null(z_a)) {
    out$z_a <- z_a
    out$h_true <- h_true
    out$h_true_coarse <- h_a
  }
  out
}

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

#' Catalogue of \code{\link{simulate_ata}} scenarios
#'
#' One row per named preset (plus the no-\code{scenario} default). Columns
#' include the generative truth and the recommended \code{\link{clgam}}
#' formula / \code{orth.smooth}. Aliases: \code{"identifiability"}
#' \eqn{\rightarrow} \code{"confounding"}; \code{"matern"} \eqn{\rightarrow}
#' \code{"matern1"}.
#'
#' Cases A and C generate residual spatial truth \eqn{f_\perp=(I-P_{B(z_{\mathrm{f}})})f_{\mathrm{raw}}}
#' (same projection as \code{orth.smooth=TRUE}). Case B is unchanged
#' (coarse \eqn{h} is not in \eqn{A_{\mathrm{f}}}). Confounding remains the
#' \eqn{\rho}-mixing stress test. Pass \code{identifying="unrestricted"} to
#' reproduce the older A/C DGP \eqn{\eta=f_{\mathrm{raw}}+g}. Misspecification
#' DGPs (Matern, jump) are spatial-only.
#'
#' @return A \code{data.frame}.
#' @seealso \code{\link{simulate_ata}}, \code{\link{clgam}}, \code{\link{s}}
#' @export
#' @examples
#' simulate_ata_scenarios()[, c("scenario", "formula", "orth.smooth")]
simulate_ata_scenarios <- function() {
  data.frame(
    scenario = c(
      "default", "A", "B", "C", "confounding", "matern1", "matern2", "jump"
    ),
    n_coarse = c(12L, 40L, 40L, 40L, 40L, 40L, 40L, 40L),
    n_fine_per = c(8L, 10L, 10L, 10L, 10L, 10L, 10L, 10L),
    n_fine = c(96L, 400L, 400L, 400L, 400L, 400L, 400L, 400L),
    spatial_truth = c(
      "pspline", "pspline", "pspline", "pspline", "pspline",
      "matern", "matern", "jump"
    ),
    covariate = c(
      NA_character_, "fine", "coarse", "both", "fine",
      NA_character_, NA_character_, NA_character_
    ),
    nl_fun = c(
      "quadratic", "sine", "sine", "tanh", "sine",
      NA_character_, NA_character_, NA_character_
    ),
    spatial_amp = c(0.45, 0.8, 0.8, 0.75, 0.8, 0.8, 0.8, 0.8),
    nl_amp = c(
      "0.85", "1.2", "1.2", "1.2, 1.2", "comp_sd=0.45",
      NA_character_, NA_character_, NA_character_
    ),
    rho = c(NA_real_, NA_real_, NA_real_, NA_real_, 0, NA_real_, NA_real_, NA_real_),
    matern_nu = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, 1, 2, NA_real_),
    matern_range = c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, 0.3, 0.3, NA_real_
    ),
    jump_amp = c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, 1.6
    ),
    formula = c(
      "y ~ s(x1, x2)",
      "y ~ s(x1, x2) + s(z_f)",
      "y ~ s(x1, x2) + s(z_a, level = \"coarse\")",
      "y ~ s(x1, x2) + s(z_f) + s(z_a, level = \"coarse\")",
      "y ~ s(x1, x2) + s(z_f)",
      "y ~ s(x1, x2)",
      "y ~ s(x1, x2)",
      "y ~ s(x1, x2)"
    ),
    orth.smooth = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    notes = c(
      "No scenario=; trigonometric spatial field. Covariate off unless include_covariate=TRUE (then quadratic g(z); identifying='perp' if fine).",
      "Case A: eta = f_perp + g(z_f), z independent Uniform, sine. Fit orth.smooth=TRUE.",
      "Case B: eta = f + h(z_a), coarse sine. Restriction is a no-op; fit orth.smooth=TRUE.",
      "Case C: eta = f_perp + g(z_f) + h(z_a); f_perp vs B(z_f) only. Fit orth.smooth=TRUE.",
      "Identifiability DGP: z = rho * f_std + sqrt(1-rho^2) * eps; eta = a_f * f_perp + a_g * g(z). Override rho. Alias: identifiability.",
      "Misspecification: Matern GP, nu=1, range=0.3. Not a P-spline. Alias of scenario='matern'.",
      "Misspecification: Matern GP, nu=2, range=0.3.",
      "Misspecification: piecewise-constant jump between two coarse groups (median split on coarse centroid x1)."
    ),
    stringsAsFactors = FALSE
  )
}

#' Paper Monte Carlo presets for \code{\link{simulate_ata}}
#' @noRd
.clgam_paper_spec <- function(scenario) {
  scenario <- match.arg(
    scenario,
    c("A", "B", "C", "confounding", "identifiability",
      "matern", "matern1", "matern2", "jump")
  )
  if (identical(scenario, "identifiability")) scenario <- "confounding"
  if (identical(scenario, "matern")) scenario <- "matern1"
  base <- list(
    n_coarse = 40L,
    n_fine_per = 10L,
    include_covariate = TRUE,
    spatial_truth = "pspline"
  )
  misspec <- list(
    n_coarse = 40L,
    n_fine_per = 10L,
    include_covariate = FALSE,
    spatial_amp = 0.8,
    rho = NULL
  )
  switch(
    scenario,
    A = utils::modifyList(base, list(
      covariate_level = "fine", spatial_amp = 0.8, nl_amp = 1.2,
      nl_fun = "sine", rho = NULL
    )),
    B = utils::modifyList(base, list(
      covariate_level = "coarse", spatial_amp = 0.8, nl_amp = 1.2,
      nl_fun = "sine", rho = NULL
    )),
    C = utils::modifyList(base, list(
      covariate_level = "both", spatial_amp = 0.75, nl_amp = c(1.2, 1.2),
      nl_fun = "tanh", rho = NULL
    )),
    confounding = utils::modifyList(base, list(
      covariate_level = "fine", spatial_amp = 0.8, nl_amp = 1.2,
      nl_fun = "sine", rho = 0
    )),
    matern1 = utils::modifyList(misspec, list(
      spatial_truth = "matern", matern_nu = 1, matern_range = 0.3
    )),
    matern2 = utils::modifyList(misspec, list(
      spatial_truth = "matern", matern_nu = 2, matern_range = 0.3
    )),
    jump = utils::modifyList(misspec, list(
      spatial_truth = "jump"
    ))
  )
}

#' Named covariate-truth presets (manuscript Cases A--C / confounding)
#' @noRd
.clgam_nl_preset <- function(name) {
  key <- tolower(as.character(name)[1L])
  switch(
    key,
    sine = ,
    sin = function(z) sin(2 * pi * z),
    tanh = ,
    tanh_asc = function(z) tanh(2.5 * (z - 0.5)),
    tanh_desc = function(z) -tanh(2.5 * (z - 0.5)),
    quadratic = ,
    quad = function(z) (z - 0.5)^2,
    stop(
      "Unknown nl_fun preset '", name, "'. Use 'sine', 'tanh', 'tanh_desc', ",
      "'tanh_asc', 'quadratic', a function, or a list.",
      call. = FALSE
    )
  )
}

#' Resolve \code{nl_fun} to a list of functions of length \code{n}
#' @noRd
.clgam_resolve_nl_fun <- function(nl_fun, n) {
  n <- as.integer(n)
  if (is.null(nl_fun)) {
    nl_fun <- if (n >= 2L) "tanh" else "quadratic"
  }
  if (is.character(nl_fun)) {
    if (length(nl_fun) == 1L && n >= 2L &&
        tolower(nl_fun) %in% c("tanh", "tanh_pair")) {
      return(list(
        .clgam_nl_preset("tanh_desc"),
        .clgam_nl_preset("tanh_asc")
      )[seq_len(n)])
    }
    if (length(nl_fun) == 1L) nl_fun <- rep(nl_fun, n)
    if (length(nl_fun) < n) {
      stop("nl_fun character vector is shorter than the number of smooths.",
           call. = FALSE)
    }
    return(lapply(nl_fun[seq_len(n)], .clgam_nl_preset))
  }
  if (is.function(nl_fun)) {
    return(rep(list(nl_fun), n))
  }
  if (is.list(nl_fun)) {
    if (length(nl_fun) < n) {
      stop("nl_fun list is shorter than the number of smooths.", call. = FALSE)
    }
    out <- vector("list", n)
    for (i in seq_len(n)) {
      fi <- nl_fun[[i]]
      out[[i]] <- if (is.character(fi)) .clgam_nl_preset(fi) else fi
      if (!is.function(out[[i]])) {
        stop("nl_fun[[", i, "]] must be a function or preset name.",
             call. = FALSE)
      }
    }
    return(out)
  }
  stop("nl_fun must be a function, character preset, or list.", call. = FALSE)
}

#' Independent smooth field used as the confounding residual (script 32)
#' @noRd
.clgam_eps_field <- function(x1, x2) {
  sin(4 * pi * x1) * cos(2 * pi * x2) + 0.4 * cos(3 * pi * (x1 + x2))
}

#' Standardize to mean 0, sd 1
#' @noRd
.clgam_scale_vec <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s < 1e-12) return(x - mean(x))
  as.numeric(scale(x))
}

#' Affine map to [0, 1]
#' @noRd
.clgam_unit01 <- function(x) {
  x <- as.numeric(x)
  rng <- diff(range(x))
  (x - min(x)) / (rng + 1e-12)
}

#' Stationary isotropic Matern covariance (Stein / SPDE scaling)
#'
#' \eqn{C(h)=\sigma^2\, 2^{1-\nu}\Gamma(\nu)^{-1} (\sqrt{2\nu}\, h/\varphi)^\nu
#' K_\nu(\sqrt{2\nu}\, h/\varphi)}, with \eqn{C(0)=\sigma^2}.
#' @noRd
.clgam_matern_cov <- function(h, nu = 1, range = 0.3, sigma = 1) {
  nu <- as.numeric(nu)[1L]
  range <- as.numeric(range)[1L]
  sigma <- as.numeric(sigma)[1L]
  if (!is.finite(nu) || nu <= 0) {
    stop("matern_nu must be positive (presets use 1 or 2).", call. = FALSE)
  }
  if (!is.finite(range) || range <= 0) {
    stop("matern_range must be positive.", call. = FALSE)
  }
  h <- pmax(as.numeric(h), 0)
  out <- rep(sigma^2, length(h))
  pos <- h > 0
  if (!any(pos)) return(out)
  sc <- sqrt(2 * nu) * h[pos] / range
  const <- 2^(1 - nu) / gamma(nu)
  kv <- besselK(sc, nu)
  val <- sigma^2 * const * (sc^nu) * kv
  val[!is.finite(val)] <- sigma^2
  out[pos] <- val
  out
}

#' Draw a mean-zero Matern Gaussian field at centroids
#' @noRd
.clgam_matern_field <- function(x1, x2, nu = 1, range = 0.3, sigma = 1) {
  n <- length(x1)
  D <- as.matrix(stats::dist(cbind(as.numeric(x1), as.numeric(x2))))
  S <- matrix(
    .clgam_matern_cov(as.numeric(D), nu = nu, range = range, sigma = 1),
    n, n
  )
  jitter <- 1e-8
  U <- tryCatch(
    chol(S + diag(jitter, n)),
    error = function(e) NULL
  )
  if (is.null(U)) {
    U <- chol(S + diag(1e-5, n))
  }
  as.numeric(sigma * crossprod(U, stats::rnorm(n)))
}

#' Two coarse-group labels split by median coarse centroid in x1
#' @noRd
.clgam_jump_groups <- function(x1, coarse_id) {
  coarse_id <- as.integer(coarse_id)
  cx <- vapply(split(as.numeric(x1), coarse_id), mean, numeric(1))
  ids <- as.integer(names(cx))
  grp <- as.integer(cx > stats::median(cx))
  if (length(unique(grp)) < 2L) {
    ord <- order(ids)
    half <- max(length(ord) %/% 2L, 1L)
    grp <- integer(length(cx))
    grp[ord[seq_len(half)]] <- 0L
    grp[ord[-seq_len(half)]] <- 1L
  }
  names(grp) <- as.character(ids)
  as.integer(unname(grp[as.character(coarse_id)]))
}

#' Piecewise-constant spatial field with a jump between two coarse groups
#' @noRd
.clgam_jump_field <- function(x1, coarse_id, jump_amp = 1.6) {
  g_fine <- .clgam_jump_groups(x1, coarse_id)
  jump_amp <- as.numeric(jump_amp)[1L]
  list(
    field = (jump_amp / 2) * (2 * g_fine - 1),
    group_fine = g_fine
  )
}

#' Fine B-spline design for confounding \eqn{\kappa} / \eqn{f_\perp}
#' @noRd
.clgam_Bz <- function(z, ndxnl = 10L, bdegnl = 3L, pordnl = 2L) {
  z <- as.numeric(z)
  mm_basis(
    x = z, xl = min(z) - 0.01, xr = max(z) + 0.01,
    ndx = ndxnl, bdeg = bdegnl, pord = pordnl, decom = 2L
  )$B
}

#' Residual spatial field \eqn{f_\perp=(I-P_{B(z)})f} for identified A/C DGPs
#' @noRd
.clgam_f_perp_from_z <- function(eta_spatial, z_f, ndxnl = 10L, bdegnl = 3L) {
  f_raw <- as.numeric(eta_spatial) - mean(as.numeric(eta_spatial))
  Af <- .clgam_Bz(z_f, ndxnl = ndxnl, bdegnl = bdegnl)
  f_perp <- as.numeric(.orth_cols(cbind(f_raw), Af)[, 1L])
  kappa <- .kappa_overlap(f_raw, Af, center = TRUE)
  list(f_raw = f_raw, f_perp = f_perp, kappa = kappa)
}

#' Mix a standardized spatial field with an independent smooth residual
#' @noRd
.clgam_confound_z <- function(f_raw, x1, x2, rho) {
  rho <- as.numeric(rho)
  if (!is.finite(rho) || rho < 0 || rho > 1) {
    stop("rho must be in [0, 1].", call. = FALSE)
  }
  f_std <- .clgam_scale_vec(f_raw)
  eps <- .clgam_scale_vec(.clgam_eps_field(x1, x2))
  z <- rho * f_std + sqrt(pmax(1 - rho * rho, 0)) * eps
  .clgam_unit01(.clgam_scale_vec(z))
}

#' Simulate ATA counts on nested areal polygons (Voronoi)
#'
#' Generative story (fine \eqn{\rightarrow} coarse):
#' \enumerate{
#'   \item Build a \strong{smooth} latent field \eqn{\eta(s)} evaluated at fine
#'     centroids (piecewise constant on fine polygons): spatial surface plus
#'     optional nonlinear \eqn{g(z)}.
#'   \item Fine intensity \eqn{\gamma = e_f \odot \exp(\eta)}, with exposures
#'     \eqn{e_f} proportional to fine area.
#'   \item Simulate fine counts \eqn{y_f \sim \mathrm{Poisson}(\gamma)}
#'     (or use the mean if \code{family="none"}).
#'   \item Aggregate the \emph{raw} coarse map \eqn{y = C y_f}.
#' }
#' For independent Poisson fine counts, \eqn{y\sim\mathrm{Poisson}(C\gamma)}
#' marginally, matching the composite-link mean in \code{clgam()}.
#'
#' Geometry: coarse Voronoi of the unit square; each coarse cell is itself
#' partitioned into exactly \code{n_fine_per} fine Voronoi subpolygons
#' (a nested partition; \code{C} is 0-1 with one 1 per column).
#'
#' Named \code{scenario} values fill manuscript Monte Carlo defaults (still
#' overridable). Call \code{\link{simulate_ata_scenarios}} for the full table
#' (geometry, truth, recommended \code{\link{clgam}} formula and
#' \code{orth.smooth}).
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
#'   \eqn{z_a} so \eqn{\mu=\exp(\log(C\gamma)+h(z_a))} when \eqn{h} is within-unit
#'   constant.
#' @param spatial_amp amplitude of the spatial field on fine centroids
#'   (Matern: marginal SD; jump: half the group-level gap unless
#'   \code{jump_amp} is set)
#' @param spatial_truth spatial DGP: \code{"pspline"} (default trigonometric
#'   field, in the span of a rich anisotropic P-spline), \code{"matern"}
#'   (Gaussian process, not a P-spline), or \code{"jump"} (discontinuous
#'   across a coarse-group boundary)
#' @param matern_nu Matern smoothness \eqn{\nu} (presets 1 and 2)
#' @param matern_range Matern scale \eqn{\varphi} in the Stein parameterization
#'   \eqn{\sqrt{2\nu}\,|h|/\varphi} (unit square; default \eqn{0.3})
#' @param jump_amp difference between the two coarse-group levels. Default
#'   \code{2 * spatial_amp} so the groups sit at \eqn{\pm} \code{spatial_amp}.
#' @param nl_amp amplitude of nonlinear effect(s); length 1 or 2 for Case C
#'   as \code{c(g, h)}. Ignored under the confounding DGP (use \code{comp_sd}
#'   / \code{a_f} / \code{a_g} instead).
#' @param nl_fun covariate truth: a function of \eqn{z\in[0,1]}, a length-2
#'   \code{list} for Case C, or a preset name \code{"sine"}, \code{"tanh"}
#'   (Case C: descending then ascending), \code{"tanh_desc"},
#'   \code{"tanh_asc"}, \code{"quadratic"}. Default (no \code{scenario}) is
#'   a centred quadratic \eqn{(z-1/2)^2}.
#' @param scenario optional manuscript preset \code{"A"}, \code{"B"},
#'   \code{"C"}, \code{"confounding"} (alias \code{"identifiability"}),
#'   \code{"matern1"}, \code{"matern2"}, or \code{"jump"}.
#'   Sets geometry, amplitudes and \code{nl_fun} / \code{spatial_truth}
#'   unless the caller overrides those arguments.
#' @param rho spatial confounding in \eqn{[0,1]}. \code{NULL} (default) draws
#'   independent Uniform covariates. A value in \eqn{[0,1]} (including 0)
#'   uses the Case A mixing construction of the identifiability experiment.
#' @param intercept additive intercept on \eqn{\eta} (log intensity)
#' @param family \code{"poisson"} (default) draws \eqn{y_f\sim\mathrm{Poisson}(\gamma)};
#'   \code{"none"} uses the mean \eqn{y_f=\gamma} so \eqn{y=C\gamma} (no Poisson noise)
#' @param exposure_scale multiplies area-normalized fine exposures
#' @param comp_sd target SD for \eqn{a_f f_\perp} and \eqn{a_g g} in the confounding
#'   DGP (manuscript uses \eqn{0.45}). Used only when \code{a_f}/\code{a_g} are
#'   not supplied: they are calibrated on this geometry at \eqn{\rho=0} and then
#'   held fixed, matching the paper's scaling (shrinkage of \eqn{f_\perp} with
#'   \eqn{\rho} is geometric, not rescaled away).
#' @param a_f,a_g optional confounding scales. If both are \code{NULL} they
#'   are calibrated from \code{comp_sd} as described above.
#' @param ndxnl,bdegnl B-spline size for \eqn{A_f=\mathrm{span}\{B(z)\}} when
#'   computing \eqn{\kappa} and \eqn{f_\perp}
#' @param identifying \code{"perp"} (default): when a fine covariate is present
#'   (Case A/C), replace the spatial field by
#'   \eqn{f_\perp=(I-P_{B(z_{\mathrm{f}})})f_{\mathrm{raw}}} so the DGP matches
#'   \code{orth.smooth=TRUE}. \code{"unrestricted"} keeps \eqn{f_{\mathrm{raw}}}
#'   (older A--C recovery tables). Ignored for Case B and for the
#'   \code{rho} confounding DGP (already \eqn{f_\perp}).
#'
#' @return A list with \code{y}, \code{y_fine}, \code{C}, \code{x1}/\code{x2},
#'   \code{efine}, \code{gamma}, \code{sf_coarse}/\code{sf_fine}, and truth.
#'   Confounding and identified Case A/C draws also store \code{kappa},
#'   \code{f_raw}, \code{f_perp}. Misspecification draws store
#'   \code{spatial_truth} and, for \code{"jump"}, \code{jump_group}.
#'
#' @section Simulation scenarios:
#' \describe{
#'   \item{default (\code{scenario = NULL})}{\eqn{n=12}, \eqn{m=96}, trigonometric
#'     spatial field of amplitude \eqn{0.45}. No covariate unless
#'     \code{include_covariate=TRUE} (then centred quadratic \eqn{g(z)}).
#'     Fit \code{clgam(y ~ s(x1, x2), C = C, exposure = ef)}.}
#'   \item{\code{"A"}}{\eqn{n=40}, \eqn{m=400}, spatial amplitude \eqn{0.8}, fine sine
#'     of amplitude \eqn{1.2}. Spatial truth is \eqn{f_\perp} against
#'     \eqn{B(z_{\mathrm{f}})}. Fit
#'     \code{clgam(y ~ s(x1, x2) + s(z_f), ..., orth.smooth = TRUE)}.}
#'   \item{\code{"B"}}{Same sine on the coarse covariate \eqn{z_a} (constant
#'     within coarse units). Spatial truth remains unrestricted \eqn{f}.
#'     Fit \code{clgam(y ~ s(x1, x2) + s(z_a, level = "coarse"), ...,
#'     orth.smooth = TRUE)} (restriction is a no-op).}
#'   \item{\code{"C"}}{Spatial amplitude \eqn{0.75}, descending/ascending
#'     \eqn{\tanh} pair of amplitude \eqn{1.2}. Spatial truth is \eqn{f_\perp}
#'     against the fine \eqn{B(z_{\mathrm{f}})} only. Fit
#'     \code{clgam(y ~ s(x1, x2) + s(z_f) + s(z_a, level = "coarse"), ...,
#'     orth.smooth = TRUE)}.}
#'   \item{\code{"confounding"}}{Case A identifiability DGP
#'     \eqn{z=\rho f_{\mathrm{std}}+\sqrt{1-\rho^2}\varepsilon},
#'     \eqn{\eta=a_f f_\perp+a_g g(z)}. Default \eqn{\rho=0} (override with
#'     \code{rho=}). Alias \code{"identifiability"}. Fit
#'     \code{clgam(y ~ s(x1, x2) + s(z_f), C = C, exposure = ef,
#'     orth.smooth = TRUE)}.}
#'   \item{\code{"matern1"} / \code{"matern2"}}{Spatial truth is a Matern
#'     Gaussian field (\eqn{\nu=1} or \eqn{2}, range \eqn{0.3}), not a P-spline;
#'     no covariate. Alias \code{"matern"} \eqn{\rightarrow} \code{"matern1"}.
#'     Fit \code{clgam(y ~ s(x1, x2), ...)}.}
#'   \item{\code{"jump"}}{Piecewise-constant spatial truth with a jump
#'     between two groups of coarse units (median split on coarse centroid
#'     \eqn{x_1}). Fit \code{clgam(y ~ s(x1, x2), ...)}.}
#' }
#'
#' @seealso \code{\link{simulate_ata_scenarios}}, \code{\link{clgam}},
#'   \code{\link{s}}
#' @export
#' @examples
#' \donttest{
#' # Manuscript Case A (n=40, m=400, sine)
#' dat <- simulate_ata(scenario = "A", seed = 1)
#' # Confounding at rho = 0.6
#' dat_c <- simulate_ata(scenario = "confounding", rho = 0.6, seed = 1)
#' # Spatial truth outside the P-spline span
#' dat_m <- simulate_ata(scenario = "matern1", seed = 1)
#' dat_j <- simulate_ata(scenario = "jump", seed = 1)
#' }
simulate_ata <- function(n_coarse = 12L,
                         n_fine_per = 8L,
                         seed = 1L,
                         include_covariate = FALSE,
                         covariate_level = c("fine", "coarse", "both"),
                         spatial_amp = 0.45,
                         nl_amp = 0.85,
                         nl_fun = NULL,
                         scenario = NULL,
                         rho = NULL,
                         intercept = 0,
                         family = c("poisson", "none"),
                         exposure_scale = 1,
                         comp_sd = 0.45,
                         a_f = NULL,
                         a_g = NULL,
                         ndxnl = 10L,
                         bdegnl = 3L,
                         identifying = c("perp", "unrestricted"),
                         spatial_truth = c("pspline", "matern", "jump"),
                         matern_nu = 1,
                         matern_range = 0.3,
                         jump_amp = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("simulate_ata() requires the 'sf' package. install.packages(\"sf\")",
         call. = FALSE)
  }
  family <- match.arg(family)
  identifying <- match.arg(identifying)
  if (!is.null(scenario)) {
    spec <- .clgam_paper_spec(scenario)
    if (missing(n_coarse)) n_coarse <- spec$n_coarse
    if (missing(n_fine_per)) n_fine_per <- spec$n_fine_per
    if (missing(include_covariate)) include_covariate <- spec$include_covariate
    if (missing(covariate_level)) covariate_level <- spec$covariate_level
    if (missing(spatial_amp)) spatial_amp <- spec$spatial_amp
    if (missing(nl_amp)) nl_amp <- spec$nl_amp
    if (missing(nl_fun)) nl_fun <- spec$nl_fun
    if (missing(rho)) rho <- spec$rho
    if (missing(spatial_truth) && !is.null(spec$spatial_truth)) {
      spatial_truth <- spec$spatial_truth
    }
    if (missing(matern_nu) && !is.null(spec$matern_nu)) {
      matern_nu <- spec$matern_nu
    }
    if (missing(matern_range) && !is.null(spec$matern_range)) {
      matern_range <- spec$matern_range
    }
    if (missing(jump_amp) && !is.null(spec$jump_amp)) {
      jump_amp <- spec$jump_amp
    }
  }
  if (!is.null(rho) && missing(nl_fun)) {
    nl_fun <- "sine"
  }
  if (is.null(nl_fun)) {
    nl_fun <- function(z) (z - 0.5)^2
  }
  covariate_level <- match.arg(covariate_level, c("fine", "coarse", "both"))
  spatial_truth <- match.arg(spatial_truth, c("pspline", "matern", "jump"))
  confound <- !is.null(rho)
  if (confound) {
    include_covariate <- TRUE
    if (!identical(covariate_level, "fine")) {
      stop("rho / confounding DGP is defined for Case A (covariate_level='fine').",
           call. = FALSE)
    }
  }
  intercept <- as.numeric(intercept)[1L]
  exposure_scale <- as.numeric(exposure_scale)[1L]
  stopifnot(is.finite(intercept), is.finite(exposure_scale), exposure_scale > 0)
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

  jump_group <- NULL
  jump_group_fine <- NULL
  if (identical(spatial_truth, "matern")) {
    eta_spatial <- .clgam_matern_field(
      x1, x2, nu = matern_nu, range = matern_range, sigma = spatial_amp
    )
  } else if (identical(spatial_truth, "jump")) {
    if (is.null(jump_amp)) jump_amp <- 2 * spatial_amp
    jp <- .clgam_jump_field(x1, coarse_id, jump_amp = jump_amp)
    eta_spatial <- jp$field
    jump_group_fine <- jp$group_fine
    jump_group <- vapply(
      split(jump_group_fine, coarse_id),
      function(v) as.integer(v[1L]),
      integer(1)
    )
  } else {
    eta_spatial <- spatial_amp * sin(2 * pi * x1) * cos(2 * pi * x2) -
      0.15 * (x1 - 0.5)
  }
  eta <- eta_spatial
  nlcovfine <- NULL
  z <- NULL
  z_f <- NULL
  z_a <- NULL
  g_true <- NULL
  h_true <- NULL
  h_a <- NULL
  case <- "spatial"

  .as_amps <- function(nl_amp, n) {
    a <- rep_len(as.numeric(nl_amp), n)
    stopifnot(length(a) == n, all(is.finite(a)))
    a
  }

  kappa <- NA_real_
  rho_emp <- NA_real_
  f_raw <- NULL
  f_perp <- NULL
  a_f_used <- NA_real_
  a_g_used <- NA_real_

  if (isTRUE(include_covariate) && confound) {
    case <- "A"
    funs <- .clgam_resolve_nl_fun(nl_fun, 1L)
    f_raw <- as.numeric(eta_spatial) - mean(eta_spatial)
    z_f <- .clgam_confound_z(f_raw, x1, x2, rho)
    z <- z_f
    Af <- .clgam_Bz(z_f, ndxnl = ndxnl, bdegnl = bdegnl)
    f_perp0 <- as.numeric(.orth_cols(cbind(f_raw), Af)[, 1L])
    g0 <- as.numeric(funs[[1]](z_f))
    g0 <- g0 - mean(g0)
    kappa <- .kappa_overlap(f_raw, Af, center = TRUE)
    rho_emp <- stats::cor(z_f, f_raw)

    if (is.null(a_f) || is.null(a_g)) {
      z_ref <- .clgam_confound_z(f_raw, x1, x2, 0)
      Af_ref <- .clgam_Bz(z_ref, ndxnl = ndxnl, bdegnl = bdegnl)
      f_ref <- as.numeric(.orth_cols(cbind(f_raw), Af_ref)[, 1L])
      g_ref <- as.numeric(funs[[1]](z_ref))
      g_ref <- g_ref - mean(g_ref)
      sd_f_ref <- stats::sd(f_ref)
      sd_g_ref <- stats::sd(g_ref)
      if (is.null(a_f)) {
        a_f <- as.numeric(comp_sd) / max(sd_f_ref, 1e-8)
      }
      if (is.null(a_g)) {
        a_g <- as.numeric(comp_sd) / max(sd_g_ref, 1e-8)
      }
    }
    a_f_used <- as.numeric(a_f)[1L]
    a_g_used <- as.numeric(a_g)[1L]
    f_perp <- a_f_used * f_perp0
    g_true <- a_g_used * g0
    eta_spatial <- f_perp
    nlcovfine <- cbind(z_f = z_f)
    eta <- eta_spatial + g_true
  } else if (isTRUE(include_covariate)) {
    if (identical(covariate_level, "fine")) {
      case <- "A"
      funs <- .clgam_resolve_nl_fun(nl_fun, 1L)
      amps <- .as_amps(nl_amp, 1L)
      z_f <- stats::runif(n_fine, 0, 1)
      z <- z_f
      g_true <- as.numeric(amps[1] * funs[[1]](z_f))
      g_true <- g_true - mean(g_true)
      nlcovfine <- cbind(z_f = z_f)
      if (identical(identifying, "perp")) {
        proj <- .clgam_f_perp_from_z(eta_spatial, z_f, ndxnl, bdegnl)
        f_raw <- proj$f_raw
        f_perp <- proj$f_perp
        kappa <- proj$kappa
        eta_spatial <- f_perp
      }
      eta <- eta_spatial + g_true
    } else if (identical(covariate_level, "coarse")) {
      case <- "B"
      funs <- .clgam_resolve_nl_fun(nl_fun, 1L)
      amps <- .as_amps(nl_amp, 1L)
      z_a <- stats::runif(n_coarse, 0, 1)
      h_a <- as.numeric(amps[1] * funs[[1]](z_a))
      h_a <- h_a - mean(h_a)
      z <- z_a[coarse_id]
      h_true <- h_a[coarse_id]
      g_true <- h_true
      nlcovfine <- cbind(z_a = z)
      eta <- eta_spatial + h_true
    } else {
      # Case C: g(z_f) + h(z_a). Project f only vs fine B(z_f).
      case <- "C"
      funs <- .clgam_resolve_nl_fun(nl_fun, 2L)
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
      if (identical(identifying, "perp")) {
        proj <- .clgam_f_perp_from_z(eta_spatial, z_f, ndxnl, bdegnl)
        f_raw <- proj$f_raw
        f_perp <- proj$f_perp
        kappa <- proj$kappa
        eta_spatial <- f_perp
      }
      eta <- eta_spatial + g_true + h_true
    }
  }

  area <- as.numeric(sf::st_area(fine_geom))
  efine <- exposure_scale * area / mean(area)
  eta <- eta + intercept

  gamma <- efine * exp(eta)
  if (identical(family, "none")) {
    y_fine <- as.numeric(gamma)
  } else {
    y_fine <- stats::rpois(n_fine, lambda = pmax(gamma, 1e-12))
  }
  y <- as.numeric(C %*% y_fine)
  mu_coarse <- as.numeric(C %*% gamma)

  stopifnot(isTRUE(all.equal(y, as.numeric(C %*% y_fine))))
  stopifnot(all.equal(sum(y), sum(y_fine)))
  stopifnot(all.equal(gamma, efine * exp(eta)))
  if (identical(case, "B")) {
    gamma0 <- efine * exp(eta_spatial + intercept)
    mu0 <- as.numeric(C %*% gamma0)
    stopifnot(all.equal(mu_coarse, mu0 * exp(h_a), tolerance = 1e-8))
  }
  if (identical(case, "C")) {
    gamma_fg <- efine * exp(eta_spatial + g_true + intercept)
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
  if (!is.null(jump_group)) {
    sf_coarse$jump_group <- as.integer(jump_group)
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
  if (!is.null(jump_group_fine)) sf_fine$jump_group <- as.integer(jump_group_fine)

  out <- list(
    y = y, y_fine = y_fine, C = C, x1 = x1, x2 = x2, efine = efine,
    gamma = gamma, eta_true = eta, eta_spatial_true = eta_spatial,
    mu_coarse = mu_coarse, coarse_id = coarse_id,
    n_coarse = n_coarse, n_fine = n_fine,
    sf_coarse = sf_coarse, sf_fine = sf_fine,
    geometry = "voronoi",
    case = case,
    covariate_level = if (isTRUE(include_covariate)) covariate_level else NA_character_,
    nl_level = NULL,
    scenario = if (is.null(scenario)) NA_character_ else as.character(scenario)[1L],
    intercept = intercept,
    family = family,
    exposure_scale = exposure_scale,
    rho = if (confound) as.numeric(rho) else NA_real_,
    rho_emp = rho_emp,
    kappa = kappa,
    a_f = a_f_used,
    a_g = a_g_used,
    spatial_truth = spatial_truth,
    matern_nu = if (identical(spatial_truth, "matern")) as.numeric(matern_nu) else NA_real_,
    matern_range = if (identical(spatial_truth, "matern")) as.numeric(matern_range) else NA_real_,
    jump_amp = if (identical(spatial_truth, "jump")) as.numeric(jump_amp) else NA_real_,
    identifying = if (confound) "perp" else identifying
  )
  if (!is.null(jump_group)) out$jump_group <- as.integer(jump_group)
  if (!is.null(f_raw)) out$f_raw <- f_raw
  if (!is.null(f_perp)) out$f_perp <- f_perp
  if (!is.null(nlcovfine)) {
    out$nlcovfine <- nlcovfine
    out$z <- z
    if (!is.null(z_f)) out$z_f <- z_f
    if (!is.null(g_true)) out$g_true <- g_true
    out$nl_level <- switch(
      case,
      A = "fine",
      B = "coarse",
      C = c("fine", "coarse"),
      rep("fine", ncol(nlcovfine))
    )
  }
  if (!is.null(z_a)) {
    out$z_a <- z_a
    out$h_true <- h_true
    out$h_true_coarse <- h_a
  }
  out
}

# Malone / dissever downscaling — Method 1 + optional revive intermediates.
#
# See experiments/improvements/MALONE_ADAPTATIONS.md

#' Fine init for nested 0/1 C.
#' @param type `"equal"` = C' (y/n_i) (Method 1); `"exposure"` = y_i * e_j / sum_{j in i} e
#'   (revive intermediate — not Method 1).
clgam_malone_init <- function(y_mun, C, e_fine = NULL, type = c("equal", "exposure")) {
  type <- match.arg(type)
  C <- as.matrix(C)
  stopifnot(all(abs(colSums(C) - 1) < 1e-8))
  if (identical(type, "equal")) {
    n_per <- rowSums(C)
    stopifnot(all(n_per > 0))
    return(as.numeric(t(C) %*% (y_mun / n_per)))
  }
  if (is.null(e_fine)) stop("e_fine required for init='exposure'")
  e_fine <- as.numeric(e_fine)
  e_mun <- as.numeric(C %*% e_fine)
  stopifnot(all(e_mun > 0))
  # y_j = e_j * (y_i / e_i) for j in i
  as.numeric(e_fine * as.vector(t(C) %*% (y_mun / e_mun)))
}

#' Mass-balance rescale: y_f[j] *= y_mun[i] / sum_{j' in i} y_f[j'].
#' (Equivalent to rate rescale r <- r * y_i/sum(e*r) when y = e*r.)
clgam_malone_rescale <- function(y_fine, y_mun, C, eps = 1e-12) {
  C <- as.matrix(C)
  agg <- as.numeric(C %*% y_fine)
  w_mun <- y_mun / pmax(agg, eps)
  as.numeric(y_fine * as.vector(t(C) %*% w_mun))
}

#' Malone Method 1 (+ optional intermediate revive flags).
#'
#' @param init `"equal"` (Method 1) or `"exposure"` (revive: mun SMR × e_f).
#' @param deliver `"fitted"` (GAM mu / eta) or `"balanced"` (post-rescale map
#'   with C%*%y = y_mun exactly — revive intermediate for evaluation).
#' @param stop_on `"map"` = mean|μ_t-μ_{t-1}|; `"balanced_map"` = mean|bal_t-bal_{t-1}|;
#'   `"mass"` = max|C%*%μ - y| (rarely met for unconstrained GAM).
#' @param family,round_counts,spatial See MALONE_ADAPTATIONS.md
clgam_malone_fit <- function(
  y_mun,
  e_fine,
  C,
  cov_fine,
  lon = NULL,
  lat = NULL,
  cov_smooth = c("ageing", "unemployed"),
  spatial = c("none", "additive", "te"),
  family = c("poisson", "quasipoisson"),
  round_counts = TRUE,
  init = c("equal", "exposure"),
  deliver = c("fitted", "balanced"),
  stop_on = c("map", "balanced_map", "mass"),
  max_iter = 300L,
  tol = 1e-3,
  tol_mass = 1,
  k_spatial = 80L,
  k_cov = 10L,
  trace = TRUE
) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for Malone downscaling.")
  }
  spatial <- match.arg(spatial)
  family <- match.arg(family)
  init <- match.arg(init)
  deliver <- match.arg(deliver)
  stop_on <- match.arg(stop_on)
  if (identical(family, "poisson") && !isTRUE(round_counts)) {
    warning("Poisson+REML needs integer counts; setting round_counts = TRUE")
    round_counts <- TRUE
  }
  if (spatial != "none" && (is.null(lon) || is.null(lat))) {
    stop("lon/lat required when spatial != 'none'")
  }

  fam <- if (identical(family, "poisson")) {
    stats::poisson(link = "log")
  } else {
    stats::quasipoisson(link = "log")
  }

  C <- as.matrix(C)
  cov_fine <- as.data.frame(cov_fine)
  e_fine <- as.numeric(e_fine)
  m <- ncol(C)
  stopifnot(
    length(y_mun) == nrow(C),
    length(e_fine) == m,
    nrow(cov_fine) == m
  )
  miss <- setdiff(cov_smooth, names(cov_fine))
  if (length(miss)) stop("Missing cov columns: ", paste(miss, collapse = ", "))

  offset <- log(pmax(e_fine, .Machine$double.xmin))
  df <- data.frame(
    offset = offset,
    cov_fine[, cov_smooth, drop = FALSE],
    check.names = FALSE
  )
  if (spatial != "none") {
    df$lon <- lon
    df$lat <- lat
  }

  rhs_cov <- paste(sprintf("s(%s, k = %d)", cov_smooth, as.integer(k_cov)), collapse = " + ")
  rhs_sp <- switch(
    spatial,
    none = NULL,
    additive = sprintf("s(lon, k = %d) + s(lat, k = %d)", as.integer(k_cov), as.integer(k_cov)),
    te = sprintf("s(lon, lat, k = %d)", as.integer(k_spatial))
  )
  rhs <- paste(c(rhs_cov, rhs_sp), collapse = " + ")
  fml <- stats::as.formula(paste("y_work ~", rhs, "+ offset(offset)"))

  adaptations <- c(
    "areal nested composition C (Malone 2012: regular raster)",
    "Poisson/quasi log-link + offset log(e_fine) (Malone 2012: continuous target)"
  )
  if (identical(init, "equal")) {
    adaptations <- c(adaptations, "init equal-split C^+ y (CL-GAMM Method 1)")
  } else {
    adaptations <- c(
      adaptations,
      "REVIVE: init exposure-proportional y_j = e_j*(y_i/e_i) (not Method 1)"
    )
  }
  if (isTRUE(round_counts)) {
    adaptations <- c(adaptations, "round working counts (CL-GAMM Method 1)")
  }
  if (identical(family, "quasipoisson")) {
    adaptations <- c(adaptations, "REVIVE: quasipoisson, typically without round")
  }
  if (identical(deliver, "balanced")) {
    adaptations <- c(
      adaptations,
      "REVIVE: deliver mass-balanced map (C%*%y=y_mun), not raw GAM fitted"
    )
  }
  if (identical(stop_on, "balanced_map")) {
    adaptations <- c(adaptations, "REVIVE: stop on successive balanced maps")
  }
  if (identical(stop_on, "mass")) {
    adaptations <- c(adaptations, "REVIVE: stop on max|C%*%mu - y| (mass of fitted)")
  }
  if (identical(spatial, "additive")) {
    adaptations <- c(adaptations, "spatial s(lon)+s(lat) (Codes3)")
  }
  if (identical(spatial, "te")) {
    adaptations <- c(adaptations, "spatial s(lon,lat) (revive option)")
  }
  bugfixes <- "mass-balance via t(C) %*% w (order-safe)"

  y_work <- clgam_malone_init(y_mun, C, e_fine = e_fine, type = init)
  if (round_counts) y_work <- pmax(0, round(y_work))

  fit_once <- function(y) {
    d <- df
    d$y_work <- if (round_counts) pmax(0, round(y)) else pmax(as.numeric(y), 0)
    mgcv::gam(fml, family = fam, data = d, method = "REML")
  }

  if (trace) {
    message(sprintf(
      "Malone: family=%s round=%s init=%s deliver=%s stop_on=%s spatial=%s",
      family, round_counts, init, deliver, stop_on, spatial
    ))
  }

  model <- fit_once(y_work)
  mu_prev <- as.numeric(stats::fitted(model))
  bal_prev <- clgam_malone_rescale(mu_prev, y_mun, C)
  eta <- as.numeric(model$linear.predictors - offset)

  n_iter <- 0L
  crit_map <- crit_balance <- crit_fitgap <- crit_bal_map <- mass_fit <- NA_real_
  hist <- list()
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    n_iter <- iter
    y_bal <- clgam_malone_rescale(mu_prev, y_mun, C)
    crit_balance <- mean(abs(y_bal - mu_prev))

    model <- fit_once(y_bal)
    mu_new <- as.numeric(stats::fitted(model))
    bal_new <- clgam_malone_rescale(mu_new, y_mun, C)

    crit_map <- mean(abs(mu_new - mu_prev))
    crit_fitgap <- mean(abs(mu_new - y_bal))
    crit_bal_map <- mean(abs(bal_new - bal_prev))
    mass_fit <- max(abs(as.numeric(C %*% mu_new) - y_mun))

    if (trace) {
      message(sprintf(
        "  iter %03d  |Δmu|=%.4g  |Δbal|=%.4g  |mu-bal|=%.4g  max|C mu-y|=%.4g",
        iter, crit_map, crit_bal_map, crit_fitgap, mass_fit
      ))
    }
    hist[[iter]] <- data.frame(
      iter = iter,
      crit_map = crit_map,
      crit_bal_map = crit_bal_map,
      crit_balance = crit_balance,
      crit_fitgap = crit_fitgap,
      mass_fitted = mass_fit
    )

    mu_prev <- mu_new
    bal_prev <- bal_new
    eta <- as.numeric(model$linear.predictors - offset)

    hit <- switch(
      stop_on,
      map = is.finite(crit_map) && crit_map < tol,
      balanced_map = is.finite(crit_bal_map) && crit_bal_map < tol,
      mass = is.finite(mass_fit) && mass_fit < tol_mass
    )
    if (isTRUE(hit)) {
      converged <- TRUE
      break
    }
  }

  y_bal_final <- clgam_malone_rescale(mu_prev, y_mun, C)
  # GAM eta always available; balanced eta only where y_bal > 0
  eta_gam <- as.numeric(model$linear.predictors - offset)
  eta_bal <- log(
    pmax(y_bal_final, .Machine$double.xmin) / pmax(e_fine, .Machine$double.xmin)
  )
  # Prefer GAM linear predictor on tracts where balanced count is tiny
  # (log(y_bal/e) is unstable); keep balanced counts for mass-exact delivery.
  tiny <- y_bal_final < 1e-8
  eta_bal[tiny] <- eta_gam[tiny]

  if (identical(deliver, "balanced")) {
    eta_out <- eta_bal
    fitted_out <- y_bal_final
  } else {
    eta_out <- eta_gam
    fitted_out <- mu_prev
  }

  list(
    eta = eta_out,
    fitted = fitted_out,
    fitted_gam = mu_prev,
    fitted_balanced = y_bal_final,
    eta_gam = eta_gam,
    eta_balanced = eta_bal,
    model = model,
    n_iter = n_iter,
    crit = crit_map,
    crit_map = crit_map,
    crit_bal_map = crit_bal_map,
    crit_balance = crit_balance,
    crit_fitgap = crit_fitgap,
    tol = tol,
    tol_mass = tol_mass,
    converged = converged,
    mass_error_fitted = max(abs(as.numeric(C %*% mu_prev) - y_mun)),
    mass_error_balanced = max(abs(as.numeric(C %*% y_bal_final) - y_mun)),
    mass_error_delivered = max(abs(as.numeric(C %*% fitted_out) - y_mun)),
    spatial = spatial,
    family = family,
    round_counts = round_counts,
    init = init,
    deliver = deliver,
    stop_on = stop_on,
    cov_smooth = cov_smooth,
    formula = fml,
    history = if (length(hist)) do.call(rbind, hist) else NULL,
    adaptations = adaptations,
    bugfixes = bugfixes,
    reference = "Malone et al. (2012); CL-GAMM Method 1; revive intermediates flagged REVIVE"
  )
}

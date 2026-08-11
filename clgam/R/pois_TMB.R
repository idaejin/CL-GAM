#' Poisson CL-GAMM via TMB Laplace approximation (\code{pois_TMB})
#'
#' Same mixed-model design as \code{\link{pois_SOP}} (anisotropic spatial
#' P-splines + optional fine-scale univariate smooths), but integrates the
#' random effects with a Laplace approximation to the marginal Poisson
#' composite-link likelihood via \pkg{TMB}. Intended as a robustness check
#' against PIRLS+SOP (PQL-style) for simulations and coverage.
#'
#' Requires \pkg{TMB} (Suggests). The template lives in \code{inst/TMB/} and is
#' compiled at runtime (separate DLL from the Rcpp SOP kernels).
#'
#' @inheritParams pois_SOP
#' @param silent suppress TMB / optimizer chatter
#' @param control passed to \code{\link[stats]{nlminb}}
#' @return An object of class \code{"clgam"} with \code{family="spatial_tmb"},
#'   including \code{eta}, \code{sd.eta} (from \code{sdreport} / ADREPORT),
#'   \code{var.comp} (\(\tau^2\)), and \code{tmb} fit bits.
#' @export
pois_TMB <- function(
  y, x1, x2,
  efine = NULL,
  lcovfine = NULL,
  nlcovfine = NULL,
  C,
  x1lim = NULL,
  x2lim = NULL,
  ndx = c(15, 15),
  bdeg = c(3, 3),
  pord = c(2, 2),
  decom = 1,
  parold = c(1, 1),
  paroldnl = NULL,
  ndxnl = 15,
  bdegnl = 3,
  pordnl = 2,
  nl.basis = c("pspline", "legacy"),
  orth.smooth = NULL,
  sparse.backend = "auto",
  silent = TRUE,
  control = list(iter.max = 200, eval.max = 400)
) {
  nl.basis <- match.arg(nl.basis)
  if (is.null(orth.smooth)) orth.smooth <- identical(nl.basis, "pspline")
  orth.smooth <- isTRUE(orth.smooth)
  start.all <- proc.time()[3]

  des <- .clgam_mm_design(
    x1 = x1, x2 = x2, efine = efine,
    lcovfine = lcovfine, nlcovfine = nlcovfine, C = C,
    x1lim = x1lim, x2lim = x2lim,
    ndx = ndx, bdeg = bdeg, pord = pord, decom = decom,
    ndxnl = ndxnl, bdegnl = bdegnl, pordnl = pordnl,
    nl.basis = nl.basis, orth.smooth = orth.smooth,
    sparse.backend = sparse.backend
  )

  if (is.null(des$groups)) {
    stop("pois_TMB currently requires a nested 0-1 partition C ",
         "(one fine unit per coarse cell membership).", call. = FALSE)
  }

  n_tau <- 2L + if (!is.null(nlcovfine)) ncol(as.matrix(nlcovfine)) else 0L
  la0 <- as.numeric(parold)
  if (length(la0) < 2L) la0 <- c(la0, rep(1, 2L - length(la0)))
  if (n_tau > 2L) {
    if (is.null(paroldnl)) paroldnl <- rep(1, n_tau - 2L)
    la0 <- c(la0[1:2], as.numeric(paroldnl))
  } else {
    la0 <- la0[1:2]
  }

  ensure_tmb_dll("clgam_pois")

  data <- list(
    y = as.numeric(y),
    X = as.matrix(des$X),
    Z = as.matrix(des$Z),
    e = as.numeric(des$efine),
    groups = as.integer(des$groups) - 1L,
    n_coarse = as.integer(des$n_coarse),
    G1inv_n = as.numeric(des$G1inv_n),
    G2inv_n = as.numeric(des$G2inv_n),
    Gnl_n = as.numeric(des$Gnl_n)
  )
  parameters <- list(
    beta = rep(0, ncol(des$X)),
    u = rep(0, ncol(des$Z)),
    log_tau2 = log(pmax(la0, 1e-8))
  )

  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    random = "u",
    DLL = "clgam_pois",
    silent = isTRUE(silent)
  )

  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = control,
    # Keep log(tau^2) in a stable range (avoid Laplace boundary collapse)
    lower = ifelse(names(obj$par) == "log_tau2", log(1e-4), -Inf),
    upper = ifelse(names(obj$par) == "log_tau2", log(1e4), Inf)
  )

  # Hyperparameters / fixed effects live in opt$par; random effects in last.par.best
  b.fixed <- as.numeric(opt$par[names(opt$par) == "beta"])
  log_tau2 <- as.numeric(opt$par[names(opt$par) == "log_tau2"])
  la <- exp(log_tau2)
  # Ensure inner RE solve is at optimized hyperparameters before sdreport
  obj$fn(opt$par)
  b.random <- as.numeric(obj$env$last.par.best[names(obj$env$last.par.best) == "u"])
  eta <- as.numeric(des$X %*% b.fixed + des$Z %*% b.random)
  eta <- pmin(pmax(eta, -20), 20)

  # sdreport at optimized hyperparameters
  rep <- tryCatch(
    TMB::sdreport(obj, getJointPrecision = FALSE),
    error = function(e) {
      warning("pois_TMB: sdreport failed: ", conditionMessage(e))
      NULL
    }
  )
  sd.eta <- rep(NA_real_, length(eta))
  if (!is.null(rep)) {
    sdr_sum <- summary(rep)
    eta_rows <- which(rownames(sdr_sum) == "eta")
    if (length(eta_rows) == length(eta)) {
      sd.eta <- as.numeric(sdr_sum[eta_rows, "Std. Error"])
    } else {
      warning("pois_TMB: ADREPORT(eta) length mismatch in sdreport()")
    }
  }

  gamma <- as.numeric(des$efine * exp(eta))
  mu <- as.numeric(rowsum(
    gamma,
    group = factor(des$groups, levels = seq_len(des$n_coarse))
  ))
  # Coarse deviance
  eps <- .Machine$double.xmin
  mu_s <- pmax(mu, eps)
  y_n <- as.numeric(y)
  sat <- ifelse(y_n > 0, y_n * log(y_n) - y_n, 0)
  fitll <- y_n * log(mu_s) - mu_s
  dev <- 2 * sum(sat - fitll)

  # Rough ED via trace of random-effect conditional precision is expensive;
  # report length(beta) + effective from TMB as NA for now (AIC not primary).
  ed <- NA_real_
  aic <- NA_real_
  bic <- NA_real_

  end.all <- proc.time()[3]
  out <- list(
    ndx = ndx, bdeg = bdeg, pord = pord,
    knots1 = des$MM1$knots, knots2 = des$MM2$knots,
    y = as.numeric(y), x1 = x1, x2 = x2, efine = des$efine,
    lcovfine = lcovfine, nlcovfine = nlcovfine,
    eta = eta, gamma = gamma, mu = mu,
    var.comp = la, edf = NA_real_, niter = as.integer(opt$iterations),
    elapsed.time = end.all - start.all,
    dev = dev, b.fixed = b.fixed, b.random = b.random,
    matlist = list(
      X = des$X, Z = des$Z, C = des$C,
      G1inv_n = des$G1inv_n, G2inv_n = des$G2inv_n, Gnl_n = des$Gnl_n
    ),
    leffects = NULL, sdleffects = NULL,
    nleffects = NULL, sdnleffects = NULL,
    nl.basis = nl.basis, orth.smooth = orth.smooth,
    ed = ed, aic = aic, bic = bic,
    sd.eta = sd.eta,
    sd.exp.eta = sd.eta * exp(eta),
    tmb = list(opt = opt, obj = obj, sdreport = rep, convergence = opt$convergence),
    engine = "pois_TMB"
  )
  .as_clgam(out, call = match.call(), family = "spatial")
}

#' Build mixed-model design (shared by SOP / TMB)
#' @keywords internal
.clgam_mm_design <- function(
  x1, x2, efine = NULL, lcovfine = NULL, nlcovfine = NULL, C,
  x1lim = NULL, x2lim = NULL,
  ndx = c(15, 15), bdeg = c(3, 3), pord = c(2, 2), decom = 1,
  ndxnl = 15, bdegnl = 3, pordnl = 2,
  nl.basis = "pspline", orth.smooth = TRUE,
  sparse.backend = "auto"
) {
  if (is.null(x1lim)) x1lim <- .clgam_xlim(x1)
  if (is.null(x2lim)) x2lim <- .clgam_xlim(x2)
  dimfine <- length(x1)
  efine <- .clgam_exposure(efine, dimfine)

  MM1 <- mm_basis(x = x1, xl = x1lim[1], xr = x1lim[2], ndx = ndx[1],
                  bdeg = bdeg[1], pord = pord[1], decom = decom)
  MM2 <- mm_basis(x = x2, xl = x2lim[1], xr = x2lim[2], ndx = ndx[2],
                  bdeg = bdeg[2], pord = pord[2], decom = decom)
  X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; c1 <- MM1$m
  X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; c2 <- MM2$m

  g1u <- rep(d1, times = pord[2])
  g2u <- rep(d2, each = pord[1])
  g1b <- rep(d1, times = (c2 - pord[2]))
  g2b <- rep(d2, each = (c1 - pord[1]))

  nk <- 0L
  Xxk <- Zxk <- Bk <- dk <- NULL
  ck <- integer(0)
  if (!is.null(nlcovfine)) {
    nlcovfine <- as.matrix(nlcovfine)
    nk <- ncol(nlcovfine)
    if (length(ndxnl) == 1) ndxnl <- rep(ndxnl, nk)
    if (length(bdegnl) == 1) bdegnl <- rep(bdegnl, nk)
    if (length(pordnl) == 1) pordnl <- rep(pordnl, nk)
    Xxk <- vector("list", nk)
    Zxk <- vector("list", nk)
    Bk <- vector("list", nk)
    dk <- vector("list", nk)
    ck <- rep(0L, nk)
    for (k in seq_len(nk)) {
      xk <- nlcovfine[, k]
      xklim <- c(min(xk) - 0.01, max(xk) + 0.01)
      decom_k <- if (identical(nl.basis, "pspline")) 2L else decom
      MMk <- mm_basis(
        x = xk, xl = xklim[1], xr = xklim[2], ndx = ndxnl[k],
        bdeg = bdegnl[k], pord = pordnl[k], decom = decom_k
      )
      Xxk[[k]] <- MMk$X
      Zxk[[k]] <- MMk$Z
      Bk[[k]] <- MMk$B
      dk[[k]] <- MMk$d
      ck[k] <- MMk$m
    }
  }

  np <- c(
    prod(pord),
    (c2 - pord[2]) * pord[1],
    (c1 - pord[1]) * pord[2],
    (c1 - pord[1]) * (c2 - pord[2])
  )
  if (!is.null(lcovfine)) np[1] <- np[1] + ncol(as.matrix(lcovfine))
  if (nk > 0L) {
    if (identical(nl.basis, "pspline")) {
      for (k in seq_len(nk)) np[1] <- np[1] + max(ncol(Xxk[[k]]) - 1L, 0L)
    } else {
      np[1] <- np[1] + nk
    }
    for (k in seq_len(nk)) np <- c(np, (ck[k] - pordnl[k]))
  }

  G1inv_n <- c(rep(0, np[2]), g1u, g1b)
  G2inv_n <- c(g2u, rep(0, np[3]), g2b)
  Gnl_n <- rep(0, sum(np[-1]))
  if (nk > 0L) {
    G1inv_n <- c(G1inv_n, rep(0, sum(np[-(1:4)])))
    G2inv_n <- c(G2inv_n, rep(0, sum(np[-(1:4)])))
    for (k in seq_len(nk)) {
      low <- sum(np[2:(k + 3)]) + 1L
      sup <- low + np[k + 4] - 1L
      Gnl_n[low:sup] <- Gnl_n[low:sup] + dk[[k]]
    }
  }

  X <- rten2(X2, X1)
  Z <- cbind(rten2(Z2, X1), rten2(X2, Z1), rten2(Z2, Z1))
  n_sp_fixed <- ncol(X)
  n_sp_random <- ncol(Z)

  if (!is.null(lcovfine)) X <- cbind(X, as.matrix(lcovfine))

  if (nk > 0L) {
    if (isTRUE(orth.smooth)) {
      Bsm <- do.call(cbind, Bk)
      if (n_sp_fixed > 1L) {
        X[, 2:n_sp_fixed] <- .orth_cols(X[, 2:n_sp_fixed, drop = FALSE], Bsm)
      }
      Z[, seq_len(n_sp_random)] <- .orth_cols(
        Z[, seq_len(n_sp_random), drop = FALSE], Bsm
      )
    }
    if (identical(nl.basis, "pspline")) {
      for (k in seq_len(nk)) {
        Xk <- Xxk[[k]]
        if (ncol(Xk) > 1L) Xk <- Xk[, -1L, drop = FALSE] else Xk <- Xk[, FALSE, drop = FALSE]
        if (ncol(Xk) > 0L) X <- cbind(X, Xk)
        Z <- cbind(Z, Zxk[[k]])
      }
    } else {
      X <- cbind(X, nlcovfine)
      for (k in seq_len(nk)) Z <- cbind(Z, Zxk[[k]])
    }
  }

  stopifnot(length(G1inv_n) == ncol(Z), length(G2inv_n) == ncol(Z),
            length(Gnl_n) == ncol(Z))

  C <- .as_comp_C(C, backend = sparse.backend)
  groups <- if (.is_partition_C(C)) .partition_groups(C) else NULL
  n_coarse <- nrow(C)

  list(
    X = X, Z = Z, efine = efine, C = C, groups = groups, n_coarse = n_coarse,
    G1inv_n = G1inv_n, G2inv_n = G2inv_n, Gnl_n = Gnl_n,
    np = np, MM1 = MM1, MM2 = MM2
  )
}

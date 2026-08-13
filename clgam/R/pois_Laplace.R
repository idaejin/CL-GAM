#' Mixed-model design for Laplace (same X/Z/Ginv layout as \code{pois_SOP})
#' @keywords internal
.clgam_mixed_design <- function(
  x1, x2, efine = NULL, lcovfine = NULL, nlcovfine = NULL, C,
  x1lim = NULL, x2lim = NULL,
  ndx = c(15, 15), bdeg = c(3, 3), pord = c(2, 2), decom = 1,
  ndxnl = 15, bdegnl = 3, pordnl = 2,
  nl.basis = "pspline", orth.smooth = TRUE, nl.level = NULL,
  sparse.backend = "auto",
  re = "none"
) {
  if (is.null(x1lim)) x1lim <- .clgam_xlim(x1)
  if (is.null(x2lim)) x2lim <- .clgam_xlim(x2)
  dimfine <- length(x1)
  efine <- .clgam_exposure(efine, dimfine)
  re <- .clgam_resolve_re(re)

  MM1 <- mm_basis(
    x = x1, xl = x1lim[1], xr = x1lim[2], ndx = ndx[1],
    bdeg = bdeg[1], pord = pord[1], decom = decom
  )
  MM2 <- mm_basis(
    x = x2, xl = x2lim[1], xr = x2lim[2], ndx = ndx[2],
    bdeg = bdeg[2], pord = pord[2], decom = decom
  )
  X1 <- MM1$X
  Z1 <- MM1$Z
  d1 <- MM1$d
  c1 <- MM1$m
  X2 <- MM2$X
  Z2 <- MM2$Z
  d2 <- MM2$d
  c2 <- MM2$m

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
  Gk <- NULL
  if (nk > 0L) {
    G1inv_n <- c(G1inv_n, rep(0, sum(np[-(1:4)])))
    G2inv_n <- c(G2inv_n, rep(0, sum(np[-(1:4)])))
    Gk <- matrix(0, length(G1inv_n), nk)
    for (k in seq_len(nk)) {
      low <- sum(np[2:(k + 3)]) + 1L
      sup <- low + np[k + 4] - 1L
      Gk[low:sup, k] <- dk[[k]]
    }
  }

  X <- rten2(X2, X1)
  Z <- cbind(rten2(Z2, X1), rten2(X2, Z1), rten2(Z2, Z1))
  n_sp_fixed <- ncol(X)
  n_sp_random <- ncol(Z)

  nl_level_resolved <- character(0)
  if (nk > 0L) nl_level_resolved <- .resolve_nl_level(nk, nl.level)
  Af <- .build_orth_Af(
    Bk = if (nk > 0L) Bk else NULL,
    nl.level = nl_level_resolved,
    lcovfine = lcovfine
  )
  if (isTRUE(orth.smooth) && !is.null(Af) && ncol(Af) > 0L) {
    if (n_sp_fixed > 1L) {
      X[, 2:n_sp_fixed] <- .orth_cols(X[, 2:n_sp_fixed, drop = FALSE], Af)
    }
    Z[, seq_len(n_sp_random)] <- .orth_cols(
      Z[, seq_len(n_sp_random), drop = FALSE], Af
    )
  }

  if (!is.null(lcovfine)) X <- cbind(X, as.matrix(lcovfine))

  if (nk > 0L) {
    if (identical(nl.basis, "pspline")) {
      for (k in seq_len(nk)) {
        Xk <- Xxk[[k]]
        if (ncol(Xk) > 1L) {
          Xk <- Xk[, -1L, drop = FALSE]
        } else {
          Xk <- Xk[, FALSE, drop = FALSE]
        }
        if (ncol(Xk) > 0L) X <- cbind(X, Xk)
        Z <- cbind(Z, Zxk[[k]])
      }
    } else {
      X <- cbind(X, nlcovfine)
      for (k in seq_len(nk)) Z <- cbind(Z, Zxk[[k]])
    }
  }

  C <- .as_comp_C(C, backend = sparse.backend)
  groups <- if (.is_partition_C(C)) .partition_groups(C) else NULL
  n_coarse <- nrow(C)
  n_re <- 0L
  Z_re <- NULL
  if (identical(re, "coarse")) {
    if (is.null(groups)) {
      stop(
        "re='coarse' requires a nested 0-1 partition C.",
        call. = FALSE
      )
    }
    n_re <- n_coarse
    Z_re <- .clgam_Z_re(groups, n_coarse)
    q0 <- ncol(Z)
    Z <- cbind(Z, Z_re)
    G1inv_n <- c(G1inv_n, rep(0, n_re))
    G2inv_n <- c(G2inv_n, rep(0, n_re))
    if (!is.null(Gk)) {
      Gk <- rbind(Gk, matrix(0, n_re, ncol(Gk)))
    }
    np <- c(np, n_re)
  }

  q <- ncol(Z)
  n_tau <- 2L + nk + as.integer(n_re > 0L)
  Ginv_n <- matrix(0, q, n_tau)
  Ginv_n[, 1L] <- G1inv_n
  Ginv_n[, 2L] <- G2inv_n
  if (nk > 0L) {
    Ginv_n[, 2L + seq_len(nk)] <- Gk
  }
  if (n_re > 0L) {
    Ginv_n[q - n_re + seq_len(n_re), n_tau] <- 1
  }
  stopifnot(nrow(Ginv_n) == ncol(Z), ncol(Ginv_n) == n_tau)

  list(
    X = X, Z = Z, Z_re = Z_re, efine = efine, C = C,
    groups = groups, n_coarse = n_coarse,
    Ginv_n = Ginv_n, n_re = n_re, nk = nk,
    np = np, MM1 = MM1, MM2 = MM2,
    nl.level = nl_level_resolved
  )
}

#' Poisson CL-GAMM via TMB Laplace approximation
#'
#' Same mixed-model design as \code{\link{pois_SOP}} (anisotropic P-splines,
#' optional smooths, optional \code{re="coarse"}), with random effects
#' integrated by Laplace's method. Requires \pkg{TMB}. Prefer
#' \code{\link{clgam}(..., method = "Laplace")}.
#'
#' Quasi-Poisson is the same post-hoc scaling as SOP: the TMB objective
#' remains the Poisson Laplace likelihood (point estimates unchanged);
#' Pearson \eqn{\hat\phi} is stored in \code{fit$phi} and SEs from
#' \code{TMB::sdreport} are multiplied by \eqn{\sqrt{\hat\phi}}.
#'
#' @inheritParams pois_SOP
#' @param re \code{"none"} (default) or \code{"coarse"}
#' @param elements if \code{TRUE} (default), run \code{TMB::sdreport} for SEs
#' @param silent passed to \code{TMB::MakeADFun}
#' @param control list for \code{stats::nlminb}
#' @return A \code{"clgam"} object with \code{method = "Laplace"}
#' @export
#' @importFrom stats nlminb
pois_Laplace <- function(
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
  nl.level = NULL,
  sparse.backend = "auto",
  family = stats::poisson(),
  re = NULL,
  elements = TRUE,
  silent = TRUE,
  control = list(iter.max = 200, eval.max = 400)
) {
  fam <- .resolve_clgam_family(family)
  nl.basis <- match.arg(nl.basis)
  if (is.null(orth.smooth)) orth.smooth <- identical(nl.basis, "pspline")
  orth.smooth <- isTRUE(orth.smooth)
  re <- .clgam_resolve_re(re)
  if (identical(re, "coarse") && isTRUE(fam$quasi)) {
    warning(
      "re='coarse' and family=quasipoisson() both model overdispersion.",
      call. = FALSE
    )
  }
  start.all <- proc.time()[3]

  des <- .clgam_mixed_design(
    x1 = x1, x2 = x2, efine = efine,
    lcovfine = lcovfine, nlcovfine = nlcovfine, C = C,
    x1lim = x1lim, x2lim = x2lim,
    ndx = ndx, bdeg = bdeg, pord = pord, decom = decom,
    ndxnl = ndxnl, bdegnl = bdegnl, pordnl = pordnl,
    nl.basis = nl.basis, orth.smooth = orth.smooth, nl.level = nl.level,
    sparse.backend = sparse.backend, re = re
  )
  if (is.null(des$groups)) {
    stop(
      "pois_Laplace currently requires a nested 0-1 partition C.",
      call. = FALSE
    )
  }

  nk <- des$nk
  n_re <- des$n_re
  n_tau <- 2L + nk + as.integer(n_re > 0L)
  la0 <- as.numeric(parold)
  if (length(la0) < 2L) la0 <- c(la0, rep(1, 2L - length(la0)))
  la0 <- la0[1:2]
  if (nk > 0L) {
    if (is.null(paroldnl)) paroldnl <- rep(1, nk)
    la0 <- c(la0, as.numeric(paroldnl)[seq_len(nk)])
  }
  if (n_re > 0L) la0 <- c(la0, 1)

  .clgam_ensure_tmb_dll("clgam_pois")

  data <- list(
    y = as.numeric(y),
    X = as.matrix(des$X),
    Z = as.matrix(des$Z),
    e = as.numeric(des$efine),
    groups = as.integer(des$groups) - 1L,
    n_coarse = as.integer(des$n_coarse),
    Ginv_n = as.matrix(des$Ginv_n),
    n_re = as.integer(n_re)
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
    lower = ifelse(names(obj$par) == "log_tau2", log(1e-4), -Inf),
    upper = ifelse(names(obj$par) == "log_tau2", log(1e4), Inf)
  )

  b.fixed <- as.numeric(opt$par[names(opt$par) == "beta"])
  log_tau2 <- as.numeric(opt$par[names(opt$par) == "log_tau2"])
  la <- exp(log_tau2)
  obj$fn(opt$par)
  b.random <- as.numeric(
    obj$env$last.par.best[names(obj$env$last.par.best) == "u"]
  )
  eta_full <- as.numeric(des$X %*% b.fixed + des$Z %*% b.random)
  eta_full <- pmin(pmax(eta_full, -20), 20)
  re_hat <- NULL
  eta <- eta_full
  if (n_re > 0L) {
    re_hat <- b.random[(length(b.random) - n_re + 1L):length(b.random)]
    eta <- eta_full - as.numeric(des$Z_re %*% re_hat)
  }

  sd.eta <- rep(NA_real_, length(eta))
  if (isTRUE(elements)) {
    rep <- tryCatch(
      TMB::sdreport(obj, getJointPrecision = FALSE),
      error = function(e) {
        warning("pois_Laplace: sdreport failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(rep)) {
      sdr_sum <- summary(rep)
      eta_rows <- which(rownames(sdr_sum) == "eta_struct")
      if (length(eta_rows) == length(eta)) {
        sd.eta <- as.numeric(sdr_sum[eta_rows, "Std. Error"])
      }
    }
  } else {
    rep <- NULL
  }

  gamma <- as.numeric(des$efine * exp(eta_full))
  mu <- as.numeric(rowsum(
    gamma,
    group = factor(des$groups, levels = seq_len(des$n_coarse))
  ))
  eps <- .Machine$double.xmin
  mu_s <- pmax(mu, eps)
  y_n <- as.numeric(y)
  sat <- ifelse(y_n > 0, y_n * log(y_n) - y_n, 0)
  fitll <- y_n * log(mu_s) - mu_s
  dev <- 2 * sum(sat - fitll)

  vc_names <- c("spatial.x1", "spatial.x2")
  if (nk > 0L) {
    vc_names <- c(vc_names, .clgam_mat_colnames(nlcovfine, "s"))
  }
  if (n_re > 0L) vc_names <- c(vc_names, "re.coarse")
  if (length(vc_names) == length(la)) names(la) <- vc_names

  end.all <- proc.time()[3]
  out <- list(
    ndx = ndx, bdeg = bdeg, pord = pord,
    knots1 = des$MM1$knots, knots2 = des$MM2$knots,
    y = as.numeric(y), x1 = x1, x2 = x2, efine = des$efine,
    lcovfine = lcovfine, nlcovfine = nlcovfine,
    eta = eta, eta.full = eta_full, gamma = gamma, mu = mu,
    re = re_hat, var.comp = la, edf = NA_real_,
    niter = as.integer(opt$iterations),
    elapsed.time = end.all - start.all,
    diverged = !isTRUE(opt$convergence == 0),
    converged = isTRUE(opt$convergence == 0),
    dev = dev, b.fixed = b.fixed, b.random = b.random,
    matlist = list(X = des$X, Z = des$Z, C = des$C, Ginv_n = des$Ginv_n),
    leffects = NULL, sdleffects = NULL,
    nleffects = NULL, sdnleffects = NULL,
    nl.basis = nl.basis, orth.smooth = orth.smooth,
    nl.level = if (length(des$nl.level)) des$nl.level else nl.level,
    method = "Laplace",
    re_term = re,
    ed = NA_real_, aic = NA_real_, bic = NA_real_,
    sd.eta = sd.eta,
    sd.exp.eta = sd.eta * exp(eta),
    tmb = list(opt = opt, obj = obj, sdreport = rep, convergence = opt$convergence)
  )
  phi <- .pearson_phi(y_n, mu, ncol(des$X) + 2)
  out$phi <- phi
  if (isTRUE(fam$quasi)) {
    out <- .scale_se_phi(out, phi)
  }
  df_res <- max(length(y_n) - ncol(des$X), 1)
  out$df.residual <- df_res
  out$deviance_df <- as.numeric(dev) / df_res
  .as_clgam(out, call = match.call(), type = "spatial", family = fam$name)
}

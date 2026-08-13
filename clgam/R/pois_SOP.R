#' Poisson CL-GAMM via PIRLS + SOP
#'
#' Prefer \code{\link{clgam}} in new code. Fits a composite-link Poisson GAMM
#' with anisotropic spatial P-splines and optional fine-scale smooth covariates,
#' estimating variance components by separation of overlapping precision
#' matrices (SOP). \code{pois_SAP} is retained as a compatibility alias.
#'
#' @param y coarse counts (length \eqn{n})
#' @param x1,x2 fine-scale spatial coordinates (length \eqn{m})
#' @param efine fine-scale exposure / expected counts (length \eqn{m}; default 1)
#' @param lcovfine optional fine-scale linear covariates (\eqn{m \times p})
#' @param nlcovfine optional fine-scale covariates for univariate P-spline
#'   smooths (\eqn{m \times K})
#' @param C composition matrix (\eqn{n \times m})
#' @param x1lim,x2lim optional coordinate limits
#' @param ndx,bdeg,pord spatial P-spline settings
#' @param decom mixed-model decomposition for spatial bases
#' @param thr convergence tolerances \code{c(eta, tau2)}
#' @param maxit max outer / inner iterations
#' @param parold,paroldnl initial variance components
#' @param bold optional starting coefficients
#' @param trace print iteration banner
#' @param elements compute AIC/BIC and SEs
#' @param ndxnl,bdegnl,pordnl settings for each smooth covariate
#' @param sparse.backend composition-matrix backend (see \code{\link{as_comp_C}})
#' @param nl.basis how smooth covariates enter the mixed model:
#'   \code{"pspline"} uses the P-spline null space \code{Xxk} plus penalized
#'   \code{Zxk} from \code{mm_basis}; \code{"legacy"} (default here) matches
#'   archival SMiMR code (raw \code{nlcovfine} unpenalized in \code{X} +
#'   \code{Zxk}). Prefer \code{\link{clgam}}, which defaults to
#'   \code{nl.basis="pspline"}.
#' @param orth.smooth if \code{TRUE}, project spatial bases (except the intercept)
#'   orthogonal to the space of \emph{fine-scale} covariate effects so
#'   \eqn{f(s)} cannot absorb \eqn{g(z)}. The projection uses the raw fine
#'   B-spline bases \eqn{B(z)} (not their mixed-model reparameterization) and,
#'   when present, linear fine covariates \code{lcovfine}. Coarse-scale smooths
#'   marked via \code{nl.level="coarse"} are excluded. Defaults to
#'   \code{TRUE} for \code{nl.basis="pspline"}, \code{FALSE} for
#'   \code{"legacy"}.
#' @param nl.level optional length-\eqn{K} labels \code{"fine"} / \code{"coarse"}
#'   for columns of \code{nlcovfine} (recycled if length 1). \code{NULL} treats
#'   all nonlinear columns as fine (Case A). Use \code{"coarse"} for Case B and
#'   \code{c("fine","coarse")} for Case C so aggregated \eqn{h(z_a)} does not
#'   enter the spatial identifiability projection.
#' @param family \code{stats::poisson()} (default) or \code{stats::quasipoisson()}.
#'   Quasi-Poisson keeps the same PIRLS+SOP point estimates as Poisson; a
#'   Pearson dispersion \eqn{\hat\phi} is estimated at convergence and all
#'   standard errors (including effect SEs) are multiplied by
#'   \eqn{\sqrt{\hat\phi}}. Negative binomial is not supported.
#' @param re \code{NULL}/\code{"none"} (default) or \code{"coarse"}: an iid
#'   Gaussian random effect per coarse observation, so
#'   \eqn{\mu_i = [C(e\odot\exp\eta)]_i\exp(u_i)} with
#'   \eqn{u_i\sim N(0,\sigma_u^2)}. Requires a nested 0-1 partition
#'   \code{C}. The stored \code{eta} is the structured fine surface
#'   (without \eqn{u}); \code{re} holds \eqn{\hat u}. Combining
#'   \code{re="coarse"} with \code{quasipoisson()} is allowed but both
#'   target overdispersion (a warning is issued).
#' @param parold_re starting \eqn{\tau^2} for \code{re="coarse"}
#' @return A \code{"clgam"} object with fine-scale \code{eta}, fitted means,
#'   variance components, and optional SEs / AIC. Quasi-Poisson fits also
#'   store \code{phi} and \code{deviance_df}.
#' @details Inner SOP inverts the Schur complement \eqn{S} with a dense
#'   \code{solve} (Rcpp when compiled). Set
#'   \code{options(clgam.sop.backend = "kron_hybrid")} to use an exact
#'   B3-vs-rest block inverse of the same \eqn{S} (same estimator; see
#'   \code{.sop_Sinv_kron}). The default \code{"dense"} is unchanged.
#'   Contrast / stacked spatial models (\code{\link{pois_incat_SOP}}) ignore
#'   this option.
#' @export
#' @seealso \code{\link{clgam}}, \code{\link{pois_incat_SOP}},
#'   \code{\link{kappa_diagnostic}}
pois_SOP <- function(y, x1, x2, efine = NULL, lcovfine = NULL, nlcovfine = NULL, C, x1lim = NULL, x2lim = NULL, ndx = c(15, 15), bdeg = c(3, 3), pord = c(2, 2), decom = 1, thr = c(1e-06, 1e-06), maxit = c(100, 100), parold = c(1, 1), bold = NULL, trace = FALSE, elements = FALSE, ndxnl = 15, bdegnl = 3, pordnl = 2, paroldnl = NULL, sparse.backend = "auto", nl.basis = c("legacy", "pspline"), orth.smooth = NULL, nl.level = NULL, family = stats::poisson(), re = NULL, parold_re = 1) {
  fam <- .resolve_clgam_family(family)
  re <- .clgam_resolve_re(re)
  if (identical(re, "coarse") && isTRUE(fam$quasi)) {
    warning(
      "re='coarse' and family=quasipoisson() both model overdispersion.",
      call. = FALSE
    )
  }
  nl.basis <- match.arg(nl.basis)
  if (is.null(orth.smooth)) {
    orth.smooth <- identical(nl.basis, "pspline")
  }
  orth.smooth <- isTRUE(orth.smooth)
  start.all <- proc.time()[3]
  if (is.null(x1lim)) x1lim <- .clgam_xlim(x1)
  if (is.null(x2lim)) x2lim <- .clgam_xlim(x2)
  dimfine <- length(x1)
  efine <- .clgam_exposure(efine, dimfine)

  # Setup for mixed model representation (spatial part)
  MM1 <- mm_basis(x = x1, xl = x1lim[1], xr = x1lim[2], ndx = ndx[1], bdeg = bdeg[1], pord = pord[1], decom = decom)
  MM2 <- mm_basis(x = x2, xl = x2lim[1], xr = x2lim[2], ndx = ndx[2], bdeg = bdeg[2], pord = pord[2], decom = decom)

  X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; c1 <- MM1$m
  X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; c2 <- MM2$m

  # Build diagonal matrices for the inverse of G (Ginv, spatial part)
  g1u <- rep(d1, times = pord[2])
  g2u <- rep(d2, each = pord[1])
  g1b <- rep(d1, times = (c2 - pord[2]))
  g2b <- rep(d2, each = (c1 - pord[1]))

  # Setup for mixed model representation (non-linear explanatory variables)
  if (!is.null(nlcovfine)) {
    nk <- ncol(nlcovfine)
    if (length(ndxnl) == 1) {
      ndxnl <- rep(ndxnl, nk)
    }
    if (length(bdegnl) == 1) {
      bdegnl <- rep(bdegnl, nk)
    }
    if (length(pordnl) == 1) {
      pordnl <- rep(pordnl, nk)
    }
    Xxk <- vector("list", nk)
    Zxk <- vector("list", nk)
    Bk <- vector("list", nk)
    dk <- vector("list", nk)
    ck <- rep(0, nk)
    for (k in 1:nk) {
      xk <- nlcovfine[, k]
      xklim <- c(min(xk) - 0.01, max(xk) + 0.01)
      # Univariate smooths: decom=2 -> X = [1, x, ..., x^{pord-1}] (clear null space).
      # Spatial bases keep the caller's `decom`.
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

  # Number of parameters in each part
  # np <-  c(prod(pord), (c1 - pord[1])*pord[2], (c2 - pord[2])*pord[1], (c1 - pord[1])*(c2 - pord[2]))
  np <-  c(prod(pord), (c2 - pord[2])*pord[1], (c1 - pord[1])*pord[2], (c1 - pord[1])*(c2 - pord[2]))

  if (!is.null(lcovfine)) {
    np[1] <- np[1] + ncol(lcovfine)
  }

  if (!is.null(nlcovfine)) {
    if (identical(nl.basis, "pspline")) {
      # Intercept column of Xxk is dropped later (spatial already has one)
      for (k in 1:nk) {
        np[1] <- np[1] + max(ncol(Xxk[[k]]) - 1L, 0L)
      }
    } else {
      np[1] <- np[1] + ncol(nlcovfine)
    }
    for (k in 1:nk) {
      np <- c(np, (ck[k] - pordnl[k]))
    }
  }

  # 0	0
  # 0 I
  D <- diag(c(rep(0, np[1]), rep(1, sum(np[-1]))))

  # Build diagonal matrices G1inv.n and G2inv.n
  # * Note: Ginv = lam1*G1inv.n + lam2*G2inv.n
  # G1inv.n <- c(g1u, rep(0, np[3]), g1b)
  # G2inv.n <- c(rep(0, np[2]), g2u, g2b)
  G1inv.n <- c(rep(0, np[2]), g1u, g1b)
  G2inv.n <- c(g2u, rep(0, np[3]), g2b)

  if (!is.null(nlcovfine)) {
    G1inv.n <- c(G1inv.n, rep(0, sum(np[-(1:4)])))
    G2inv.n <- c(G2inv.n, rep(0, sum(np[-(1:4)])))
    Gkinv.n <- NULL
    for (k in 1:nk){
      Gkinv.n[[k]] <- rep(0, sum(np[-1]))
      low <- sum(np[2:(k+3)]) + 1
      sup <- low + np[k+4] - 1
      Gkinv.n[[k]][low:sup] <- dk[[k]]
    }
  }

  n_nl_random <- if (!is.null(nlcovfine)) {
    as.integer(sum(np[-(1:4)]))
  } else {
    0L
  }

  # Build spatial mixed model matrices
  X <- rten2(X2, X1)
  Z <- cbind(rten2(Z2, X1), rten2(X2, Z1), rten2(Z2, Z1))
  n_sp_fixed <- ncol(X)
  n_sp_random <- ncol(Z)

  # Identifiability: spatial bases perp fine-scale covariate space A_f before
  # appending covariate columns. A_f = span{B(z) for fine smooths} union {lcovfine}.
  # Coarse-scale smooths (nl.level="coarse") are excluded -- they enter after C.
  # A_f is built whenever fine covariates exist so kappa (eq. 21) is defined
  # even if orth.smooth is FALSE.
  nl_level_resolved <- character(0)
  if (!is.null(nlcovfine)) {
    nl_level_resolved <- .resolve_nl_level(nk, nl.level)
  }
  Af <- .build_orth_Af(
    Bk = if (!is.null(nlcovfine)) Bk else NULL,
    nl.level = nl_level_resolved,
    lcovfine = lcovfine
  )
  orth_info <- list(
    applied = FALSE,
    nl.level = nl_level_resolved,
    max_abs_QX = NA_real_,
    max_abs_QZ = NA_real_,
    kappa = NA_real_
  )
  if (orth.smooth && !is.null(Af) && ncol(Af) > 0L) {
    if (n_sp_fixed > 1L) {
      X[, 2:n_sp_fixed] <- .orth_cols(X[, 2:n_sp_fixed, drop = FALSE], Af)
    }
    Z[, seq_len(n_sp_random)] <- .orth_cols(
      Z[, seq_len(n_sp_random), drop = FALSE], Af
    )
    orth_info$applied <- TRUE
    orth_info$max_abs_QX <- if (n_sp_fixed > 1L) {
      .orth_cross_max(X[, 2:n_sp_fixed, drop = FALSE], Af)
    } else {
      0
    }
    orth_info$max_abs_QZ <- .orth_cross_max(
      Z[, seq_len(n_sp_random), drop = FALSE], Af
    )
  }

  # Include linear fine-scale covariates into X (after spatial orthogonalization)
  if (!is.null(lcovfine)) {
    X <- cbind(X, lcovfine)
  }

  # Include non-linear covariates into X & Z
  nl_fixed_idx <- NULL
  if (!is.null(nlcovfine)) {
    nl_fixed_idx <- vector("list", nk)
    if (identical(nl.basis, "pspline")) {
      # Null space without constant: spatial tensor already has an intercept.
      # For pord=2 this is the linear trend in z; Zxk carries the smooth deviation.
      for (k in 1:nk) {
        Xk <- Xxk[[k]]
        if (ncol(Xk) > 1L) {
          Xk <- Xk[, -1L, drop = FALSE]
        } else {
          Xk <- Xk[, FALSE, drop = FALSE]
        }
        Xxk[[k]] <- Xk
        if (ncol(Xk) > 0L) {
          j0 <- ncol(X)
          X <- cbind(X, Xk)
          nl_fixed_idx[[k]] <- (j0 + 1L):ncol(X)
        } else {
          nl_fixed_idx[[k]] <- integer(0)
        }
        Z <- cbind(Z, Zxk[[k]])
      }
    } else {
      # Archival SMiMR: raw covariate unpenalized + Zxk only
      j0 <- ncol(X)
      X <- cbind(X, nlcovfine)
      for (k in 1:nk) {
        nl_fixed_idx[[k]] <- j0 + k
        Z <- cbind(Z, Zxk[[k]])
      }
    }
  }

  # Sparse composition (Matrix / spam / dense); partition -> rowsum groups
  C <- .as_comp_C(C, backend = sparse.backend)
  C_groups <- if (.is_partition_C(C)) .partition_groups(C) else NULL

  n_re <- 0L
  Z_re <- NULL
  re_idx <- integer(0)
  Greinv.n <- NULL
  if (identical(re, "coarse")) {
    if (is.null(C_groups)) {
      stop(
        "re='coarse' requires a nested 0-1 partition C ",
        "(one fine unit per coarse cell).",
        call. = FALSE
      )
    }
    n_re <- nrow(C)
    Z_re <- .clgam_Z_re(C_groups, n_re)
    q0 <- ncol(Z)
    Z <- cbind(Z, Z_re)
    np <- c(np, n_re)
    G1inv.n <- c(G1inv.n, rep(0, n_re))
    G2inv.n <- c(G2inv.n, rep(0, n_re))
    if (!is.null(nlcovfine)) {
      for (k in seq_len(nk)) {
        Gkinv.n[[k]] <- c(Gkinv.n[[k]], rep(0, n_re))
      }
    }
    Greinv.n <- c(rep(0, q0), rep(1, n_re))
    re_idx <- q0 + seq_len(n_re)
    n_nl_random <- n_nl_random + n_re
  }

  kron.meta <- tryCatch(
    .sop_kron_meta(MM1, MM2, pord = pord, n_nl_random = n_nl_random),
    error = function(e) NULL
  )
  if (is.null(bold)) {
    bold <- rep(0, sum(np))
  }

  # Compute eta, gamma, and mu vectors
  eta <- c(X %*% bold[1:np[1]] + Z %*% bold[-(1:np[1])])
  gamma <- c(efine*exp(eta))
  mu <- if (is.null(C_groups)) {
    as.numeric(C %*% gamma)
  } else {
    as.numeric(rowsum(gamma, group = factor(C_groups, levels = seq_len(nrow(C)))))
  }

  # Set starting variance components
  la <- parold
  if (!is.null(nlcovfine)) {
    if (is.null(paroldnl)) {
      paroldnl <- rep(1, ncol(nlcovfine))
    }
    la <- c(la, paroldnl)
  }
  if (n_re > 0L) {
    la <- c(la, as.numeric(parold_re)[1L])
  }

  # Optimization procedure
  diverged <- FALSE
  pirls_elements <- NULL
  tol <- NA_real_
  ed_re <- NA_real_
  for (i in 1:(maxit[1])) {
    # Set the clock for SOP
    start.SOP <- proc.time()[3]
    # Compute working vector
    z <- .clmm_working_z(C, gamma, eta, mu, y, groups = C_groups)
      # Compute CLMM matrices
    mat <- clmm_mat(C, gamma, X, Z, z, mu, groups = C_groups)
    sop_cache <- NULL

    for (it in 1:(maxit[2])) {
      # Compute penalty matrix: block diagonal matrix
      # Ginv <- c((1/la[1])*g1u, (1/la[2])*g2u, (1/la[1])*g1b + (1/la[2])*g2b)
      Ginv <- c((1/la[2])*g2u, (1/la[1])*g1u, (1/la[2])*g2b + (1/la[1])*g1b)
      if (!is.null(nlcovfine)) {
        for (k in 1:nk) {
          Ginv <- c(Ginv, (1/la[(k+2)])*dk[[k]])
        }
      }
      if (n_re > 0L) {
        Ginv <- c(Ginv, rep(1 / la[length(la)], n_re))
      }
      G <- 1/Ginv

      sop <- .sop_solve_schur(
        XtX = mat$XtX, ZtX = mat$ZtX, ZtZ = mat$ZtZ,
        ZtXtZ = mat$ZtXtZ, u = mat$u, G = G, cache = sop_cache,
        meta = kron.meta
      )
      sop_cache <- sop$cache
      b.fixed <- sop$b.fixed
      b.random <- sop$b.random
      dZtNZ <- sop$dZtNZ

      # Tau 1. ED is floored at a small positive constant (not just guarded
      # against exact 0) so a numerically-near-zero-but-nonzero ED cannot
      # still send tau1 to Inf/NaN.
      G1inv.d <- (1/la[1])*G1inv.n
      ed1 <- max(sum(dZtNZ*(G1inv.d*G^2)), 1e-50)
      tau1 <- sum(b.random^2 * G1inv.n) / ed1
      tau1 <- max(tau1, 1e-50)

      # Tau 2
      G2inv.d <- (1 / la[2]) * G2inv.n
      ed2 <- max(sum(dZtNZ * (G2inv.d * G^2)), 1e-50)
      tau2 <- sum(b.random^2 * G2inv.n) / ed2
      tau2 <- max(tau2, 1e-50)

      # Rest of tau's corresponding to each non- linear explanatory variable
      if (!is.null(nlcovfine)) {
        tauk <- rep(0, nk)
        edk <- tauk
        for (k in 1:nk) {
          Gkinv.d <- (1 / la[(k + 2)]) * Gkinv.n[[k]]
          edk[k] <- max(sum(dZtNZ * (Gkinv.d * G^2)), 1e-50)
          tauk[k] <- sum(b.random^2 * Gkinv.n[[k]]) / edk[k]
          tauk[k] <- max(tauk[k], 1e-50)
        }
      }

      # New variance components and convergence check
      lanew <- c(tau1, tau2)
      if (!is.null(nlcovfine)) {
        lanew <- c(lanew, tauk)
      }
      if (n_re > 0L) {
        Greinv.d <- (1 / la[length(la)]) * Greinv.n
        ed_re <- max(sum(dZtNZ * (Greinv.d * G^2)), 1e-50)
        tau_re <- sum(b.random^2 * Greinv.n) / ed_re
        tau_re <- max(tau_re, 1e-50)
        lanew <- c(lanew, tau_re)
      }
      dla <- mean(abs(la - lanew))
      # Early exit when relative change is already tiny (avoids grinding to maxit=100)
      dla_rel <- mean(abs(la - lanew) / pmax(abs(la), abs(lanew), 1e-8))

      if (trace) {
        if (!is.null(nlcovfine)) {
          cat(sprintf("%1$3d %2$10.6f", it, dla))
          cat(sprintf("%8.3f", c(ed1, ed2, edk)), "\n")
        } else {
          cat(sprintf("%1$3d %2$10.6f", it, dla))
          cat(sprintf("%8.3f", c(ed1, ed2)), "\n")
        }
      }
      # Check finiteness of the CANDIDATE update before committing it to `la`.
      # Previously `la <- lanew` ran unconditionally before this check, so a
      # divergent tau (Inf/NaN) could be written into `la` and silently
      # propagate into the next outer PIRLS iteration's penalty. Now: on a
      # non-finite update we keep the last finite `la` (and the b.fixed /
      # b.random already computed with it, which are still valid), warn once,
      # and stop both loops.
      if (!all(is.finite(lanew)) || !is.finite(dla) || !is.finite(dla_rel)) {
        warning(
          "clgam: non-finite variance-component update at outer iteration ",
          i, ", inner iteration ", it,
          "; keeping the last finite estimate and stopping (fit not converged).",
          call. = FALSE
        )
        diverged <- TRUE
        break
      }
      la <- lanew
      if (dla < thr[2] || (it >= 8L && dla_rel < thr[2])) break
    }
    if (diverged) break

    # Stop the clock for SOP
    end.SOP <- proc.time()[3]

    # Print time
    if (trace) {
      cat("Elapsed time in seconds:", c(end.SOP - start.SOP),"\n")
    }

    # Hold old eta
    eta.old <- eta

    # Obtain new vectors
    eta <- c(X %*% b.fixed + Z %*% b.random)
    gamma <- c(efine*exp(eta))
    mu <- if (is.null(C_groups)) {
      as.numeric(C %*% gamma)
    } else {
      as.numeric(rowsum(gamma, group = factor(C_groups, levels = seq_len(nrow(C)))))
    }

    # Convergence criterion
    tol <- sum((eta - eta.old)^2)/sum(eta^2)

    if (trace) {cat("Convergence criterion: ", tol, "\n")}
    # Cache PIRLS-weighted covariance blocks once at the final outer step
    # (reused below for AIC/BIC/SEs when elements=TRUE).
      if (isTRUE(elements) && (tol < thr[1] || i == maxit[1])) {
      ginvsp_pe <- .ginvsp_from_la(la, g2u, g1u, g2b, g1b,
                                   dk = if (!is.null(nlcovfine)) dk else NULL,
                                   n_re = n_re)
      pirls_elements <- .pirls_bayes_blocks(
        C, gamma, eta, mu, y, X, Z, ginvsp_pe, groups = C_groups
      )
    }
    if (tol < (thr[1])) break
  }

  # Compute deviance (same formula for Poisson and quasi-Poisson)
  dev <- 2*sum(y*log(ifelse(y == 0, 1, y/mu)) - (y - mu))

  eta.full <- eta
  re_hat <- NULL
  q_struct <- ncol(Z) - n_re
  if (n_re > 0L) {
    re_hat <- as.numeric(b.random[re_idx])
    eta <- eta.full - as.numeric(Z_re %*% re_hat)
  }

  # Obtain inverse of G (store diagonal; expand only if callers need matrix)
  ginvsp <- .ginvsp_from_la(la, g2u, g1u, g2b, g1b,
                            dk = if (!is.null(nlcovfine)) dk else NULL,
                            n_re = n_re)
  # Stored only for inspection (no internal code path reads matlist$Ginv);
  # ginvsp is structurally diagonal, so a dense q x q diag() here wastes
  # O(q^2) memory (q = number of penalized coefficients, can be in the
  # hundreds+) for what is a length-q vector. Matrix::Diagonal keeps the same
  # `[i,j]`/`%*%`-compatible matrix semantics without materialising it densely.
  Ginv <- Matrix::Diagonal(x = ginvsp)

  # Effective degrees of freedom
  edf <- c(ed1, ed2)
  if (!is.null(nlcovfine)) {
    edf <- c(edf, edk)
  }
  if (n_re > 0L) {
    edf <- c(edf, ed_re)
  }
  vc_names <- c("spatial.x1", "spatial.x2")
  if (!is.null(nlcovfine)) {
    vc_names <- c(vc_names, .clgam_mat_colnames(nlcovfine, "s"))
  }
  if (n_re > 0L) {
    vc_names <- c(vc_names, "re.coarse")
  }
  if (length(vc_names) == length(la)) names(la) <- vc_names
  if (length(vc_names) == length(edf)) names(edf) <- vc_names
  # Fallback ED for Pearson phi when elements=FALSE: fixed dim + random EDs
  ed_fallback <- np[1] + sum(edf)

  # Compute linear- and non-linear-effect POINT ESTIMATES (always; cheap
  # matrix-vector products, independent of `elements`).
  if (!is.null(lcovfine)) {
    leffects <- matrix(0, dimfine, ncol(lcovfine))
    for (k in 1:ncol(lcovfine)) {
      leffects[, k] <- c(X[, (k + 4)] * b.fixed[(k + 4)])
    }
  } else {
    leffects <- NULL
  }
  if (!is.null(nlcovfine)) {
    nleffects <- matrix(0, dimfine, nk)
    for (k in 1:nk) {
      low <- sum(np[2:(k + 3)]) + 1
      sup <- low + np[k + 4] - 1
      idx <- nl_fixed_idx[[k]]
      fixed_part <- if (length(idx)) {
        as.numeric(X[, idx, drop = FALSE] %*% b.fixed[idx])
      } else {
        rep(0, dimfine)
      }
      nleffects[, k] <- fixed_part + as.numeric(Zxk[[k]] %*% b.random[low:sup])
    }
  } else {
    nleffects <- NULL
  }

  sm_names <- .clgam_mat_colnames(nlcovfine, "s")
  lin_names <- .clgam_mat_colnames(lcovfine, "x")
  if (!is.null(nleffects) && length(sm_names)) {
    colnames(nleffects) <- sm_names
  }
  if (!is.null(leffects) && length(lin_names)) {
    colnames(leffects) <- lin_names
  }

  f_sp <- .clgam_spatial_from_parts(eta, nleffects, leffects)
  if (!is.null(Af) && ncol(Af) > 0L) {
    orth_info$kappa <- .kappa_overlap(f_sp, Af)
    parts_k <- .build_orth_Af_parts(
      Bk = if (!is.null(nlcovfine)) Bk else NULL,
      nl.level = nl_level_resolved,
      lcovfine = lcovfine,
      smooth_names = sm_names,
      linear_names = lin_names
    )
    if (length(parts_k)) {
      orth_info$kappa_by <- vapply(
        parts_k, function(A) .kappa_overlap(f_sp, A), numeric(1)
      )
    }
  }

  # Standard errors. When elements=TRUE, the block below (PIRLS-weighted
  # Bayesian covariance M1$S) is strictly better and computed once for BOTH
  # linear and smooth effects; the unweighted-Gram fallback here is then
  # unnecessary. Building it anyway (as earlier versions did, independently
  # for leffects and again for nleffects) meant factorising the same dense
  # (p+q)x(p+q) matrix crossprod(Cmat)+Gmat up to twice per fit and
  # discarding the result when elements=TRUE -- computed here only when
  # actually needed, and only once.
  sdleffects <- NULL
  sdnleffects <- NULL
  if ((!is.null(lcovfine) || !is.null(nlcovfine)) && !isTRUE(elements)) {
    Cmat <- cbind(X, Z)
    Gmat <- diag(c(rep(0, np[1]), ginvsp))
    Rmat <- tryCatch(
      solve(crossprod(Cmat) + Gmat),
      error = function(e) MASS::ginv(as.matrix(crossprod(Cmat) + Gmat))
    )
    if (!is.null(lcovfine)) {
      sdleffects <- matrix(0, dimfine, ncol(lcovfine))
      for (k in 1:ncol(lcovfine)) {
        # Var(X[,j]*beta_j) = X[,j]^2 * Var(beta_j); the column selector used
        # to have only one nonzero entry, so forming the full dimfine x
        # dimfine covariance matrix (as an earlier version did, via
        # Cmat %*% diag(onesk) %*% Rmat %*% t(...)) just to read off its
        # diagonal was O(dimfine^2) for an O(dimfine) quantity.
        j <- k + 4
        sdleffects[, k] <- abs(X[, j]) * sqrt(max(Rmat[j, j], 0))
      }
    }
    if (!is.null(nlcovfine)) {
      sdnleffects <- matrix(0, dimfine, nk)
      for (k in 1:nk) {
        idx <- nl_fixed_idx[[k]]
        low <- sum(np[2:(k + 3)]) + 1
        sup <- low + np[k + 4] - 1
        sel <- c(idx, np[1] + low:sup)
        Cg <- Cmat[, sel, drop = FALSE]
        Sg <- Rmat[sel, sel, drop = FALSE]
        # Sum-to-zero / centred SE (matches centred g plots; drops level uncertainty)
        sdnleffects[, k] <- .se_smooth_centered(Cg, Sg)
      }
    }
  }

  end.all <- proc.time()[3]
  comp.time <- end.all - start.all
  .clgam_report(trace, i, la, tol, comp.time)

  out <- list(
    ndx = ndx, bdeg = bdeg, pord = pord, decom = decom,
    ndxnl = ndxnl, bdegnl = bdegnl, pordnl = pordnl,
    knots1 = MM1$knots, knots2 = MM2$knots,
    y = y, x1 = x1, x2 = x2, efine = efine,
    lcovfine = lcovfine, nlcovfine = nlcovfine,
    eta = eta, eta.full = eta.full, gamma = gamma, mu = mu,
    re = re_hat, var.comp = la, edf = edf, niter = i, elapsed.time = comp.time,
    diverged = diverged,
    tol = tol, maxit = maxit, thr = thr,
    converged = isTRUE(!diverged && is.finite(tol) && tol < thr[1]),
    dev = dev, b.fixed = b.fixed, b.random = b.random,
    matlist = list(
      B1 = MM1$B, B2 = MM2$B, D1 = MM1$D, D2 = MM2$D,
      X = X, Z = Z, C = C, Ginv = Ginv
    ),
    leffects = leffects, sdleffects = sdleffects,
    nleffects = nleffects, sdnleffects = sdnleffects,
    nl.basis = nl.basis, orth.smooth = orth.smooth,
    nl.level = if (length(nl_level_resolved)) nl_level_resolved else nl.level,
    orth.info = orth_info,
    method = "SOP",
    re_term = re
  )

  if (isTRUE(elements)) {
    if (is.null(pirls_elements)) {
      pirls_elements <- .pirls_bayes_blocks(
        C, gamma, eta, mu, y, X, Z, ginvsp, groups = C_groups
      )
    }
    opt.mat <- pirls_elements$mat
    M1 <- pirls_elements$M1
    M2 <- pirls_elements$M2
    ed <- trprod(M1$S, M2)
    aic <- dev + 2 * ed
    bic <- dev + log(length(y)) * ed
    sd.eta <- sqrt(
      .quad_diag(X, M1$S11) +
        2 * rowSums((X %*% M1$S12[, seq_len(q_struct), drop = FALSE]) *
                      Z[, seq_len(q_struct), drop = FALSE]) +
        .quad_diag(
          Z[, seq_len(q_struct), drop = FALSE],
          M1$S22[seq_len(q_struct), seq_len(q_struct), drop = FALSE]
        )
    )
    out$ed <- ed
    out$aic <- aic
    out$bic <- bic
    out$sd.eta <- sd.eta
    out$sd.exp.eta <- sd.eta * exp(eta)

    # Standard errors from the PIRLS-weighted Bayesian covariance, computed
    # once and shared by BOTH linear and smooth effects. Previously only
    # sdnleffects was upgraded here; sdleffects (linear covariates) silently
    # kept the cruder unweighted-Gram fallback SE even when elements=TRUE,
    # an inconsistency between the two effect types with no test coverage.
    if ((!is.null(lcovfine) || !is.null(nlcovfine))) {
      Vp <- rbind(
        cbind(M1$S11, M1$S12),
        cbind(M1$S21, M1$S22)
      )
      Cmat <- cbind(X, Z)
      if (!is.null(lcovfine)) {
        sdleffects <- matrix(0, dimfine, ncol(lcovfine))
        for (k in 1:ncol(lcovfine)) {
          # See the non-elements branch above: Var(X[,j]*beta_j) =
          # X[,j]^2 * Var(beta_j) directly, no dimfine x dimfine matrix needed.
          j <- k + 4
          sdleffects[, k] <- abs(X[, j]) * sqrt(max(Vp[j, j], 0))
        }
        out$sdleffects <- sdleffects
      }
      if (!is.null(nlcovfine) && nk > 0L) {
        sdnleffects <- matrix(0, dimfine, nk)
        for (k in 1:nk) {
          idx <- nl_fixed_idx[[k]]
          low <- sum(np[2:(k + 3)]) + 1
          sup <- low + np[k + 4] - 1
          sel <- c(idx, np[1] + low:sup)
          Cg <- Cmat[, sel, drop = FALSE]
          Sg <- Vp[sel, sel, drop = FALSE]
          sdnleffects[, k] <- .se_smooth_centered(Cg, Sg)
        }
        out$sdnleffects <- sdnleffects
      }
    }
  }

  # Pearson dispersion (always report; for poisson this should be near 1 under
  # a well-specified model). Quasi-Poisson: inflate SEs by sqrt(phi); point
  # estimates and SOP tau^2 are unchanged (same PIRLS weights W = diag(mu)).
  # Use [["ed"]] -- out$ed would partial-match out$edf when elements=FALSE.
  ed_for_phi <- if (!is.null(out[["ed"]])) out[["ed"]] else ed_fallback
  phi <- .pearson_phi(y, mu, ed_for_phi)
  if (!fam$quasi) {
    # Keep poisson path bit-stable for SEs; still store diagnostic phi.
    out$phi <- phi
  } else {
    out$phi <- phi
    out <- .scale_se_phi(out, phi)
  }
  df_res <- max(length(y) - as.numeric(ed_for_phi)[1L], 1)
  out$df.residual <- df_res
  out$deviance_df <- as.numeric(dev) / df_res

  .as_clgam(out, call = match.call(), type = "spatial", family = fam$name)
}

#' @rdname pois_SOP
#' @export
pois_SAP <- pois_SOP


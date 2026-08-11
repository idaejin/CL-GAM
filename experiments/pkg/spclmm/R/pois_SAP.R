pois_SAP <- function(y, x1, x2, efine = NULL, lcovfine = NULL, nlcovfine = NULL, C, x1lim = NULL, x2lim = NULL, ndx = c(15, 15), bdeg = c(3, 3), pord = c(2, 2), decom = 1, thr = c(1e-06, 1e-06), maxit = c(100, 100), parold = c(1, 1), bold = NULL, trace = FALSE, elements = FALSE, ndxnl = 15, bdegnl = 3, pordnl = 2, paroldnl = NULL, sparse.backend = "auto") {

  # Start the clock!
  start.all <- proc.time()[3]

  # Check x1lim and x2lim arguments
  if (is.null(x1lim)) {
    ## Giancarlo's approach
    #x1ran <- max(x1) - min(x1)
    #x1lim <- c(min(x1) - 0.01*x1ran, max(x1) + 0.01*x1ran)
    ## Coté's approach
    x1lim <- c(min(x1) - 0.01, max(x1) + 0.01)
    ## Paul's approach
    #x1lim <- c(min(x1), max(x1))
  }
  if (is.null(x2lim)) {
    ## Giancarlo's approach
    #x2ran <- max(x2) - min(x2)
    #x2lim <- c(min(x2) - 0.01*x2ran, max(x2) + 0.01*x2ran)
    ## Coté's approach
    x2lim <- c(min(x2) - 0.01, max(x2) + 0.01)
    ## Paul's approach
    #x2lim <- c(min(x2), max(x2))
  }

  # Check efine argument
  dimfine <- length(x1)
  if (is.null(efine)) {
    efine <- rep(1, dimfine)
  } else if (length(efine) == 1) {
    efine <- rep(efine, dimfine)
  }

  # Setup for mixed model representation (spatial part)
  MM1 <- spclmm::mm_basis(x = x1, xl = x1lim[1], xr = x1lim[2], ndx = ndx[1], bdeg = bdeg[1], pord = pord[1], decom = decom)
  MM2 <- spclmm::mm_basis(x = x2, xl = x2lim[1], xr = x2lim[2], ndx = ndx[2], bdeg = bdeg[2], pord = pord[2], decom = decom)

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
    dk <- vector("list", nk)
    ck <- rep(0, nk)
    for (k in 1:nk) {
      xk <- nlcovfine[,k]
      xklim <- c(min(xk) - 0.01, max(xk) + 0.01)
      MMk <- spclmm::mm_basis(x = xk, xl = xklim[1], xr = xklim[2], ndx = ndxnl[k], bdeg = bdegnl[k], pord = pordnl[k], decom = decom)
      Xxk[[k]] <- MMk$X
      Zxk[[k]] <- MMk$Z
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
    np[1] <- np[1] + ncol(nlcovfine)
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

  # Set starting coefficients
  if (is.null(bold)) {
    bold <- rep(0, sum(np))
  }

  # Build spatial mixed model matrices
  X <- spclmm::rten2(X2, X1)
  Z <- cbind(spclmm::rten2(Z2, X1), spclmm::rten2(X2, Z1), spclmm::rten2(Z2, Z1))

  # Include linear fine-scale covariates into X
  if (!is.null(lcovfine)) {
    X <- cbind(X, lcovfine)
  }

  # Include non-linear fine-scale covariates into X & Z
  if (!is.null(nlcovfine)) {
    X <- cbind(X, nlcovfine)
    for (k in 1:nk) {
      Z <- cbind(Z, Zxk[[k]])
    }
  }

  # Sparse composition (Matrix / spam / dense); partition → rowsum groups
  C <- .as_comp_C(C, backend = sparse.backend)
  C_groups <- if (.is_partition_C(C)) .partition_groups(C) else NULL

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

  # Optimization procedure
  for (i in 1:(maxit[1])) {
    # Set the clock for SAP
    start.SAP <- proc.time()[3]
    # Compute working vector
    z <- .clmm_working_z(C, gamma, eta, mu, y, groups = C_groups)
      # Compute CLMM matrices
    mat <- spclmm::clmm_mat(C, gamma, X, Z, z, mu, groups = C_groups)
    sap_cache <- NULL

    for (it in 1:(maxit[2])) {
      # Compute penalty matrix: block diagonal matrix
      # Ginv <- c((1/la[1])*g1u, (1/la[2])*g2u, (1/la[1])*g1b + (1/la[2])*g2b)
      Ginv <- c((1/la[2])*g2u, (1/la[1])*g1u, (1/la[2])*g2b + (1/la[1])*g1b)
      if (!is.null(nlcovfine)) {
        for (k in 1:nk) {
          Ginv <- c(Ginv, (1/la[(k+2)])*dk[[k]])
        }
      }
      G <- 1/Ginv

      sap <- .sap_solve_schur(
        XtX = mat$XtX, ZtX = mat$ZtX, ZtZ = mat$ZtZ,
        ZtXtZ = mat$ZtXtZ, u = mat$u, G = G, cache = sap_cache
      )
      sap_cache <- sap$cache
      b.fixed <- sap$b.fixed
      b.random <- sap$b.random
      dZtNZ <- sap$dZtNZ

      # Tau 1
      G1inv.d <- (1/la[1])*G1inv.n
      ed1 <- sum(dZtNZ*(G1inv.d*G^2))
      ed1 <- ifelse(ed1 == 0, 1e-50,ed1)
      tau1 <- sum(b.random^2*G1inv.n)/ed1
      tau1 <- ifelse(tau1 == 0, 1e-50, tau1)

      # Tau 2
      G2inv.d <- (1/la[2])*G2inv.n
      ed2 <- sum(dZtNZ*(G2inv.d*G^2))
      ed2 <- ifelse(ed2 == 0, 1e-50, ed2)
      tau2 <- sum(b.random^2*G2inv.n)/ed2
      tau2 <- ifelse(tau2 == 0, 1e-50, tau2)

      # Rest of tau's corresponding to each non- linear explanatory variable
      if (!is.null(nlcovfine)) {
        tauk <- rep(0, nk)
        edk <- tauk
        for (k in 1:nk) {
          Gkinv.d <- (1/la[(k+2)])*Gkinv.n[[k]]
          edk[k] <- sum(dZtNZ*(Gkinv.d*G^2))
          edk[k] <- ifelse(edk[k] == 0, 1e-50, edk[k])
          tauk[k] <- sum(b.random^2*Gkinv.n[[k]])/edk[k]
          tauk[k] <- ifelse(tauk[k] == 0, 1e-50, tauk[k])
        }
      }

      # New variance components and convergence check
      lanew <- c(tau1, tau2)
      if (!is.null(nlcovfine)) {
        lanew <- c(lanew, tauk)
      }
      dla <- mean(abs(la - lanew))
      # Early exit when relative change is already tiny (avoids grinding to maxit=100)
      dla_rel <- mean(abs(la - lanew) / pmax(abs(la), abs(lanew), 1e-8))
      la <- lanew

      if (trace) {
        if (!is.null(nlcovfine)) {
          cat(sprintf("%1$3d %2$10.6f", it, dla))
          cat(sprintf("%8.3f", c(ed1, ed2, edk)), "\n")
        } else {
          cat(sprintf("%1$3d %2$10.6f", it, dla))
          cat(sprintf("%8.3f", c(ed1, ed2)), "\n")
        }
      }
      if (dla < thr[2] || (it >= 8L && dla_rel < thr[2])) break
    }

    # Stop the clock for SAP
    end.SAP <- proc.time()[3]

    # Print time
    if (trace) {
      cat("Elapsed time in seconds:", c(end.SAP - start.SAP),"\n")
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
    if (tol < (thr[1])) break
  }

  # Compute deviance
  dev <- 2*sum(y*log(ifelse(y == 0, 1, y/mu)) - (y - mu))

  # Obtain inverse of G (store diagonal; expand only if callers need matrix)
  ginvsp <- c((1/la[2])*g2u, (1/la[1])*g1u, (1/la[2])*g2b + (1/la[1])*g1b)
  if (!is.null(nlcovfine)) {
    for (k in 1:nk) {
      ginvsp <- c(ginvsp, (1/la[(k+2)])*dk[[k]])
    }
  }
  Ginv <- diag(ginvsp, nrow = length(ginvsp))

  # Effective degrees of freedom
  edf <- c(ed1, ed2)
  if (!is.null(nlcovfine)) {
    edf <- c(edf, edk)
  }

  # Compute linear effects and standard deviations
  if (!is.null(lcovfine)) {
    # Linear effects
    leffects <- matrix(0, dimfine, ncol(lcovfine))
    for (k in 1:ncol(lcovfine)){
      leffects[,k] <- c(X[,(k+4)]*b.fixed[(k+4)])
    }
    # Associated standard deviations (ginv if Gram is rank-deficient)
    Cmat <- cbind(X, Z)
    Gmat <- diag(c(rep(0, np[1]), ginvsp))
    Rmat <- tryCatch(
      solve(crossprod(Cmat) + Gmat),
      error = function(e) MASS::ginv(as.matrix(crossprod(Cmat) + Gmat))
    )
    sdleffects <- matrix(0, dimfine, ncol(lcovfine))
    for (k in 1:ncol(lcovfine)) {
      onesk <- rep(0, ncol(Cmat))
      # Ones of the fixed part
      onesk[(k+4)] <- 1
      Ckmat <- Cmat %*% diag(onesk)
      covkmat <- Ckmat %*% Rmat %*% t(Ckmat)
      sdleffects[, k] <- sqrt(pmax(diag(covkmat), 0))
    }
  } else {
    leffects <- NULL
    sdleffects <- NULL
  }

  # Compute non-linear effects and standard deviations
  if (!is.null(nlcovfine)) {
    # Non-linear effects
    nleffects <- matrix(0, dimfine, nk)
    for (k in 1:nk) {
      low <- sum(np[2:(k+3)]) + 1
      sup <- low + np[k+4] - 1
      if (is.null(lcovfine)) {
        nleffects[,k] <- c(X[,(k+4)]*b.fixed[(k+4)]) + c(Zxk[[k]] %*% b.random[low:sup])
      } else {
        nleffects[,k] <- c(X[,(ncol(lcovfine)+k+4)]*b.fixed[(ncol(lcovfine)+k+4)]) + c(Zxk[[k]] %*% b.random[low:sup])
      }
    }
    # Associated standard deviations (ginv if Gram is rank-deficient)
    Cmat <- cbind(X, Z)
    Gmat <- diag(c(rep(0, np[1]), ginvsp))
    Rmat <- tryCatch(
      solve(crossprod(Cmat) + Gmat),
      error = function(e) MASS::ginv(as.matrix(crossprod(Cmat) + Gmat))
    )
    sdnleffects <- matrix(0, dimfine, nk)
    for (k in 1:nk) {
      onesk <- rep(0, ncol(Cmat))
      # Ones of the random part
      #low <- sum(np[2:(k+3)]) + 1
      low <- sum(np[1:(k+3)]) + 1
      sup <- low + np[k+4] - 1
      onesk[low:sup] <- 1
      # Ones of the fixed part
      if (is.null(lcovfine)) {
        onesk[(k+4)] <- 1
      } else {
        onesk[(ncol(lcovfine)+k+4)] <- 1
      }
      Ckmat <- Cmat %*% diag(onesk)
      covkmat <- Ckmat %*% Rmat %*% t(Ckmat)
      sdnleffects[, k] <- sqrt(pmax(diag(covkmat), 0))
    }
  } else {
    nleffects <- NULL
    sdnleffects <- NULL
  }

  # Stop the clock
  end.all <- proc.time()[3]
  comp.time <- end.all - start.all

  # Print a procedure summary
  cat("Number of iterations:", i, "\n")
  cat("Optimal variance components:", la, "\n")
  cat("Convergence criterion value: ", tol, "\n")
  cat("Elapsed time of estimation procedure:", comp.time, "seconds\n")

  # Output 1
  list1 <- list(ndx = ndx, bdeg = bdeg, pord = pord, knots1 = MM1$knots, knots2 = MM2$knots, eta = eta, gamma = gamma, mu = mu, var.comp = la, edf = edf, niter = i, elapsed.time = comp.time, dev = dev, b.fixed = b.fixed, b.random = b.random, matlist = list(B1 = MM1$B, B2 = MM2$B, D1 = MM1$D, D2 = MM2$D, X = X, Z = Z, C = C, Ginv = Ginv), leffects = leffects, sdleffects = sdleffects, nleffects = nleffects, sdnleffects = sdnleffects)

  if (elements) {
    # Compute updated working vector
    z <- .clmm_working_z(C, gamma, eta, mu, y, groups = C_groups)

    # Compute effective dimension
    # * Note: ed = sum(edf)
    opt.mat <- spclmm::clmm_mat(C, gamma, X, Z, z, mu, groups = C_groups)
    ZtZpen <- opt.mat$ZtZ
    diag(ZtZpen) <- diag(ZtZpen) + ginvsp
    M1 <- spclmm::inv_bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, ZtZpen)
    M2 <- spclmm::bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, opt.mat$ZtZ)
    ed <- spclmm::trprod(M1$S, M2)

    # Compute Akaike Information Criterion
    aic <- dev + 2*ed

    # Compute Bayesian Information Criterion
    bic <- dev + log(length(y))*ed

    # Approximate SEs for eta without forming n×n covariances
    sd.eta <- sqrt(
      .quad_diag(X, M1$S11) +
        2 * rowSums((X %*% M1$S12) * Z) +
        .quad_diag(Z, M1$S22)
    )

    # Compute approximate standard errors for exp(eta)
    sd.exp.eta <- sd.eta*exp(eta)

    # Output 2
    list2 <- list(ed = ed, aic = aic, bic = bic, sd.eta = sd.eta, sd.exp.eta = sd.exp.eta)

    # Joint outputs 1 & 2
    list3 <- c(list1, list2)

  }

  ifelse(elements == FALSE, return(list1), return(list3))

}

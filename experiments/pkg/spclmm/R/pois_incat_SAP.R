pois_incat_SAP <- function(y, x1, x2, efine = NULL, lcovfine = NULL, cat, Ccat1, Ccat2, x1lim = NULL, x2lim = NULL, ndx = c(15, 15), bdeg = c(3, 3), pord = c(2, 2), thr = c(1e-06, 1e-06), maxit = c(100, 100), parold = c(1, 1, 1, 1), bold = NULL, trace = FALSE, elements = FALSE, decom = 1, sparse.backend = "auto") {

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

  # Build extended composition matrix (keep sparse — do not densify)
  C <- .as_comp_C(Matrix::bdiag(Ccat1, Ccat2), backend = sparse.backend)
  C_groups <- if (.is_partition_C(C)) .partition_groups(C) else NULL

  # Setup for mixed model representation
  MM1 <- spclmm::mm_basis(x = x1, xl = x1lim[1], xr = x1lim[2], ndx = ndx[1], bdeg = bdeg[1], pord = pord[1], decom = decom)
  MM2 <- spclmm::mm_basis(x = x2, xl = x2lim[1], xr = x2lim[2], ndx = ndx[2], bdeg = bdeg[2], pord = pord[2], decom = decom)

  X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; c1 <- MM1$m
  X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; c2 <- MM2$m

  # Build diagonal matrices for the inverse of G (Ginv)
  g1u <- rep(d1, times = pord[2])
  g2u <- rep(d2, each = pord[1])
  g1b <- rep(d1, times = (c2 - pord[2]))
  g2b <- rep(d2, each = (c1 - pord[1]))

  # Number of parameters in each part
  np.aux <-  c((c2 - pord[2])*pord[1], (c1 - pord[1])*pord[2], (c1 - pord[1])*(c2 - pord[2]))
  np <- c(4 + 4, np.aux, np.aux)

  if (!is.null(lcovfine)) {
    np[1] <- np[1] + ncol(lcovfine)
  }

  # Vector of one up to the number of parameters per category
  np.cat <- 1:sum(np.aux)

  # 0	0
  # 0 I
  D <- diag(c(rep(0, np[1]), rep(1, sum(np[-1]))))

  # Build diagonal matrices G1inv.n, G2inv.n, G3inv.n and G4inv.n
  # *Note: Ginv = lam1*G1inv.n + lam2*G2inv.n + lam3*G3inv.n + lam4*G2inv.n
  G1inv.n <- c(rep(0, np[2]), g1u, g1b, rep(0, sum(np[5:7])))
  G2inv.n <- c(g2u, rep(0, np[3]), g2b, rep(0, sum(np[5:7])))
  G3inv.n <- c(rep(0, sum(np[2:4])), rep(0, np[5]), g1u, g1b)
  G4inv.n <- c(rep(0, sum(np[2:4])), g2u, rep(0, np[6]), g2b)

  # Set starting coefficients
  if (is.null(bold)) {
    bold <- rep(0, sum(np))
  }

  # Build spatial mixed model matrices
  X.aux <- spclmm::rten2(X2, X1)
  Z.aux <- cbind(spclmm::rten2(Z2, X1), spclmm::rten2(X2, Z1), spclmm::rten2(Z2, Z1))

  # Build mixed model matrices
  X <- model.matrix(rep(1, dimfine)~cat*X.aux[,-1])
  Z <- model.matrix(rep(1, dimfine)~(Z.aux):cat-1)

  # Include fine-scale covariates linearly
  if (!is.null(lcovfine)) {
    X <- cbind(X, lcovfine)
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
      Ginv <- c((1/la[2])*g2u, (1/la[1])*g1u, (1/la[2])*g2b + (1/la[1])*g1b, (1/la[4])*g2u, (1/la[3])*g1u, (1/la[4])*g2b + (1/la[3])*g1b)
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
      ed1 <- sum(dZtNZ*(G1inv.d*(G^2)))
      ed1 <- ifelse(ed1 == 0, 1e-50, ed1)
      tau1 <- sum((b.random)^2*G1inv.n)/ed1
      tau1 <- ifelse(tau1 == 0, 1e-50, tau1)

      # Tau 2
      G2inv.d <- (1/la[2])*G2inv.n
      ed2 <- sum(dZtNZ*(G2inv.d*(G^2)))
      ed2 <- ifelse(ed2 == 0, 1e-50, ed2)
      tau2 <- sum((b.random)^2*G2inv.n)/ed2
      tau2 <- ifelse(tau2 == 0, 1e-50, tau2)

      # Tau 3
      G3inv.d <- (1/la[3])*G3inv.n
      ed3 <- sum(dZtNZ*(G3inv.d*(G^2)))
      ed3 <- ifelse(ed3 == 0, 1e-50, ed3)
      tau3 <- sum((b.random)^2*G3inv.n)/ed3
      tau3 <- ifelse(tau3 == 0, 1e-50, tau3)

      # Tau 4
      G4inv.d <- (1/la[4])*G4inv.n
      ed4 <- sum(dZtNZ*(G4inv.d*(G^2)))
      ed4 <- ifelse(ed4 == 0, 1e-50, ed4)
      tau4 <- sum((b.random)^2*G4inv.n)/ed4
      tau4 <- ifelse(tau4 == 0, 1e-50, tau4)

      # New variance components and convergence check
      lanew <- c(tau1, tau2, tau3, tau4)
      dla <- mean(abs(la - lanew))
      # Early exit when relative change is already tiny (avoids grinding to maxit=100)
      dla_rel <- mean(abs(la - lanew) / pmax(abs(la), abs(lanew), 1e-8))
      la <- lanew

      if (trace){
        cat(sprintf("%1$3d %2$10.6f", it, dla))
        cat(sprintf("%8.3f", c(ed1, ed2, ed3, ed4)), "\n")
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

  # Compute deviances
  dev.cat1 <- 2*sum(y[1:nrow(Ccat1)]*log(ifelse(y[1:nrow(Ccat1)] == 0, 1, y[1:nrow(Ccat1)]/mu[1:nrow(Ccat1)])) - (y[1:nrow(Ccat1)] - mu[1:nrow(Ccat1)]))

  dev.cat2 <- 2*sum(y[-c(1:nrow(Ccat1))]*log(ifelse(y[-c(1:nrow(Ccat1))] == 0, 1, y[-c(1:nrow(Ccat1))]/mu[-c(1:nrow(Ccat1))])) - (y[-c(1:nrow(Ccat1))] - mu[-c(1:nrow(Ccat1))]))

  # Obtain inverse of G
  ginvsp <- c((1/la[2])*g2u, (1/la[1])*g1u, (1/la[2])*g2b + (1/la[1])*g1b, (1/la[4])*g2u, (1/la[3])*g1u, (1/la[4])*g2b + (1/la[3])*g1b)
  Ginv <- diag(ginvsp, nrow = length(ginvsp))

  # Stop the clock
  end.all <- proc.time()[3]
  comp.time <- end.all - start.all

  # Print a procedure summary
  cat("Number of iterations:", i, "\n")
  cat("Optimal variance components:", la, "\n")
  cat("Convergence criterion value:", tol, "\n")
  cat("Elapsed time of estimation procedure:", comp.time, "seconds\n")

  # Output 1
  list1 <- list(ndx = ndx, bdeg = bdeg, pord = pord, knots1 = MM1$knots, knots2 = MM2$knots, eta = eta, gamma = gamma, mu = mu, var.comp = la, edf = c(ed1, ed2, ed3, ed4), niter = i, elapsed.time = comp.time, dev = c(dev.cat1, dev.cat2), b.fixed = b.fixed, b.random = b.random, matlist = list(B1 = MM1$B, B2 = MM2$B, D1 = MM1$D, D2 = MM2$D, X = X, Z = Z, C = C, Ginv = Ginv))

  if (elements) {
    # Compute updated working vector
    z <- .clmm_working_z(C, gamma, eta, mu, y, groups = C_groups)

    # Compute effective dimension
    opt.mat <- spclmm::clmm_mat(C, gamma, X, Z, z, mu, groups = C_groups)
    ZtZpen <- opt.mat$ZtZ
    diag(ZtZpen) <- diag(ZtZpen) + ginvsp
    M1 <- spclmm::inv_bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, ZtZpen)
    M2 <- spclmm::bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, opt.mat$ZtZ)
    ed <- spclmm::trprod(M1$S, M2)

    # Compute Akaike Information Criterion
    aic1 <- dev.cat1 + 2*ed
    aic2 <- dev.cat2 + 2*ed

    # Compute Bayesian Information Criterion
    bic1 <- dev.cat1 + log(length(y[c(1:nrow(Ccat1))]))*ed
    bic2 <- dev.cat2 + log(length(y[-c(1:nrow(Ccat1))]))*ed

    # Approximate SEs for eta (no n×n Gram matrix)
    sd.eta <- sqrt(
      .quad_diag(X, M1$S11) +
        2 * rowSums((X %*% M1$S12) * Z) +
        .quad_diag(Z, M1$S22)
    )

    # Difference surface SEs: Var(eta1 - eta2) via contrast rows of [X Z]
    n2 <- dimfine / 2
    B <- cbind(X, Z)
    Dmat <- B[seq_len(n2), , drop = FALSE] - B[(n2 + 1):dimfine, , drop = FALSE]
    Sfull <- M1$S
    sd.dif <- sqrt(.quad_diag(Dmat, Sfull))

    # Another approach (independence approx — kept for API compatibility)
    sd.dif2 <- sqrt((sd.eta[1:n2])^2 + (sd.eta[-(1:n2)])^2)

    # Output 2
    list2 <- list(ed = ed, aic = c(aic1, aic2), bic = c(bic1, bic2), sd.eta = sd.eta, sd.dif = sd.dif, sd.dif2 = sd.dif2)

    # Joint outputs 1 & 2
    list3 <- c(list1, list2)

  }

  ifelse(elements == FALSE, return(list1), return(list3))

}

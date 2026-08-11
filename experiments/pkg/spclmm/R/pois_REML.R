pois_REML <- function(y, x1, x2, efine = NULL, lcovfine = NULL, nlcovfine = NULL, C, x1lim = NULL, x2lim = NULL, ndx = c(15, 15), bdeg = c(3, 3), pord = c(2, 2), decom = 1, thr = 1e-06, maxit = 100, parold = c(1, 1), etaold = NULL, trace = FALSE, elements = FALSE, ndxnl = 15, bdegnl = 3, pordnl = 2, paroldnl = NULL) {

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

  # Check etaold argument
  if (is.null(etaold)) {
    Cs <- C*(1/rowSums(C))
    ynfine <- c(crossprod(Cs, y))
    eta <- log((ynfine + 1e-06)/(efine + 1e-06))
  } else {
    eta <- etaold
  }

  # Setup for mixed model representation (spatial component)
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

  # Compute gamma and mu vectors
  gamma <- c(efine*exp(eta))
  mu <- c(C %*% gamma)

  # Set starting smoothing parameters
  lam <- parold
  if (!is.null(nlcovfine)) {
    if (is.null(paroldnl)) {
      paroldnl <- rep(1, ncol(nlcovfine))
    }
    lam <- c(lam, paroldnl)
  }

  # Optimization procedure
  for (i in 1:maxit) {
    # Compute working vector
    z <- t(t((1/mu)*C)*gamma) %*% eta + (y - mu)/mu
    # Compute CLMM matrices
    mat <- spclmm::clmm_mat(C, gamma, X, Z, z, mu)
    # Assign required elements
    aI <- - sum(log(mu))
    ZtWZ <- mat$ZtZ
    XtWX <- mat$XtX
    ZtWX <- mat$ZtX
    ztWz <- mat$yty
    ZtWz <- mat$Zty
    XtWz <- mat$Xty

    REML <- function(lam) {
      # Change parameters
      lam <- mgcv::notExp2(lam)
      # Compute the inverse of G
      Ginv <- c(lam[2]*g2u, lam[1]*g1u, lam[2]*g2b + lam[1]*g1b)
      if (!is.null(nlcovfine)) {
        for (k in 1:nk) {
          Ginv <- c(Ginv, lam[(k+2)]*dk[[k]])
        }
      }
      G <- 1/Ginv

      H <- ZtWZ + diag(Ginv)
      Hinv <- chol2inv(chol(H))

      ### Part I: log|V|
      partI <- aI - sum(log(Ginv)) + spclmm::logdetsym(H)

      ### Part II: log|X'V^-1 X|
      b0 <- crossprod(Hinv, ZtWX)
      bI <- XtWX - crossprod(ZtWX, b0)
      partII <- spclmm::logdetsym(bI)

      ### Part III: z'(V^-1 - V^-1 X(X'V^-1X)^-1X'V^-1)z
      # -> z'*V^-1*z
      c0 <- crossprod(Hinv, ZtWz)
      cI <- ztWz - as.numeric(crossprod(ZtWz, c0))
      # -> X'V^-1z
      cII <- XtWz - as.vector(crossprod(ZtWX, c0))
      # -> (X' V^-1 X)^-1
      cIII <- chol2inv(chol(bI))
      partIII <- cI - as.numeric(crossprod(cII, crossprod(cIII, cII)))

      # Log-likelihood function
      # Note: optim() minimizes a quantity; therefore we use + into the function
      loglik <- 0.5*(partI + partII + partIII)

      # BETA
      b.fixed <- c(cIII %*% cII)

      # ALPHA
      apI <- ZtWz - crossprod(ZtWZ, c0)
      apII <- ZtWX - crossprod(ZtWZ, b0)
      b.random <- c(G*(apI - apII %*% b.fixed))

      attr(loglik, "b.fixed") <- b.fixed
      attr(loglik, "b.random") <- b.random

      return(loglik)
    }

    # Set the clock for optim function
    start.optim <- proc.time()[3]

    optimREML <- optim(par = mgcv::notLog2(lam), fn = REML, gr = NULL, method = "L-BFGS-B")

    # Stop the clock for optim function
    end.optim <- proc.time()[3]

    # Smoothing parameters at current iteration
    lam <- mgcv::notExp2(optimREML$par)

    # Value of loglik at current iteration
    llik <- - optimREML$value

    # Print summary
    if (trace) {
      cat("Iteration:", i, "\n")
      cat("Smoothing parameters (lon/lat):", lam[1:2], "\n")
      if (!is.null(nlcovfine)) {
        cat("Smoothing parameters (nlcovfine):", lam[-(1:2)], "\n")
      }
      cat("Value of approximate REML:", llik, "\n")
      cat("Elapsed time in seconds:", c(end.optim - start.optim),"\n")
    }

    # Obtain fixed and random coefficients at current iteration
    b.fixed <- attr(REML(optimREML$par), "b.fixed")
    b.random <- attr(REML(optimREML$par), "b.random")

    # Hold old eta
    eta.old <- eta

    # Obtain new vectors
    eta <- c(X %*% b.fixed + Z %*% b.random)
    gamma <- c(efine*exp(eta))
    mu <- c(C %*% gamma)

    # Convergence criterion
    tol <- sum((eta - eta.old)^2)/sum(eta^2)

    if (trace) {cat("Convergence criterion: ", tol, "\n")}
    if (tol < thr) break
  }

  # Compute deviance
  dev <- 2*sum(y*log(ifelse(y == 0, 1, y/mu)) - y + mu)

  # Obtain inverse of G
  ginvsp <- c(lam[2]*g2u, lam[1]*g1u, lam[2]*g2b + lam[1]*g1b)
  if (!is.null(nlcovfine)) {
    for (k in 1:nk) {
      ginvsp <- c(ginvsp, lam[(k+2)]*dk[[k]])
    }
  }
  Ginv <- diag(ginvsp)

  # Stop the clock
  end.all <- proc.time()[3]
  comp.time <- end.all - start.all

  # Print a procedure summary
  cat("Number of iterations:", i, "\n")
  cat("Optimal smoothing parameters:", lam, "\n")
  cat("Convergence criterion value:", tol, "\n")
  cat("Elapsed time of estimation procedure:", comp.time, "seconds\n")

  # Output 1
  list1 <- list(ndx = ndx, bdeg = bdeg, pord = pord, knots1 = MM1$knots, knots2 = MM2$knots, eta = eta, gamma = gamma, mu = mu, lam = lam, niter = i, elapsed.time = comp.time, dev = dev, b.fixed = b.fixed, b.random = b.random, matlist = list(B1 = MM1$B, B2 = MM2$B, D1 = MM1$D, D2 = MM2$D, X = X, Z = Z, C = C, Ginv = Ginv))

  if (elements) {
    # Compute updated working vector
    z <- t(t((1/mu)*C)*gamma) %*% eta + (y - mu)/mu

    # Compute effective dimension
    opt.mat <- spclmm::clmm_mat(C, gamma, X, Z, z, mu)
    M1 <- spclmm::inv_bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, opt.mat$ZtZ + Ginv)
    M2 <- spclmm::bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, opt.mat$ZtZ)
    ed <- spclmm::trprod(M1$S, M2)

    # Compute Akaike information criterion
    aic <- dev + 2*ed

    # Compute Bayesian information criterion
    bic <- dev + log(length(y))*ed

    # Compute approximate standard errors for eta
    sd.eta <- sqrt(diag(X %*% tcrossprod(M1$S11, X)) + 2*diag(X %*% tcrossprod(M1$S12, Z)) + diag(Z %*% tcrossprod(M1$S22, Z)))

    # Compute approximate standard errors for exp(eta)
    sd.exp.eta <- sd.eta*exp(eta)

    # Output 2
    list2 <- list(ed = ed, aic = aic, bic = bic, sd.eta = sd.eta, sd.exp.eta = sd.exp.eta)

    # Joint outputs 1 & 2
    list3 <- c(list1, list2)

  }

  ifelse(elements == FALSE, return(list1), return(list3))

}

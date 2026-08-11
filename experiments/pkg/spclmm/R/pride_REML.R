pride_REML <- function(y, x1, x2, efine = NULL, lcovfine = NULL, C, x1lim = NULL, x2lim = NULL, ndx = c(15, 15), bdeg = c(3, 3), pord = c(2, 2), decom = 1, thr = 1e-06, maxit = 100, parold = c(1, 1, 1), etaold = NULL, deltaold = NULL, trace = FALSE, elements = FALSE) {

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

  # Check deltaold argument
  if (is.null(deltaold)) {
    b.pride <- rep(0, length(y))
  } else {
    b.pride <- deltaold
  }

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

  # Build spatial mixed model matrices
  X <- spclmm::rten2(X2, X1)
  Z <- cbind(spclmm::rten2(Z2, X1), spclmm::rten2(X2, Z1), spclmm::rten2(Z2, Z1))

  # Include fine-scale covariates linearly
  if (!is.null(lcovfine)) {
    X <- cbind(X, lcovfine)
  }

  # Compute gamma, phi, and mu vectors
  gamma <- c(efine*exp(eta))
  phi <- c(C %*% gamma)
  mu <- exp(log(phi) + b.pride)

  # Set starting parameters
  lam <- parold

  # Optimization procedure
  for (i in 1:maxit) {
    # Compute working vector
    z <- t(t((1/phi)*C)*gamma) %*% eta + (y - mu)/mu + b.pride
    # Compute Xdip and Zdip matrices
    Xdip <- (1/phi)*(C %*% (gamma*X))
    Zdip <- (1/phi)*(C %*% (gamma*Z))

    REML <- function(par) {
      # Change parameters
      lam <- mgcv::notExp2(par)
      # Compute the inverse of G
      Ginv <- c(lam[2]*g2u, lam[1]*g1u, lam[2]*g2b + lam[1]*g1b)
      G <- 1/Ginv
      # Compute vector w
      w <- c(lam[3]*(1/(mu + lam[3]))*mu)
      # Assign required elements
      aI <- - sum(log(w))
      ZtWZ <- crossprod(Zdip, w*Zdip)
      XtWX <- crossprod(Xdip, w*Xdip)
      ZtWX <- crossprod(Zdip, w*Xdip)
      ztWz <- sum((z^2)*w)
      ZtWz <- c(crossprod(Zdip, w*z))
      XtWz <- c(crossprod(Xdip, w*z))

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

      # DELTA
      b.pride <- c((mu/(mu + lam[3]))*(z - Xdip %*% b.fixed - Zdip %*% b.random))

      attr(loglik, "b.fixed") <- b.fixed
      attr(loglik, "b.random") <- b.random
      attr(loglik, "b.pride") <- b.pride

      return(loglik)

    }

    # Set the clock for optim function
    start.optim <- proc.time()[3]

    optimREML <- optim(par = mgcv::notLog2(lam), fn = REML, gr = NULL, method = "L-BFGS-B")

    # Stop the clock for optim function
    end.optim <- proc.time()[3]

    # Parameters at current iteration
    lam <- mgcv::notExp2(optimREML$par)

    # Value of loglik at current iteration
    llik <- - optimREML$value

    # Print summary
    if (trace) {
      cat("Iteration:", i, "\n")
      cat("Smoothing parameters (lon/lat):", lam[1:2], "\n")
      cat("Kappa parameter:", lam[3], "\n")
      cat("Value of approximate REML:", llik, "\n")
      cat("Elapsed time in seconds:", c(end.optim - start.optim),"\n")
    }

    # Obtain fixed and random coefficients at current iteration
    b.fixed <- attr(REML(optimREML$par), "b.fixed")
    b.random <- attr(REML(optimREML$par), "b.random")

    # Obtain individual effects coefficients at current iteration
    b.pride <- attr(REML(optimREML$par), "b.pride")

    # Hold old eta
    eta.old <- eta

    # Obtain new vectors
    eta <- c(X %*% b.fixed + Z %*% b.random)
    gamma <- c(efine*exp(eta))
    phi <- c(C %*% gamma)
    mu <- exp(log(phi) + b.pride)

    # Convergence criterion
    tol <- sum((eta - eta.old)^2)/sum(eta^2)

    if (trace) {cat("Convergence criterion: ", tol, "\n")}
    if (tol < thr) break
  }

  # Compute deviance
  dev <- 2*sum(y*log(ifelse(y == 0, 1, y/mu)) - (y - mu))

  # Obtain inverse of G
  Ginv <- diag(c(lam[2]*g2u, lam[1]*g1u, lam[2]*g2b + lam[1]*g1b))

  # Stop the clock
  end.all <- proc.time()[3]
  comp.time <- end.all - start.all

  # Print a procedure summary
  cat("Number of iterations:", i, "\n")
  cat("Optimal smoothing parameters:", lam[1:2], "\n")
  cat("Optimal kappa parameter:", lam[3], "\n")
  cat("Convergence criterion value:", tol, "\n")
  cat("Elapsed time of estimation procedure:", comp.time, "seconds\n")

  # Output 1
  list1 <- list(ndx = ndx, bdeg = bdeg, pord = pord, knots1 = MM1$knots, knots2 = MM2$knots, eta = eta, gamma = gamma, phi = phi, mu = mu, lam = lam, niter = i, elapsed.time = comp.time, dev = dev, b.fixed = b.fixed, b.random = b.random, b.pride = b.pride, matlist = list(B1 = MM1$B, B2 = MM2$B, D1 = MM1$D, D2 = MM2$D, X = X, Z = Z, C = C, Ginv = Ginv))

  if (elements) {
    # Updated working vector
    z <- t(t((1/phi)*C)*gamma) %*% eta + (y - mu)/mu + b.pride

    # Compute updated Xdip and Zdip matrices
    Xdip <- (1/phi)*(C %*% (gamma*X))
    Zdip <- (1/phi)*(C %*% (gamma*Z))

    # Compute updated vector w
    w <- c(lam[3]*(1/(mu + lam[3]))*mu)

    # Compute updated elements
    WX <- mu*Xdip
    WZ <- mu*Zdip
    ZtWZ <- crossprod(Zdip, WZ)
    XtWX <- crossprod(Xdip, WX)
    ZtWX <- crossprod(Zdip, WX)

    # Compute effective dimension
    B1 <- spclmm::bblock3(XtWX, t(ZtWX), t(WX),
                  ZtWX, ZtWZ + Ginv, t(WZ),
                  WX, WZ, diag(mu + lam[3]))
    inv.B1 <- chol2inv(chol(B1))
    B2 <- spclmm::bblock3(XtWX, t(ZtWX), t(WX),
                          ZtWX, ZtWZ, t(WZ),
                          WX, WZ, diag(mu))
    ed <- spclmm::trprod(inv.B1, B2)

    # Compute Akaike Information Criterion
    aic <- dev + 2*ed

    # Compute Bayesian Information Criterion
    bic <- dev + log(length(y))*ed

    # Compute approximate standard errors for eta
    ZtwZ <- crossprod(Zdip, w*Zdip)
    XtwX <- crossprod(Xdip, w*Xdip)
    ZtwX <- crossprod(Zdip, w*Xdip)
    B3 <- spclmm::inv_bblock2(XtwX, t(ZtwX), ZtwX, ZtwZ + Ginv)
    sd.eta <- sqrt(diag(X %*% tcrossprod(B3$S11, X)) + 2*diag(X %*% tcrossprod(B3$S12, Z)) + diag(Z %*% tcrossprod(B3$S22, Z)))

    # Output 2
    list2 <- list(ed = ed, aic = aic, bic = bic, sd.eta = sd.eta)

    # Joint outputs 1 & 2
    list3 <- c(list1, list2)

  }

  ifelse(elements == FALSE, return(list1), return(list3))

}

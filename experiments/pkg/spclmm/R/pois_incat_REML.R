pois_incat_REML <- function(y, x1, x2, efine = NULL, lcovfine = NULL, cat, Ccat1, Ccat2, x1lim = NULL, x2lim = NULL, ndx = c(15, 15), bdeg = c(3, 3), pord = c(2, 2), thr = 1e-06, maxit = 100, parold = c(1, 1, 1, 1), etaold = NULL, trace = FALSE, elements = FALSE) {

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

  # Build extended composition matrix
  C <- as.matrix(Matrix::bdiag(Ccat1, Ccat2))

  # Check etaold argument
  if (is.null(etaold)) {
    Cs <- C*(1/rowSums(C))
    ynfine <- c(crossprod(Cs, y))
    eta <- log((ynfine + 1e-06)/(efine + 1e-06))
  } else {
    eta <- etaold
  }

  # Setup for mixed model representation
  MM1 <- spclmm::mm_basis(x = x1, xl = x1lim[1], xr = x1lim[2], ndx = ndx[1], bdeg = bdeg[1], pord = pord[1], decom = 1)
  MM2 <- spclmm::mm_basis(x = x2, xl = x2lim[1], xr = x2lim[2], ndx = ndx[2], bdeg = bdeg[2], pord = pord[2], decom = 1)

  X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; c1 <- MM1$m
  X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; c2 <- MM2$m

  # Build diagonal matrices for the inverse of G (Ginv)
  g1u <- rep(d1, times = pord[2])
  g2u <- rep(d2, each = pord[1])
  g1b <- rep(d1, times = (c2 - pord[2]))
  g2b <- rep(d2, each = (c1 - pord[1]))

  # Build spatial mixed model matrices
  X.aux <- spclmm::rten2(X2, X1)
  Z.aux <- cbind(spclmm::rten2(Z2, X1), spclmm::rten2(X2, Z1), spclmm::rten2(Z2, Z1))

  # Build mixed model matrices
  X <- model.matrix(rep(1, dimfine)~cat*X.aux[,-1])
  Z <- model.matrix(rep(1, dimfine)~(Z.aux):cat-1)
  #X <- cbind(rep(1, length(dimfine)), x1, x2, x1*x2, cat, cat*x1, cat*x2, cat*x1*x2)
  #Z <- cbind(Z.aux, cat*Z.aux)

  # Include fine-scale covariates linearly
  if (!is.null(lcovfine)) {
    X <- cbind(X, lcovfine)
  }

  # Compute gamma and mu vectors
  gamma <- c(efine*exp(eta))
  mu <- c(C %*% gamma)

  # Set starting smoothing parameters
  lam <- parold

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

    REML <- function(par) {
      # Change parameters
      lam <- mgcv::notExp2(par)

      # Compute the inverse of G
      Ginv <- c(lam[2]*g2u, lam[1]*g1u, lam[2]*g2b + lam[1]*g1b, lam[4]*g2u, lam[3]*g1u, lam[4]*g2b + lam[3]*g1b)
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
      cat("Smoothing parameters:", lam, "\n")
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

  # Compute deviances
  dev.cat1 <- 2*sum(y[1:nrow(Ccat1)]*log(ifelse(y[1:nrow(Ccat1)] == 0, 1, y[1:nrow(Ccat1)]/mu[1:nrow(Ccat1)])) - (y[1:nrow(Ccat1)] - mu[1:nrow(Ccat1)]))

  dev.cat2 <- 2*sum(y[-c(1:nrow(Ccat1))]*log(ifelse(y[-c(1:nrow(Ccat1))] == 0, 1, y[-c(1:nrow(Ccat1))]/mu[-c(1:nrow(Ccat1))])) - (y[-c(1:nrow(Ccat1))] - mu[-c(1:nrow(Ccat1))]))

  # Obtain inverse of G
  Ginv <- diag(c(lam[2]*g2u, lam[1]*g1u, lam[2]*g2b + lam[1]*g1b, lam[4]*g2u, lam[3]*g1u, lam[4]*g2b + lam[3]*g1b))

  # Stop the clock
  end.all <- proc.time()[3]
  comp.time <- end.all - start.all

  # Print a procedure summary
  cat("Number of iterations:", i, "\n")
  cat("Optimal smoothing parameters:", lam, "\n")
  cat("Convergence criterion value:", tol, "\n")
  cat("Elapsed time of estimation procedure:", comp.time, "seconds\n")

  # Output 1
  list1 <- list(ndx = ndx, bdeg = bdeg, pord = pord, knots1 = MM1$knots, knots2 = MM2$knots, eta = eta, gamma = gamma, mu = mu, lam = lam, niter = i, elapsed.time = comp.time, dev = c(dev.cat1, dev.cat2), b.fixed = b.fixed, b.random = b.random, matlist = list(B1 = MM1$B, B2 = MM2$B, D1 = MM1$D, D2 = MM2$D, X = X, Z = Z, C = C, Ginv = Ginv))

  if (elements) {
    # Compute updated working vector
    z <- t(t((1/mu)*C)*gamma) %*% eta + (y - mu)/mu

    # Compute effective dimension
    opt.mat <- spclmm::clmm_mat(C, gamma, X, Z, z, mu)
    M1 <- spclmm::inv_bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, opt.mat$ZtZ + Ginv)
    M2 <- spclmm::bblock2(opt.mat$XtX, opt.mat$XtZ, opt.mat$ZtX, opt.mat$ZtZ)
    ed <- spclmm::trprod(M1$S, M2)

    # Compute Akaike Information Criterion
    aic1 <- dev.cat1 + 2*ed
    aic2 <- dev.cat2 + 2*ed

    # Compute Bayesian Information Criterion
    bic1 <- dev.cat1 + log(length(y[c(1:nrow(Ccat1))]))*ed
    bic2 <- dev.cat2 + log(length(y[-c(1:nrow(Ccat1))]))*ed

    # Compute approximate standard errors for eta
    sd.eta <- sqrt(diag(X %*% tcrossprod(M1$S11, X)) + 2*diag(X %*% tcrossprod(M1$S12, Z)) + diag(Z %*% tcrossprod(M1$S22, Z)))

    # Matrix of differences
    I <- diag(dimfine/2)
    D <- cbind(I, -I) %*% cbind(X, Z)
    cov.mat.dif <- tcrossprod(tcrossprod(D, M1$S), D)
    sd.dif <- sqrt(diag(cov.mat.dif))

    # Another approach
    sd.dif2 <- sqrt((sd.eta[1:(dimfine/2)])^2 + (sd.eta[-(1:(dimfine/2))])^2)

    # Output 2
    list2 <- list(ed = ed, aic = c(aic1, aic2), bic = c(bic1, bic2), sd.eta = sd.eta, sd.dif = sd.dif, sd.dif2 = sd.dif2)

    # Joint outputs 1 & 2
    list3 <- c(list1, list2)

  }

  ifelse(elements == FALSE, return(list1), return(list3))

}

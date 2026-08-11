gaus_REML <- function(y, x1, x2, C, x1lim = NULL, x2lim = NULL, ndx = c(15, 15), bdeg = c(3, 3), pord = c(2, 2), decom = 1, parold = c(1, 1, 1)) {

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

  # Compute (normalized) Xdip and Zdip
  Xdip <- C %*% X
  Xdip <- (1/rowSums(C))*Xdip
  Zdip <- C %*% Z
  Zdip <- (1/rowSums(C))*Zdip

  # Compute fixed matrices
  XtXdip <- crossprod(Xdip)
  ZtZdip <- crossprod(Zdip)
  ZtXdip <- crossprod(Zdip, Xdip)

  n <- length(y)
  yty <- as.numeric(crossprod(y))
  Ztydip <- crossprod(Zdip, y)
  Xtydip <- crossprod(Xdip, y)

  # Set starting smoothing parameters
  lam <- parold

  REML <- function(par){

    lam <- mgcv::notExp2(par)

    Ginv <- c(lam[2]*g2u, lam[1]*g1u, lam[2]*g2b + lam[1]*g1b)
    invGplusZtZdip <- diag(Ginv) + (1/lam[3])*ZtZdip
    H <- chol2inv(chol(invGplusZtZdip))

    ###  partI:
    #   -> log|V|:
    partI <- n*log(lam[3]) - sum(log(Ginv)) + spclmm::logdetsym(invGplusZtZdip)

    ### partII:
    #   -> log| X' V^-1 X |:
    bI <- (1/lam[3])*XtXdip - ((1/lam[3])^2)*crossprod(ZtXdip, crossprod(H, ZtXdip))
    partII <- spclmm::logdetsym(bI)

    # partIII:
    #   -> y'(V^-1 - V^-1 X(X' V^-1 X)^-1X'V^-1)y
    c0 <- c(crossprod(H, Ztydip))
    # y'*V^-1*y
    cI <- (1/lam[3])*yty - ((1/lam[3])^2)*as.numeric(crossprod(Ztydip, c0))
    # X'V^-1y
    cII <- (1/lam[3])*Xtydip - ((1/lam[3])^2)*c(crossprod(ZtXdip, c0))
    cIII <- chol2inv(chol(bI))
    partIII <- cI - as.numeric(crossprod(cII, crossprod(cIII, cII)))

    # NOTE: optim() minimizes a quantity; therefore we use + in this function
    loglik <- 0.5*(partI + partII + partIII)

    # BETA
    b.fixed <- c(cIII %*% cII)

    # ALPHA
    apI <- (1/lam[3])*Ztydip - ((1/lam[3])^2)*crossprod(ZtZdip, crossprod(H, Ztydip))
    apII <- (1/lam[3])*ZtXdip - ((1/lam[3])^2)*crossprod(ZtZdip, crossprod(H, ZtXdip))
    b.random <- (1/Ginv)*c(apI - c(apII %*% b.fixed))

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
  cat("Smoothing parameters:", lam, "\n")
  cat("Value of approximate REML:", llik, "\n")
  cat("Elapsed time in seconds:", c(end.optim - start.optim),"\n")

  # Obtain fixed and random coefficients at current iteration
  b.fixed <- attr(REML(optimREML$par), "b.fixed")
  b.random <- attr(REML(optimREML$par), "b.random")

  # Obtain new vectors
  eta <- c(X %*% b.fixed + Z %*% b.random)
  mu <- c(((1/rowSums(C))*C) %*% eta)

  # Compute deviance
  dev <- sum((y - mu)^2)

  # Compute inverse of G
  Ginv <- diag(c(lam[2]*g2u, lam[1]*g1u, lam[2]*g2b + lam[1]*g1b))

  # Compute effective dimension
  M1 <- spclmm::inv_bblock2(XtXdip, t(ZtXdip), ZtXdip, ZtZdip + Ginv)
  M2 <- spclmm::bblock2(XtXdip, t(ZtXdip), ZtXdip, ZtZdip)
  ed <- spclmm::trprod(M1$S, M2)

  # Compute Akaike Information Criterion
  aic <- dev + 2*ed

  # Compute Bayesian Information Criterion
  bic <- dev + log(n)*ed

  # Approximate standard errors for eta
  sd.eta <- sqrt(diag(X %*% tcrossprod(M1$S11, X)) + 2*diag(X %*% tcrossprod(M1$S12, Z)) + diag(Z %*% tcrossprod(M1$S22, Z)))

  # Output
  output <- list(ndx = ndx, bdeg = bdeg, pord = pord, MM1 = MM1, MM2 = MM2, eta = eta, mu = mu, lam = lam, elapsed.time = end.optim - start.optim, dev = dev, b.fixed = b.fixed, b.random = b.random, dev = dev, Ginv = Ginv, ed = ed, aic = aic, bic = bic, sd.eta = sd.eta)

  return(output)

}

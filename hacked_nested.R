hacked_pez <- function (formula, data = list(), family = "gaussian", sp = NULL, 
          site = NULL, random.effects = list(), REML = TRUE, s2.init = NULL, 
          B.init = NULL, reltol = 10^-6, maxit = 500, tol.pql = 10^-6, 
          maxit.pql = 200, verbose = FALSE) 
{
  if (!is.factor(sp)) 
    stop("'sp' must be a factor")
  if (!is.factor(site)) 
    stop("'site' must be a factor")
  if (family == "gaussian") 
    z <- hacked_gau(formula = formula, data = data, 
                                 sp = sp, site = site, random.effects = random.effects, 
                                 REML = REML, s2.init = s2.init, B.init = B.init, 
                                 reltol = reltol, maxit = maxit, verbose = verbose)
  if (family == "binomial") {
    if (is.null(s2.init)) 
      s2.init <- 0.25
    z <- communityPGLMM.binary(formula = formula, data = data, 
                               sp = sp, site = site, random.effects = random.effects, 
                               REML = REML, s2.init = s2.init, B.init = B.init, 
                               reltol = reltol, maxit = maxit, tol.pql = tol.pql, 
                               maxit.pql = maxit.pql, verbose = verbose)
  }
  if (!is.element(family, c("gaussian", "binomial"))) 
    cat("\nSorry, but only binomial (binary) and gaussian options exist at this time")
  return(z)
}

hacked_gau <- function (formula, data = list(), family = "gaussian", sp = NULL, 
          site = NULL, random.effects = list(), REML = TRUE, s2.init = NULL, 
          B.init = NULL, reltol = 10^-8, maxit = 500, verbose = FALSE) 
{
  plmm.LL <- function(par, X, Y, Zt, St, nested = NULL, REML, 
                      verbose) {
    n <- dim(X)[1]
    p <- dim(X)[2]
    if (!is.null(St)) {
      q.nonNested <- dim(St)[1]
      sr <- Re(par[1:q.nonNested])
      iC <- sr[1] * St[1, ]
      if (length(sr) > 1) 
        for (i in 2:q.nonNested) {
          iC <- iC + sr[i] * St[i, ]
        }
      iC <- as(diag(iC), "dsCMatrix")
      Ut <- iC %*% Zt
      U <- t(Ut)
    }
    else {
      q.nonNested <- 0
      sr <- NULL
    }
    if (is.null(nested[[1]])) {
      q.Nested <- 0
    }
    else {
      q.Nested <- length(nested)
    }
    if (q.Nested == 0) {
      sn <- NULL
    }
    else {
      sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
    }
    if (q.Nested == 0) {
      iA <- as(diag(n), "dsCMatrix")
      Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
      Ut.iA.U <- Ut %*% U
      iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
    }
    else {
      A <- as(diag(n), "dsCMatrix")
      for (j in 1:q.Nested) {
        A <- A + sn[j]^2 * nested[[j]]
      }
      iA <- solve(A)
      if (q.nonNested > 0) {
        Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
        Ut.iA.U <- Ut %*% iA %*% U
        iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% 
          Ut %*% iA
      }
      else {
        iV <- iA
      }
    }
    denom <- t(X) %*% iV %*% X
    num <- t(X) %*% iV %*% Y
    B <- solve(denom, num)
    B <- as.matrix(B)
    H <- Y - X %*% B
    if (q.Nested == 0) {
      logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
      if (is.infinite(logdetV)) 
        logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
    }
    else {
      logdetV <- -determinant(iV)$modulus[1]
      if (is.infinite(logdetV)) 
        logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
      if (is.infinite(logdetV)) 
        return(10^10)
    }
    if (REML == TRUE) {
      s2.conc <- t(H) %*% iV %*% H/(n - p)
      LL <- 0.5 * ((n - p) * log(s2.conc) + logdetV + (n - 
                                                         p) + determinant(t(X) %*% iV %*% X)$modulus[1])
    }
    else {
      s2.conc <- t(H) %*% iV %*% H/n
      LL <- 0.5 * (n * log(s2.conc) + logdetV + n)
    }
    if (verbose == T) 
      show(c(as.numeric(LL), par))
    return(as.numeric(LL))
  }
  if (is.null(sp) | is.null(site)) 
    stop("Categorical variables for 'sp' and 'site' must be specified")
  nspp <- nlevels(sp)
  nsite <- nlevels(site)
  mf <- model.frame(formula = formula, data = data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(mf)
  re <- random.effects
  q <- length(re)
  Ztt <- list(NULL)
  St.lengths <- array(0, q)
  nested <- list(NULL)
  ii <- 0
  jj <- 0
  for (i in 1:q) {
    re.i <- re[[i]]
    if (length(re.i) == 3) {
      counter <- 0
      Z.i <- matrix(0, nrow = nspp * nsite, ncol = nlevels(re.i[[2]]))
      for (i.levels in levels(re.i[[2]])) {
        counter <- counter + 1
        Z.i[, counter] <- re.i[[1]] * as.numeric(i.levels == 
                                                   re.i[[2]])
      }
      # Zt.i <- chol(re.i[[3]]) %*% t(Z.i)
      Z.i <- matrix(0,nrow=nspp*nsite,ncol=nspp*nsite)
      ss <- function(x,i){
        rr <- (x-1)*28 + i
        cc <- (i-1)*20 + x
        return(c(rr,cc))
      }
      
      for(w in 1:20){
        for(r in 1:28){
          Z.i[ss(w,r)[1],ss(w,r)[2]] <- 1
        }
      }
      
      Zt.i <- chol(re.i[[3]]) %*% Z.i
      ii <- ii + 1
      Ztt[[ii]] <- Zt.i
      St.lengths[ii] <- nlevels(re.i[[2]])
    }
    if (length(re.i) == 4) {
      if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == 
                                                         sp)) {
        if (length(re.i[[1]]) > 1) 
          stop("Nested terms can only be for intercepts")
        nestedsp.j <- re.i[[3]]
        nestedsite.j <- diag(nsite)
        nested.j <- as(kronecker(nestedsite.j, nestedsp.j), 
                       "dgCMatrix")
      }
      if (setequal(levels(re.i[[2]]), levels(site)) && 
          all(re.i[[2]] == site)) {
        if (length(re.i[[1]]) > 1) 
          stop("Nested terms can only be for intercepts")
        nestedsp.j <- diag(nspp)
        nestedsite.j <- re.i[[3]]
        nested.j <- as(kronecker(nestedsite.j, nestedsp.j), 
                       "dgCMatrix")
      }
      jj <- jj + 1
      nested[[jj]] <- nested.j
    }
  }
  q.nonNested <- ii
  q.Nested <- jj
  if (q.nonNested > 0) {
    St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
#     Zt <- matrix(0, nrow = sum(St.lengths), ncol = nspp * 
#                    nsite)
    Zt <- matrix(0,nrow=nspp*nsite,ncol=nspp*nsite)
    count <- 1
    for (i in 1:q.nonNested) {
      St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, 
                                                         nrow = 1, ncol = St.lengths[i])
      St.lengths[i] <- nspp*nsite
      Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
      count <- count + St.lengths[i]
    }
    Zt <- as(Zt, "dgTMatrix")
    St <- matrix(1,nrow=1,ncol=nspp*nsite)
    St <- as(St, "dgTMatrix")
  }
  else {
    Zt <- NULL
    St <- NULL
  }
  p <- ncol(X)
  n <- nrow(X)
  if (!is.null(B.init) & length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
  }
  if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != 
                          p)) & !is.null(s2.init)) {
    B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, 
                       ncol = p))
  }
  if (!is.null(B.init) & is.null(s2.init)) {
    s2.init <- var(lm(formula = formula, data = data)$residuals)/q
  }
  if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != 
                          p)) & is.null(s2.init)) {
    B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, 
                       ncol = p))
    s2.init <- var(lm(formula = formula, data = data)$residuals)/q
  }
  B <- B.init
  s <- as.vector(array(s2.init^0.5, dim = c(1, q)))
  if (q > 1) {
    opt <- optim(fn = plmm.LL, par = s, X = X, Y = Y, Zt = Zt, 
                 St = St, nested = nested, REML = REML, verbose = verbose, 
                 method = "Nelder-Mead", control = list(maxit = maxit, 
                                                        reltol = reltol))
  }
  else {
    opt <- optim(fn = plmm.LL, par = s, X = X, Y = Y, Zt = Zt, 
                 St = St, nested = nested, REML = REML, verbose = verbose, 
                 method = "L-BFGS-B", control = list(maxit = maxit))
  }
  par <- abs(Re(opt$par))
  LL <- opt$value
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC <- sr[1] * St[1, ]
    if (length(sr) > 1) 
      for (i in 2:q.nonNested) {
        iC <- iC + sr[i] * St[i, ]
      }
    iC <- as(diag(iC), "dsCMatrix")
    Ut <- iC %*% Zt
    U <- t(Ut)
  }
  else {
    q.nonNested <- 0
    sr <- NULL
  }
  if (is.null(nested[[1]])) {
    q.Nested <- 0
  }
  else {
    q.Nested <- length(nested)
  }
  if (q.Nested == 0) {
    sn <- NULL
  }
  else {
    sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
  }
  if (q.Nested == 0) {
    iA <- as(diag(n), "dsCMatrix")
    Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
    Ut.iA.U <- Ut %*% U
    iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
  }
  else {
    A <- as(diag(n), "dsCMatrix")
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    iA <- solve(A)
    if (q.nonNested > 0) {
      Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
      Ut.iA.U <- Ut %*% iA %*% U
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% 
        Ut %*% iA
    }
    else {
      iV <- iA
    }
  }
  denom <- t(X) %*% iV %*% X
  num <- t(X) %*% iV %*% Y
  B <- solve(denom, num)
  B <- as.matrix(B)
  H <- Y - X %*% B
  if (q.Nested == 0) {
    logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
    if (is.infinite(logdetV)) 
      logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
  }
  else {
    logdetV <- -determinant(iV)$modulus[1]
    if (is.infinite(logdetV)) 
      logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
    if (is.infinite(logdetV)) 
      return(10^10)
  }
  if (REML == TRUE) {
    s2resid <- as.numeric(t(H) %*% iV %*% H/(n - p))
  }
  else {
    s2resid <- as.numeric(t(H) %*% iV %*% H/n)
  }
  s2r <- s2resid * sr^2
  s2n <- s2resid * sn^2
  ss <- c(sr, sn, s2resid^0.5)
  iV <- iV/s2resid
  B.cov <- solve(t(X) %*% iV %*% X)
  B.se <- as.matrix(diag(B.cov))^0.5
  B.zscore <- B/B.se
  B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)
  if (REML == TRUE) {
    logLik <- -0.5 * (n - p) * log(2 * pi) + 0.5 * determinant(t(X) %*% 
                                                                 X)$modulus[1] - LL
  }
  else {
    logLik <- -0.5 * n * log(2 * pi) - LL
  }
  k <- p + q + 1
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + k * (log(n) - log(pi))
  results <- list(formula = formula, data = data, family = family, 
                  random.effects = random.effects, B = B, B.se = B.se, 
                  B.cov = B.cov, B.zscore = B.zscore, B.pvalue = B.pvalue, 
                  ss = ss, s2n = s2n, s2r = s2r, s2resid = s2resid, logLik = logLik, 
                  AIC = AIC, BIC = BIC, REML = REML, s2.init = s2.init, 
                  B.init = B.init, Y = Y, X = X, H = H, iV = iV, mu = NULL, 
                  nested = nested, sp = sp, site = site, Zt = Zt, St = St, 
                  convcode = opt$convergence, niter = opt$counts)
  class(results) <- "communityPGLMM"
  results
}


hacked_pez_bad <- function (formula, data = list(), family = "gaussian", sp = NULL, 
                        site = NULL, random.effects = list(), REML = TRUE, s2.init = NULL, 
                        B.init = NULL, reltol = 10^-6, maxit = 500, tol.pql = 10^-6, 
                        maxit.pql = 200, verbose = FALSE) 
{
  if (!is.factor(sp)) 
    stop("'sp' must be a factor")
  if (!is.factor(site)) 
    stop("'site' must be a factor")
  if (family == "gaussian") 
    z <- hacked_gau_bad(formula = formula, data = data, 
                    sp = sp, site = site, random.effects = random.effects, 
                    REML = REML, s2.init = s2.init, B.init = B.init, 
                    reltol = reltol, maxit = maxit, verbose = verbose)
  if (family == "binomial") {
    if (is.null(s2.init)) 
      s2.init <- 0.25
    z <- communityPGLMM.binary(formula = formula, data = data, 
                               sp = sp, site = site, random.effects = random.effects, 
                               REML = REML, s2.init = s2.init, B.init = B.init, 
                               reltol = reltol, maxit = maxit, tol.pql = tol.pql, 
                               maxit.pql = maxit.pql, verbose = verbose)
  }
  if (!is.element(family, c("gaussian", "binomial"))) 
    cat("\nSorry, but only binomial (binary) and gaussian options exist at this time")
  return(z)
}

hacked_gau_bad <- function (formula, data = list(), family = "gaussian", sp = NULL, 
                        site = NULL, random.effects = list(), REML = TRUE, s2.init = NULL, 
                        B.init = NULL, reltol = 10^-8, maxit = 500, verbose = FALSE) 
{
  plmm.LL <- function(par, X, Y, Zt, St, nested = NULL, REML, 
                      verbose) {
    n <- dim(X)[1]
    p <- dim(X)[2]
    if (!is.null(St)) {
      q.nonNested <- dim(St)[1]
      sr <- Re(par[1:q.nonNested])
      iC <- sr[1] * St[1, ]
      if (length(sr) > 1) 
        for (i in 2:q.nonNested) {
          iC <- iC + sr[i] * St[i, ]
        }
      iC <- as(diag(iC), "dsCMatrix")
      Ut <- iC %*% Zt
      U <- t(Ut)
    }
    else {
      q.nonNested <- 0
      sr <- NULL
    }
    if (is.null(nested[[1]])) {
      q.Nested <- 0
    }
    else {
      q.Nested <- length(nested)
    }
    if (q.Nested == 0) {
      sn <- NULL
    }
    else {
      sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
    }
    if (q.Nested == 0) {
      iA <- as(diag(n), "dsCMatrix")
      Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
      Ut.iA.U <- Ut %*% U
      iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
    }
    else {
      A <- as(diag(n), "dsCMatrix")
      for (j in 1:q.Nested) {
        A <- A + sn[j]^2 * nested[[j]]
      }
      iA <- solve(A)
      if (q.nonNested > 0) {
        Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
        Ut.iA.U <- Ut %*% iA %*% U
        iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% 
          Ut %*% iA
      }
      else {
        iV <- iA
      }
    }
    denom <- t(X) %*% iV %*% X
    num <- t(X) %*% iV %*% Y
    B <- solve(denom, num)
    B <- as.matrix(B)
    H <- Y - X %*% B
    if (q.Nested == 0) {
      logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
      if (is.infinite(logdetV)) 
        logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
    }
    else {
      logdetV <- -determinant(iV)$modulus[1]
      if (is.infinite(logdetV)) 
        logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
      if (is.infinite(logdetV)) 
        return(10^10)
    }
    if (REML == TRUE) {
      s2.conc <- t(H) %*% iV %*% H/(n - p)
      LL <- 0.5 * ((n - p) * log(s2.conc) + logdetV + (n - 
                                                         p) + determinant(t(X) %*% iV %*% X)$modulus[1])
    }
    else {
      s2.conc <- t(H) %*% iV %*% H/n
      LL <- 0.5 * (n * log(s2.conc) + logdetV + n)
    }
    if (verbose == T) 
      show(c(as.numeric(LL), par))
    return(as.numeric(LL))
  }
  if (is.null(sp) | is.null(site)) 
    stop("Categorical variables for 'sp' and 'site' must be specified")
  nspp <- nlevels(sp)
  nsite <- nlevels(site)
  mf <- model.frame(formula = formula, data = data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(mf)
  re <- random.effects
  q <- length(re)
  Ztt <- list(NULL)
  St.lengths <- array(0, q)
  nested <- list(NULL)
  ii <- 0
  jj <- 0
  for (i in 1:q) {
    re.i <- re[[i]]
    if (length(re.i) == 3) {
      counter <- 0
      Z.i <- matrix(0, nrow = nspp * nsite, ncol = nlevels(re.i[[2]]))
      for (i.levels in levels(re.i[[2]])) {
        counter <- counter + 1
        Z.i[, counter] <- re.i[[1]] * as.numeric(i.levels == 
                                                   re.i[[2]])
      }
      # Zt.i <- chol(re.i[[3]]) %*% t(Z.i)
      Z.i <- matrix(0,nrow=nspp*nsite,ncol=nspp*nsite)
      ss <- function(x,i){
        rr <- (x-1)*28 + i
        cc <- (i-1)*20 + x
        return(c(rr,cc))
      }
      
      for(w in 1:20){
        for(r in 1:28){
          Z.i[ss(w,r)[1],ss(w,r)[2]] <- 1
        }
      }
      
      Zt.i <- chol(re.i[[3]]) #%*% Z.i
      ii <- ii + 1
      Ztt[[ii]] <- Zt.i
      St.lengths[ii] <- nlevels(re.i[[2]])
    }
    if (length(re.i) == 4) {
      if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == 
                                                         sp)) {
        if (length(re.i[[1]]) > 1) 
          stop("Nested terms can only be for intercepts")
        nestedsp.j <- re.i[[3]]
        nestedsite.j <- diag(nsite)
        nested.j <- as(kronecker(nestedsite.j, nestedsp.j), 
                       "dgCMatrix")
      }
      if (setequal(levels(re.i[[2]]), levels(site)) && 
          all(re.i[[2]] == site)) {
        if (length(re.i[[1]]) > 1) 
          stop("Nested terms can only be for intercepts")
        nestedsp.j <- diag(nspp)
        nestedsite.j <- re.i[[3]]
        nested.j <- as(kronecker(nestedsite.j, nestedsp.j), 
                       "dgCMatrix")
      }
      jj <- jj + 1
      nested[[jj]] <- nested.j
    }
  }
  q.nonNested <- ii
  q.Nested <- jj
  if (q.nonNested > 0) {
    St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
    #     Zt <- matrix(0, nrow = sum(St.lengths), ncol = nspp * 
    #                    nsite)
    Zt <- matrix(0,nrow=nspp*nsite,ncol=nspp*nsite)
    count <- 1
    for (i in 1:q.nonNested) {
      St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, 
                                                         nrow = 1, ncol = St.lengths[i])
      St.lengths[i] <- nspp*nsite
      Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
      count <- count + St.lengths[i]
    }
    Zt <- as(Zt, "dgTMatrix")
    St <- matrix(1,nrow=1,ncol=nspp*nsite)
    St <- as(St, "dgTMatrix")
  }
  else {
    Zt <- NULL
    St <- NULL
  }
  p <- ncol(X)
  n <- nrow(X)
  if (!is.null(B.init) & length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
  }
  if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != 
                          p)) & !is.null(s2.init)) {
    B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, 
                       ncol = p))
  }
  if (!is.null(B.init) & is.null(s2.init)) {
    s2.init <- var(lm(formula = formula, data = data)$residuals)/q
  }
  if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != 
                          p)) & is.null(s2.init)) {
    B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, 
                       ncol = p))
    s2.init <- var(lm(formula = formula, data = data)$residuals)/q
  }
  B <- B.init
  s <- as.vector(array(s2.init^0.5, dim = c(1, q)))
  if (q > 1) {
    opt <- optim(fn = plmm.LL, par = s, X = X, Y = Y, Zt = Zt, 
                 St = St, nested = nested, REML = REML, verbose = verbose, 
                 method = "Nelder-Mead", control = list(maxit = maxit, 
                                                        reltol = reltol))
  }
  else {
    opt <- optim(fn = plmm.LL, par = s, X = X, Y = Y, Zt = Zt, 
                 St = St, nested = nested, REML = REML, verbose = verbose, 
                 method = "L-BFGS-B", control = list(maxit = maxit))
  }
  par <- abs(Re(opt$par))
  LL <- opt$value
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC <- sr[1] * St[1, ]
    if (length(sr) > 1) 
      for (i in 2:q.nonNested) {
        iC <- iC + sr[i] * St[i, ]
      }
    iC <- as(diag(iC), "dsCMatrix")
    Ut <- iC %*% Zt
    U <- t(Ut)
  }
  else {
    q.nonNested <- 0
    sr <- NULL
  }
  if (is.null(nested[[1]])) {
    q.Nested <- 0
  }
  else {
    q.Nested <- length(nested)
  }
  if (q.Nested == 0) {
    sn <- NULL
  }
  else {
    sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
  }
  if (q.Nested == 0) {
    iA <- as(diag(n), "dsCMatrix")
    Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
    Ut.iA.U <- Ut %*% U
    iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
  }
  else {
    A <- as(diag(n), "dsCMatrix")
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    iA <- solve(A)
    if (q.nonNested > 0) {
      Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
      Ut.iA.U <- Ut %*% iA %*% U
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% 
        Ut %*% iA
    }
    else {
      iV <- iA
    }
  }
  denom <- t(X) %*% iV %*% X
  num <- t(X) %*% iV %*% Y
  B <- solve(denom, num)
  B <- as.matrix(B)
  H <- Y - X %*% B
  if (q.Nested == 0) {
    logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
    if (is.infinite(logdetV)) 
      logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
  }
  else {
    logdetV <- -determinant(iV)$modulus[1]
    if (is.infinite(logdetV)) 
      logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
    if (is.infinite(logdetV)) 
      return(10^10)
  }
  if (REML == TRUE) {
    s2resid <- as.numeric(t(H) %*% iV %*% H/(n - p))
  }
  else {
    s2resid <- as.numeric(t(H) %*% iV %*% H/n)
  }
  s2r <- s2resid * sr^2
  s2n <- s2resid * sn^2
  ss <- c(sr, sn, s2resid^0.5)
  iV <- iV/s2resid
  B.cov <- solve(t(X) %*% iV %*% X)
  B.se <- as.matrix(diag(B.cov))^0.5
  B.zscore <- B/B.se
  B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)
  if (REML == TRUE) {
    logLik <- -0.5 * (n - p) * log(2 * pi) + 0.5 * determinant(t(X) %*% 
                                                                 X)$modulus[1] - LL
  }
  else {
    logLik <- -0.5 * n * log(2 * pi) - LL
  }
  k <- p + q + 1
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + k * (log(n) - log(pi))
  results <- list(formula = formula, data = data, family = family, 
                  random.effects = random.effects, B = B, B.se = B.se, 
                  B.cov = B.cov, B.zscore = B.zscore, B.pvalue = B.pvalue, 
                  ss = ss, s2n = s2n, s2r = s2r, s2resid = s2resid, logLik = logLik, 
                  AIC = AIC, BIC = BIC, REML = REML, s2.init = s2.init, 
                  B.init = B.init, Y = Y, X = X, H = H, iV = iV, mu = NULL, 
                  nested = nested, sp = sp, site = site, Zt = Zt, St = St, 
                  convcode = opt$convergence, niter = opt$counts)
  class(results) <- "communityPGLMM"
  results
}
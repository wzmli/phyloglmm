glmmTMBhacked <- function (formula, data = NULL, family = gaussian(), ziformula = ~0, 
          dispformula = ~1, weights = NULL, offset = NULL, se = TRUE, phyloZ = NULL
          , phylonm = NULL ,
          verbose = FALSE, doFit = TRUE, control = glmmTMBControl()) 
{
  call <- mf <- mc <- match.call()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (!all(c("family", "link") %in% names(family))) 
    stop("'family' must contain at least 'family' and 'link' components")
  if (grepl("^quasi", family$family)) 
    stop("\"quasi\" families cannot be used in glmmTMB")
  link <- family$link
  environment(formula) <- parent.frame()
  call$formula <- mc$formula <- formula
  environment(ziformula) <- environment(formula)
  call$ziformula <- ziformula
  environment(dispformula) <- environment(formula)
  call$dispformula <- dispformula
  m <- match(c("data", "subset", "weights", "na.action", "offset"), 
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  if (glmmTMB:::inForm(ziformula, quote(.))) {
    ziformula <- update(RHSForm(drop.special2(formula), as.form = TRUE), 
                        ziformula)
  }
  formList <- list(formula, ziformula, dispformula)
  formList <- lapply(formList, function(x) noSpecials(lme4:::subbars(x), 
                                                      delete = FALSE))
  combForm <- do.call(addForm, formList)
  environment(combForm) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) 
      assign(i, get(i, parent.frame()), environment(combForm))
  }
  mf$formula <- combForm
  fr <- eval(mf, envir = environment(formula), enclos = parent.frame())
  nobs <- nrow(fr)
  weights <- as.vector(model.weights(fr))
  if (!is.null(weights) & !glmmTMB:::okWeights(family$family)) {
    stop("'weights' are not available for this family. See `dispformula` instead.")
  }
  if (is.null(weights)) 
    weights <- rep(1, nobs)
  respCol <- attr(terms(fr), "response")
  names(respCol) <- names(fr)[respCol]
  y <- fr[, respCol]
  if (is.matrix(y)) {
    if (family$family == "binomial") {
      warning("binomial models with N>1 are preferably specified as proportion~..., weights=N")
    }
    else {
      stop("matrix-valued responses are not allowed")
    }
  }
  etastart <- start <- mustart <- NULL
  if (!is.null(family$initialize)) {
    eval(family$initialize)
  }
  y <- as.numeric(y)
  if (grepl("^truncated", family$family) & (any(y < 1)) & (ziformula == 
                                                           ~0)) 
    stop(paste0("'", names(respCol), "'", " contains zeros (or values below the allowable range). ", 
                "Zeros are compatible with a trucated distribution only when zero-inflation is added."))
  TMBStruc <- mkTMBStruchacked(formula, ziformula, dispformula, combForm, 
                         mf, fr, yobs = y, respCol, offset, weights, family = family, 
                         se = se, phyloZ = phyloZ, phylonm = phylonm ,call = call, verbose = verbose)
  TMBStruc$control <- lapply(control, eval, envir = TMBStruc)
  if (!doFit) 
    return(TMBStruc)
  # res <- fit_TMBstruc(TMBStruc)
  res <- glmmTMB:::fitTMB(TMBStruc)
  return(res)
}

mkTMBStruchacked <- function (formula, ziformula, dispformula, combForm, mf, fr, 
          yobs, respCol, offset, weights, family, se = NULL, phyloZ = phyloZ, phylonm = phylonm, 
          call = NULL, verbose = NULL, ziPredictCode = "corrected", doPredict = 0, 
          whichPredict = integer(0)) 
{
  mapArg <- NULL
  if (glmmTMB:::usesDispersion(family$family) && (dispformula == ~0)) {
    if (family != "gaussian") 
      stop("~0 dispersion not implemented for ", sQuote(family$family), 
           " family")
    betad_init <- log(sqrt(.Machine$double.eps))
    dispformula[] <- ~1
    mapArg <- list(betad = factor(NA))
  }
  else {
    betad_init <- 0
  }
  if (!glmmTMB:::usesDispersion(family$family)) {
    dispformula[] <- ~0
  }
  condList <- getXReTrmshacked(formula, mf, fr, phyloZ=phyloZ, phylonm=phylonm)
  ziList <- glmmTMB:::getXReTrms(ziformula, mf, fr)
  dispList <- glmmTMB:::getXReTrms(dispformula, mf, fr, ranOK = FALSE, 
                         "dispersion")
  ziReStruc <- with(ziList, getReStruc(reTrms, ss))
  grpVar <- with(condList, getGrpVar(reTrms$flist))
  condReStruc <- with(condList, getReStruc(reTrms, ss))
  nobs <- nrow(fr)
  if (is.null(offset <- model.offset(fr))) 
    offset <- rep(0, nobs)
  if (is.null(weights <- fr[["(weights)"]])) 
    weights <- rep(1, nobs)
  data.tmb <- lme4:::namedList(X = condList$X, Z = condList$Z, Xzi = ziList$X, 
                        Zzi = ziList$Z, Xd = dispList$X, yobs, respCol, offset, 
                        weights, terms = condReStruc, termszi = ziReStruc, family = glmmTMB:::.valid_family[family$family], 
                        link = glmmTMB:::.valid_link[family$link], ziPredictCode = glmmTMB:::.valid_zipredictcode[ziPredictCode], 
                        doPredict = doPredict, whichPredict = whichPredict)
  getVal <- function(obj, component) vapply(obj, function(x) x[[component]], 
                                            numeric(1))
  beta_init <- if (family$link %in% c("identity", "inverse")) 
    1
  else 0
  numThetaFamily <- (family$family == "tweedie")
  parameters <- with(data.tmb, list(beta = rep(beta_init, ncol(X)), 
                                    betazi = rep(0, ncol(Xzi)), b = rep(beta_init, ncol(Z)), 
                                    bzi = rep(0, ncol(Zzi)), betad = rep(betad_init, ncol(Xd)), 
                                    theta = rep(0, sum(getVal(condReStruc, "blockNumTheta"))), 
                                    thetazi = rep(0, sum(getVal(ziReStruc, "blockNumTheta"))), 
                                    thetaf = rep(0, numThetaFamily)))
  randomArg <- c(if (ncol(data.tmb$Z) > 0) "b", if (ncol(data.tmb$Zzi) > 
                                                    0) "bzi")
  return(lme4:::namedList(data.tmb, parameters, mapArg, randomArg, 
                   grpVar, condList, ziList, dispList, condReStruc, ziReStruc, 
                   family, respCol, allForm = lme4:::namedList(combForm, formula, 
                                                        ziformula, dispformula), fr, se, call, verbose))
}

getXReTrmshacked <- function (formula, mf, fr, ranOK = TRUE, type = "", phyloZ=phyloZ, phylonm=phylonm) 
{
  fixedform <- formula
  glmmTMB:::RHSForm(fixedform) <- lme4:::nobars(glmmTMB:::RHSForm(fixedform))
  nobs <- nrow(fr)
  if (identical(glmmTMB:::RHSForm(fixedform), ~0) || identical(glmmTMB:::RHSForm(fixedform), 
                                                     ~-1)) {
    X <- NULL
  }
  else {
    mf$formula <- fixedform
    terms_fixed <- terms(eval(mf, envir = environment(fixedform)))
    X <- model.matrix(fixedform, fr, contrasts)
    terms <- list(fixed = terms(terms_fixed))
  }
  ranform <- formula
  if (is.null(lme4:::findbars(ranform))) {
    reTrms <- NULL
    Z <- new("dgCMatrix", Dim = c(as.integer(nobs), 0L))
    ss <- integer(0)
  }
  else {
    if (!ranOK) 
      stop("no random effects allowed in ", type, " term")
    glmmTMB:::RHSForm(ranform) <- lme4:::subbars(glmmTMB:::RHSForm(glmmTMB:::reOnly(formula)))
    mf$formula <- ranform
    reTrms <- mkReTrms(lme4:::findbars(glmmTMB:::RHSForm(formula)), fr,phylonm,phyloZ)
    ss <- splitForm(formula)
    ss <- unlist(ss$reTrmClasses)
    Z <- t(reTrms$Zt)
  }
  lme4:::namedList(X, Z, reTrms, ss, terms)
}

getRestruchacked <- function (reTrms, ss = NULL) 
{
  if (is.null(reTrms)) {
    list()
  }
  else {
    assign <- attr(reTrms$flist, "assign")
    nreps <- vapply(assign, function(i) length(levels(reTrms$flist[[i]])), 
                    0)
    blksize <- diff(reTrms$Gp)/nreps
    if (is.null(ss)) {
      ss <- rep("us", length(blksize))
    }
    covCode <- .valid_covstruct[ss]
    parFun <- function(struc, blksize) {
      switch(as.character(struc), `0` = blksize, `1` = blksize * 
               (blksize + 1)/2, `2` = blksize + 1, `3` = 2, 
             `4` = 2, `5` = 2, `6` = 2, `7` = 3, `8` = 2 * 
               blksize - 1)
    }
    blockNumTheta <- mapply(parFun, covCode, blksize, SIMPLIFY = FALSE)
    ans <- lapply(seq_along(ss), function(i) {
      tmp <- list(blockReps = nreps[i], blockSize = blksize[i], 
                  blockNumTheta = blockNumTheta[[i]], blockCode = covCode[i])
      if (ss[i] == "ar1") {
        if (any(reTrms$cnms[[i]][1] == "(Intercept)")) 
          warning("AR1 not meaningful with intercept")
      }
      if (ss[i] == "ou") {
        times <- parseNumLevels(reTrms$cnms[[i]])
        if (ncol(times) != 1) 
          stop("'ou' structure is for 1D coordinates only.")
        if (is.unsorted(times, strictly = TRUE)) 
          stop("'ou' is for strictly sorted times only.")
        tmp$times <- drop(times)
      }
      if (ss[i] %in% c("exp", "gau", "mat")) {
        coords <- parseNumLevels(reTrms$cnms[[i]])
        tmp$dist <- as.matrix(dist(coords))
      }
      tmp
    })
    setNames(ans, names(reTrms$Ztlist))
  }
}
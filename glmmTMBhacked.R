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
                         se = se, call = call, verbose = verbose)
  TMBStruc$control <- lapply(control, eval, envir = TMBStruc)
  if (!doFit) 
    return(TMBStruc)
  res <- fitTMB(TMBStruc)
  return(res)
}

mkTMBStruchacked <- function (formula, ziformula, dispformula, combForm, mf, fr, 
          yobs, respCol, offset, weights, family, se = NULL, call = NULL, 
          verbose = NULL, ziPredictCode = "corrected", doPredict = 0, 
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
  condList <- glmmTMB:::getXReTrms(formula, mf, fr)
  ziList <- getXReTrms(ziformula, mf, fr)
  dispList <- getXReTrms(dispformula, mf, fr, ranOK = FALSE, 
                         "dispersion")
  condReStruc <- with(condList, getReStruc(reTrms, ss))
  ziReStruc <- with(ziList, getReStruc(reTrms, ss))
  grpVar <- with(condList, getGrpVar(reTrms$flist))
  nobs <- nrow(fr)
  if (is.null(offset <- model.offset(fr))) 
    offset <- rep(0, nobs)
  if (is.null(weights <- fr[["(weights)"]])) 
    weights <- rep(1, nobs)
  data.tmb <- namedList(X = condList$X, Z = condList$Z, Xzi = ziList$X, 
                        Zzi = ziList$Z, Xd = dispList$X, yobs, respCol, offset, 
                        weights, terms = condReStruc, termszi = ziReStruc, family = .valid_family[family$family], 
                        link = .valid_link[family$link], ziPredictCode = .valid_zipredictcode[ziPredictCode], 
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
  return(namedList(data.tmb, parameters, mapArg, randomArg, 
                   grpVar, condList, ziList, dispList, condReStruc, ziReStruc, 
                   family, respCol, allForm = namedList(combForm, formula, 
                                                        ziformula, dispformula), fr, se, call, verbose))
}
#' hacked version of glmmTMB
## @params copied from glmmTMB as well ...
##' @inheritParams phylo_lmm
##' @param formula combined fixed and random effects formula, following lme4 syntax.
##' @param data optional data frame containing model variables.
##' @param family a family function, a character string naming a family function, or the result of a call to a family function (variance/link function) information. See \code{\link{family}} for a generic discussion of families or \code{\link{family_glmmTMB}} for details of \code{glmmTMB}-specific families.
##' @param ziformula a \emph{one-sided} (i.e., no response variable) formula for zero-inflation combining fixed and random effects: the default \code{~0} specifies no zero-inflation. Specifying \code{~.} sets the zero-inflation formula identical to the right-hand side of \code{formula} (i.e., the conditional effects formula); terms can also be added or subtracted. \strong{When using \code{~.} as the zero-inflation formula in models where the conditional effects formula contains an offset term, the offset term will automatically be dropped}. The zero-inflation model uses a logit link.
##' @param dispformula a \emph{one-sided} formula for dispersion containing only fixed effects: the default \code{~1} specifies the standard dispersion given any family. The argument is ignored for families that do not have a dispersion parameter. For an explanation of the dispersion parameter for each family, see \code{\link{sigma}}. The dispersion model uses a log link. In Gaussian mixed models, \code{dispformula=~0} fixes the residual variance to be 0 (actually a small non-zero value), forcing variance into the random effects. The precise value can be controlled via \code{control=glmmTMBControl(zero_dispval=...)}; the default value is \code{sqrt(.Machine$double.eps)}.
##' @param weights weights, as in \code{glm}. Not automatically scaled to have sum 1.
##' @param offset offset for conditional model (only).
##' @param contrasts an optional list, e.g., \code{list(fac1="contr.sum")}. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
##' @param na.action how to handle missing values, see \code{\link{na.action}} and \code{\link{model.frame}}. From \code{\link{lm}}: \dQuote{The default is set by the \code{\link{na.action}} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The \sQuote{factory-fresh} default is \code{\link{na.omit}}.}
##' @param se whether to return standard errors.
##' @param verbose whether progress indication should be printed to the console.
##' @param doFit whether to fit the full model, or (if FALSE) return the preprocessed data and parameter objects, without fitting the model.
##' @param control control parameters, see \code{\link{glmmTMBControl}}.
##' @param REML whether to use REML estimation rather than maximum likelihood.
##' @param start starting values, expressed as a list with possible components \code{beta}, \code{betazi}, \code{betad} (fixed-effect parameters for conditional, zero-inflation, dispersion models); \code{b}, \code{bzi} (conditional modes for conditional and zero-inflation models); \code{theta}, \code{thetazi} (random-effect parameters, on the standard deviation/Cholesky scale, for conditional and z-i models); \code{thetaf} (extra family parameters, e.g., shape for Tweedie models).
##' @param map a list specifying which parameter values should be fixed to a constant value rather than estimated. \code{map} should be a named list containing factors corresponding to a subset of the internal parameter names (see \code{start} parameter). Distinct factor values are fitted as separate parameter values, \code{NA} values are held fixed: e.g., \code{map=list(beta=factor(c(1,2,3,NA)))} would fit the first three fixed-effect parameters of the conditional model and fix the fourth parameter to its starting value. In general, users will probably want to use \code{start} to specify non-default starting values for fixed parameters. See \code{\link[TMB]{MakeADFun}} for more details.
##' @param sparseX a named logical vector containing (possibly) elements named "cond", "zi", "disp" to indicate whether fixed-effect model matrices for particular model components should be generated as sparse matrices, e.g. \code{c(cond=TRUE)}. Default is all \code{FALSE}
##' @importFrom stats make.link na.fail update as.formula terms model.weights gaussian model.matrix
##' @importFrom methods is new
##' @importFrom glmmTMB inForm extractForm dropHead noSpecials addForm splitForm getReStruc glmmTMBControl getGrpVar fitTMB
#' @export
phylo_glmmTMB <-  function(formula, data = NULL, family = gaussian(), ziformula = ~0,
                         dispformula = ~1, weights = NULL, offset = NULL, contrasts = NULL, phylo = NULL, phyloZ = NULL,
                         phylonm = NULL,
                         na.action = na.fail, se = TRUE, verbose = FALSE, doFit = TRUE,
                         control = glmmTMBControl(), REML = FALSE, map = NULL, sparseX = NULL) {

  call <- mf <- mc <- match.call()

  ## include with hack_function
  phyloZ <- get_phyloZ(phylo, phyloZ, data[[phylonm]])

  if (is.character(family)) {
    if (family == "beta") {
      family <- "beta_family"
      warning(
        "please use ", sQuote("beta_family()"), " rather than ",
        sQuote("\"beta\""), " to specify a Beta-distributed response"
      )
    }
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  fnames <- names(family)
  if (!all(c("family", "link") %in% fnames)) {
    stop("'family' must contain at least 'family' and 'link' components")
  }
  if (length(miss_comp <- setdiff(
    c("linkfun", "variance"),
    fnames
  )) > 0) {
    warning(
      "some components missing from ", sQuote("family"),
      ": downstream methods may fail"
    )
  }
  if (grepl("^quasi", family$family)) {
    stop("\"quasi\" families cannot be used in glmmTMB")
  }
  link <- family$link
  environment(formula) <- parent.frame()
  call$formula <- mc$formula <- formula
  if (!is.null(eval(substitute(offset), data, enclos = environment(formula)))) {
    formula <- glmmTMB::addForm0(formula, glmmTMB:::makeOp(substitute(offset),
      op = quote(offset)
    ))
  }
  environment(ziformula) <- environment(formula)
  call$ziformula <- ziformula
  environment(dispformula) <- environment(formula)
  call$dispformula <- dispformula
  m <- match(
    c("data", "subset", "weights", "na.action", "offset"),
    names(mf), 0L
  )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  if (inForm(ziformula, quote(.))) {
    ziformula <- update(
      glmmTMB::RHSForm(glmmTMB::drop.special(formula), as.form = TRUE),
      ziformula
    )
  }
  formList <- list(formula, ziformula, dispformula)
  for (i in seq_along(formList)) {
    f <- formList[[i]]
    f <- noSpecials(subbars(f), delete = FALSE)
    formList[[i]] <- f
  }
  combForm <- do.call(addForm, formList)
  environment(combForm) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) {
      assign(i, get(i, parent.frame()), environment(combForm))
    }
  }
  mf$formula <- combForm
  fr <- eval(mf, envir = environment(formula), enclos = parent.frame())
  nobs <- nrow(fr)
  weights <- as.vector(model.weights(fr))
  if (!is.null(weights) & !glmmTMB:::okWeights(family$family)) {
    stop("'weights' are not available for this family.")
  }
  if (is.null(weights)) {
    weights <- rep(1, nobs)
  }
  respCol <- attr(terms(fr), "response")
  names(respCol) <- names(fr)[respCol]
  y <- fr[, respCol]
  if (is.matrix(y)) {
    if (!glmmTMB:::binomialType(family$family)) {
      stop("matrix-valued responses are not allowed")
    }
  }
  etastart <- start <- mustart <- NULL
  if (!is.null(family$initialize)) {
    local(eval(family$initialize))
  }
  if (grepl("^truncated", family$family) && (!is.factor(y) &&
    any(y < 1)) & (ziformula == ~0)) {
    stop(paste0(
      "'", names(respCol), "'", " contains zeros (or values below the allowable range). ",
      "Zeros are compatible with a trucated distribution only when zero-inflation is added."
    ))
  }
  ## hack_function: start at TMBStruc <-
  ## end with control=control)
  TMBStruc <- mkTMBStrucphylo(formula, ziformula, dispformula, combForm,
    mf, fr,
    yobs = y, respCol, weights, contrasts = contrasts,
    family = family, se = se, call = call, verbose = verbose,
    REML = REML, phylonm = phylonm, phyloZ = phyloZ, map = map, sparseX = sparseX
  )
  TMBStruc$control <- lapply(control, eval, envir = TMBStruc)
  if (!doFit) {
    return(TMBStruc)
  }
  res <- glmmTMB::fitTMB(TMBStruc)
  return(res)
}

# make TMB structure
# don't export!
mkTMBStrucphylo <- function(formula, ziformula, dispformula, combForm, mf, fr,
                            yobs, respCol, weights, contrasts = contrasts, size = NULL, family, se = NULL, phyloZ = phyloZ,
                            phylonm = phylonm,
                            call = NULL, verbose = NULL, ziPredictCode = "corrected",
                            doPredict = 0, whichPredict = integer(0), REML = FALSE, map = NULL, sparseX = NULL) {
  if (!is(family, "family")) {
    if (is.list(family)) {
      warning("specifying ", sQuote("family"), " as a plain list is deprecated")
    }
    fname <- family$family
    args <- family["link"]
    ff <- try(do.call(fname, args), silent = TRUE)
    if (!inherits(ff, "try-error")) {
      family <- ff
    }
    else {
      if (is.null(family$linkfun)) {
        family <- c(family, make.link(family$link))
      }
    }
  }
  mapArg <- NULL
  dispformula.orig <- dispformula
  if (glmmTMB:::usesDispersion(family$family) && (dispformula == ~0)) {
    if (family$family != "gaussian") {
      stop(
        "~0 dispersion not implemented for ", sQuote(family$family),
        " family"
      )
    }
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
  condList <- glmmTMB:::getXReTrms(formula, mf, fr, contrasts = contrasts)
  ## hack_function: add after "condList <- getXReTerms"
  condListphylo <- getXReTrmsphylo(formula, mf, fr, contrasts = contrasts, phylonm = phylonm, phyloZ = phyloZ)
  ziList <- glmmTMB:::getXReTrms(ziformula, mf, fr, contrasts = contrasts)
  dispList <- glmmTMB:::getXReTrms(dispformula, mf, fr,
    ranOK = FALSE,
    type = "dispersion", contrasts = contrasts
  )
  condReStruc <- with(condList, getReStruc(reTrms, ss))
  ziReStruc <- with(ziList, getReStruc(reTrms, ss))
  grpVar <- with(condList, getGrpVar(reTrms$flist))
  nobs <- nrow(fr)
  if (is.null(weights)) {
    weights <- rep(1, nobs)
  }
  if (glmmTMB:::binomialType(family$family)) {
    if (is.factor(yobs)) {
      yobs <- pmin(as.numeric(yobs) - 1, 1)
      size <- rep(1, nobs)
    }
    else {
      if (is.matrix(yobs)) {
        size <- yobs[, 1] + yobs[, 2]
        yobs <- yobs[, 1]
      }
      else {
        if (all(yobs %in% c(0, 1))) {
          size <- rep(1, nobs)
        }
        else {
          yobs <- weights * yobs
          size <- weights
          weights <- rep(1, nobs)
        }
      }
    }
  }
  if (is.null(size)) {
    size <- numeric(0)
  }
  data.tmb <- lme4:::namedList(
    X = condList$X, Z = condList$Z, Xzi = ziList$X,
    Zzi = ziList$Z, Xd = dispList$X, yobs, respCol, offset = condList$offset,
    zioffset = ziList$offset, doffset = dispList$offset,
    weights, size, terms = condReStruc, termszi = ziReStruc,
    family = glmmTMB:::.valid_family[family$family], link = glmmTMB:::.valid_link[family$link],
    ziPredictCode = glmmTMB:::.valid_zipredictcode[ziPredictCode],
    doPredict = doPredict, whichPredict = whichPredict
  )
  getVal <- function(obj, component) {
    vapply(
      obj, function(x) x[[component]],
      numeric(1)
    )
  }
  beta_init <- if (family$link %in% c("identity", "inverse")) {
    1
  } else {
    0
  }
  numThetaFamily <- (family$family == "tweedie")
  parameters <- with(data.tmb, list(
    beta = rep(beta_init, ncol(X)),
    betazi = rep(0, ncol(Xzi)), b = rep(beta_init, ncol(Z)),
    bzi = rep(0, ncol(Zzi)), betad = rep(betad_init, ncol(Xd)),
    theta = rep(0, sum(getVal(condReStruc, "blockNumTheta"))),
    thetazi = rep(0, sum(getVal(ziReStruc, "blockNumTheta"))),
    thetaf = rep(0, numThetaFamily)
  ))
  randomArg <- c(if (ncol(data.tmb$Z) > 0) "b", if (ncol(data.tmb$Zzi) >
    0) {
    "bzi"
  })
  if (REML) {
    randomArg <- c(randomArg, "beta")
  }
  ## hack_function: add after 'if (REML) randomArg <- c(randomArg, "beta")'
  n.edge <- ncol(phyloZ)
  n.site <- length(unique(fr[["site"]]))
  relength <- length(condReStruc)
  if (relength == 1) {
    REname <- unlist(strsplit(names(condReStruc), " "))
    rightbar <- REname[length(REname)]
    if (rightbar %in% phylonm) {
      condReStruc[[1]]$blockReps <- n.edge
      data.tmb$terms[[1]]$blockReps <- n.edge
    }
  }
  if (relength > 1) {
    for (i in seq_along(condReStruc)) {
      REname <- unlist(strsplit(names(condReStruc)[i], " "))
      rightbar <- REname[length(REname)]
      if (rightbar %in% phylonm) {
        condReStruc[[i]]$blockReps <- n.edge
        data.tmb$terms[[i]]$blockReps <- n.edge
        if (rightbar == "sp:site") {
          data.tmb$terms[[i]]$blockReps <- n.edge * n.site ## FIXME: pull out number of site from somewhere
        }
      }
    }
  }
  condList$Z <- t(condListphylo$reTrms$Zt)
  data.tmb$Z <- t(condListphylo$reTrms$Zt)
  parameters$b <- rep(0, ncol(data.tmb$Z))
  ## hack_function: add up to here
  dispformula <- dispformula.orig
  return(lme4:::namedList(data.tmb, parameters, mapArg, randomArg,
    grpVar, condList, ziList, dispList, condReStruc, ziReStruc,
    family, contrasts, respCol,
    allForm = lme4:::namedList(
      combForm,
      formula, ziformula, dispformula
    ), fr, se, call, verbose,
    REML, map, sparseX
  ))
}

getXReTrmsphylo <- function(formula, mf, fr, ranOK = TRUE, type = "", contrasts, phyloZ = phyloZ,
                            phylonm = phylonm) {
  fixedform <- formula
  glmmTMB::RHSForm(fixedform) <- nobars(glmmTMB::RHSForm(fixedform))
  nobs <- nrow(fr)
  if (identical(glmmTMB::RHSForm(fixedform), ~0) || identical(
    glmmTMB::RHSForm(fixedform),
    ~ -1
  )) {
    X <- NULL
  }
  else {
    mf$formula <- fixedform
    terms_fixed <- terms(eval(mf, envir = environment(fixedform)))
    X <- model.matrix(glmmTMB::drop.special(fixedform), fr, contrasts)
    offset <- rep(0, nobs)
    terms <- list(fixed = terms(terms_fixed))
    if (inForm(fixedform, quote(offset))) {
      for (o in extractForm(fixedform, quote(offset))) {
        offset_nm <- deparse(o)
        if (length(offset_nm) > 1) {
          stop("trouble reconstructing offset name")
        }
        offset <- offset + fr[[offset_nm]]
      }
    }
  }
  ranform <- formula
  if (is.null(lme4::findbars(ranform))) {
    reTrms <- NULL
    Z <- new("dgCMatrix", Dim = c(as.integer(nobs), 0L))
    ss <- integer(0)
  }
  else {
    if (!ranOK) {
      stop("no random effects allowed in ", type, " term")
    }
    glmmTMB::RHSForm(ranform) <- subbars(glmmTMB::RHSForm(glmmTMB::reOnly(formula)))
    mf$formula <- ranform
    reTrms <- mkReTrms(lme4::findbars(glmmTMB::RHSForm(formula)), fr, phylonm, phyloZ)
    ss <- splitForm(formula)
    ss <- unlist(ss$reTrmClasses)
    Z <- t(reTrms$Zt)
  }
  lme4:::namedList(X, Z, reTrms, ss, terms, offset)
}

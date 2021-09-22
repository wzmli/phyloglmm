##' phylogenetic linear mixed model
##' @importFrom lme4 optimizeLmer optimizeGlmer mkMerMod mkLmerDevfun mkGlmerDevfun
##' @importFrom lme4 findbars nobars lmerControl glmerControl updateGlmerDevfun subbars
##' @param formula mixed-model formula
##' @param data data frame
##' @param phylo phylogenetic tree in \code{\link{phylo}} format
##' @param phylonm name of phylogenetic grouping variable
##' @param phyloZ phylogenetic Z-matrix (see \code{\link{phylo.to.Z}}): optional, will be computed
##' internally from \code{phylo} if not specified here. (For large phylogenies that are going to
##' be used in multiple models it may make sense to compute the Z matrix first and pass it to
##' the modeling functions to avoid recomputing it each time.)
##' to compute \code{phyloZ} fi
##' @param control control
##' @param REML use restricted max likelihood?
## unimplemented ... pass via ... ?
## @param subset subset expression
## @param weights weights
## @param na.action na.action
## @param offset offset
## @param contrasts contrasts
##' @export
phylo_lmm <- function(formula, data, phylo = NULL, phylonm = NULL, phyloZ = NULL, control, REML = FALSE) {
  phyloZ <- get_phyloZ(phylo, phyloZ, data[[phylonm]])
  lmod <- lFormula(formula = formula, data = data, control = control, REML = REML, phylonm = phylonm, phyloZ = phyloZ)
  devfun <- do.call(mkLmerDevfun, lmod)
  opt <- optimizeLmer(devfun, control = control$optCtrl)
  mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
}

## DON'T document externally
## hacked version of lme4::lFormula
## @name lme4_utils
## @rdname phylo_lmm
## @param formula mixed-model formula
## @param data data frame
## @param REML use restricted max likelihood?
## @param subset subset expression
## @param weights weights
## @param na.action na.action
## @param offset offset
## @param contrasts contrasts
## @param control control
## @param phylonm name of phylogenetic grouping variable
## @param phyloZ phylogenetic Z-matrix
## @param ... additional arguments (ignored)
##' @importFrom lme4 factorize expandDoubleVerts
lFormula <- function(formula, data = NULL, REML = TRUE, subset, weights,
                     na.action, offset, contrasts = NULL, control = lme4::lmerControl(), phylonm, phyloZ, ...) {
  control <- control$checkControl
  mf <- mc <- match.call()
  ignoreArgs <- c("start", "verbose", "devFunOnly", "control")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4:::checkArgs, c(list("lmer"), l...))
  if (!is.null(list(...)[["family"]])) {
    mc[[1]] <- quote(lme4::glFormula)
    if (missing(control)) {
      mc[["control"]] <- glmerControl()
    }
    return(eval(mc, parent.frame()))
  }
  cstr <- "check.formula.LHS"
  lme4:::checkCtrlLevels(cstr, control[[cstr]])
  denv <- lme4:::checkFormulaData(formula, data, checkLHS = control$check.formula.LHS == "stop")
  formula <- as.formula(formula, env = denv)
  lme4:::RHSForm(formula) <- lme4::expandDoubleVerts(lme4:::RHSForm(formula))
  mc$formula <- formula
  m <- match(
    c("data", "subset", "weights", "na.action", "offset"),
    names(mf), 0L
  )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  fr.form <- subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) {
      assign(i, get(i, parent.frame()), environment(fr.form))
    }
  }
  mf$formula <- fr.form
  fr <- eval(mf, parent.frame())
  fr <- factorize(fr.form, fr, char.only = TRUE)
  attr(fr, "formula") <- formula
  attr(fr, "offset") <- mf$offset
  n <- nrow(fr)
  reTrms <- mkReTrms(findbars(lme4:::RHSForm(formula)), fr, phylonm, phyloZ)
  wmsgNlev <- lme4:::checkNlevels(reTrms$flist, n = n, control)
  wmsgZdims <- lme4:::checkZdims(reTrms$Ztlist, n = n, control, allow.n = FALSE)
  if (anyNA(reTrms$Zt)) {
    stop(
      "NA in Z (random-effects model matrix): ", "please use ",
      shQuote("na.action='na.omit'"), " or ", shQuote("na.action='na.exclude'")
    )
  }
  wmsgZrank <- lme4:::checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06)
  fixedform <- formula
  lme4:::RHSForm(fixedform) <- nobars(lme4:::RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(
    fixedfr,
    "terms"
  ), "predvars")
  ranform <- formula
  lme4:::RHSForm(ranform) <- subbars(lme4:::RHSForm(lme4:::reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(
    terms(ranfr),
    "predvars"
  )
  X <- model.matrix(fixedform, fr, contrasts)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
  }
  X <- lme4:::chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- lme4:::checkScaleX(X, kind = scaleX.chk)
  list(
    fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula, phyloZ = phyloZ,
    wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)
  )
}


## hacked version of lme4::mkReTrms
## @rdname lme4_utils
##' @importFrom Matrix KhatriRao fac2sparse sparseMatrix drop0
mkReTrms <- function(bars, fr, phylonm, phyloZ, drop.unused.levels = TRUE) {
  if (!length(bars)) {
    stop("No random effects terms specified in formula",
      call. = FALSE
    )
  }
  stopifnot(is.list(bars), vapply(bars, is.language, NA), inherits(
    fr,
    "data.frame"
  ))
  names(bars) <- lme4:::barnames(bars)
  term.names <- vapply(bars, deparse1, "")
  blist <- lapply(bars, mkBlist, fr, phylonm, phyloZ, drop.unused.levels)
  nl <- vapply(blist, `[[`, 0L, "nl")
  if (any(diff(nl) > 0)) {
    ord <- rev(order(nl))
    blist <- blist[ord]
    nl <- nl[ord]
    term.names <- term.names[ord]
  }
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(rbind, Ztlist)
  names(Ztlist) <- term.names
  q <- nrow(Zt)
  cnms <- lapply(blist, `[[`, "cnms")
  nc <- lengths(cnms)
  nth <- as.integer((nc * (nc + 1)) / 2)
  nb <- nc * nl
  #   if (sum(nb) != q) {
  #     stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)",
  #                  sum(nb), q))
  #   }
  boff <- cumsum(c(0L, nb))
  thoff <- cumsum(c(0L, nth))
  Lambdat <- t(do.call(sparseMatrix, do.call(rbind, lapply(
    seq_along(blist),
    function(i) {
      mm <- matrix(seq_len(nb[i]), ncol = nc[i], byrow = TRUE)
      dd <- diag(nc[i])
      ltri <- lower.tri(dd, diag = TRUE)
      ii <- row(dd)[ltri]
      jj <- col(dd)[ltri]
      data.frame(i = as.vector(mm[, ii]) + boff[i], j = as.vector(mm[
        ,
        jj
      ]) + boff[i], x = as.double(rep.int(
        seq_along(ii),
        rep.int(nl[i], length(ii))
      ) + thoff[i]))
    }
  ))))
  thet <- numeric(sum(nth))
  ll <- list(
    Zt = drop0(Zt), theta = thet, Lind = as.integer(Lambdat@x),
    Gp = unname(c(0L, cumsum(nb)))
  )
  ll$lower <- -Inf * (thet + 1)
  ll$lower[unique(diag(Lambdat))] <- 0
  ll$theta[] <- is.finite(ll$lower)
  Lambdat@x[] <- ll$theta[ll$Lind]
  ll$Lambdat <- Lambdat
  fl <- lapply(blist, `[[`, "ff")
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  }
  else {
    asgn <- seq_along(fl)
  }
  names(fl) <- ufn
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist
  ll
}

## hacked version of mkBlist
## @rdname lme4_utils
## @inheritParams lme4_utils
mkBlist <- function(x, frloc, phylonm, phyloZ, drop.unused.levels = TRUE) {
  frloc <- factorize(x, frloc)
  if (is.null(ff <- tryCatch(eval(substitute(
    lme4:::makeFac(fac),
    list(fac = x[[3]])
  ), frloc), error = function(e) NULL))) {
    stop("couldn't evaluate grouping factor ", deparse(x[[3]]),
      " within model frame:", " try adding grouping factor to data ",
      "frame explicitly if possible",
      call. = FALSE
    )
  }
  if (all(is.na(ff))) {
    stop("Invalid grouping factor specification, ", deparse(x[[3]]),
      call. = FALSE
    )
  }
  if (drop.unused.levels) {
    ff <- factor(ff, exclude = NA)
  }
  if (phylonm[1] %in% names(frloc)) {
    phyloZ <- phyloZ[levels(frloc[, phylonm[1]]), ]
  }
  nl <- length(levels(ff))
  mm <- model.matrix(
    eval(substitute(~foo, list(foo = x[[2]]))),
    frloc
  )
  sm <- fac2sparse(ff, to = "d", drop.unused.levels = drop.unused.levels)
  if (grepl(phylonm[1], x[3])) {
    # nbranch <- ncol(phyloZ)
    nrep <- nrow(sm) / nrow(phyloZ)
    #     bmat <- matrix(0,nrow=nbranch,ncol=length(levels(ff)))
    #     colnames(bmat) <- levels(ff)
    #     for(i in rownames(phyloZ)){
    #       bmat[,which(grepl(i,colnames(bmat)))] <- phyloZ[i,]
    #     }
    lkr <- 1
    rkr <- 1
    if (nrep > 1) {
      if (strsplit(as.character(x[[3]]), ":")[[2]] == phylonm[1]) {
        lkr <- nrep
      }
      if (strsplit(as.character(x[[3]]), ":")[[2]] == phylonm[2]) {
        rkr <- nrep
      }
      # sm <- t(sm)
    }
    ## bmat <- kronecker(kronecker(diag(rkr),phyloZ),diag(lkr))
    bmat <- kronecker(diag(lkr), phyloZ)
    # bmat <- kronecker(phyloZ,diag(lkr))
    # bmat <- kronecker(kronecker(diag(lkr),phyloZ),diag(rkr))  ## This is the one that works?
    #   if(nrow(sm) == nspp*nsite){
    #    sm <- t(sm[ff,])
    #  }
    sm <- t(t(sm) %*% bmat)
    nl <- nrow(sm)
  }
  sm <- KhatriRao(sm, t(mm))
  # dimnames(sm) <- list(rep(1:nrow(sm), each = ncol(mm)), rownames(mm))
  list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
}

#' phylogenetic GLMM
#' @rdname phylo_lmm
#' @param family GLMM family
#' @export
phylo_glmm <- function(formula, data, phylo, phylonm = NULL,
                       phyloZ = NULL, control, family) {
  phyloZ <- get_phyloZ(phylo, phyloZ, data[[phylonm]])
  glmod <- glFormula(formula = formula, data = data, control = control, family = family,
                     phylonm = phylonm, phyloZ = phyloZ)
  # glmod$reTrms <- modify_phylo_retrms(glmod$reTrms,phylo,phylonm,phyloZ)
  devfun <- do.call(mkGlmerDevfun, glmod)
  opt <- optimizeGlmer(devfun)
  devfun <- updateGlmerDevfun(devfun, glmod$reTrms)
  opt <- optimizeGlmer(devfun, stage = 2)
  mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr)
}

## @param orig_fn original function to hack
## @param string_target_start regex for beginning of hacked section
## @param string_target_end regex for end of hacked section
## @param string_replace text to include in place of hacked section
## @param add_args list of additional args (names + values)
## @param after_arg name of argument after which to add new args (not first or last!)
## @param subs_fns named list of functions to :::-preface in order to grab them from other namespaces (name=package, value = character vec of f'n names)
hack_function <- function(orig_fn,
                          string_target_start,
                          string_target_end = string_target_start,
                          string_target_fixed = TRUE,
                          string_replace,
                          add_args,
                          after_arg,
                          subs_fns = NULL) {
  f_tmp <- deparse(body(orig_fn))
  beg_pos <- grep(string_target_start, f_tmp, fixed = string_target_fixed)
  end_pos <- grep(string_target_end, f_tmp, fixed = string_target_fixed)
  f_tmp[beg_pos:end_pos] <- string_replace
  ## allow for private functions/get from namespace
  if (!is.null(subs_fns)) {
    for (i in seq_along(subs_fns)) {
      for (j in seq_along(subs_fns[[i]])) {
        ## substitute
        f0 <- subs_fns[[i]][[j]]
        ## cat(names(subs_fns)[[i]], f0, "\n")
        ## negative lookbehind: don't match f0 preceded by ":::"
        f_tmp <- gsub(sprintf("(?<!:::)%s",f0),
                      sprintf("%s:::%s", names(subs_fns)[[i]], f0), f_tmp,
                      perl = TRUE)
      }
    }
  }
  res <- function() {}
  body(res) <- parse(text = f_tmp)
  ## adjust arguments
  ff <- formals(orig_fn)
  arg_pos <- match(after_arg, names(ff))
  ## FIXME: warn if after_arg == end ...
  formals(res) <- c(ff[1:arg_pos], as.pairlist(add_args), ff[(arg_pos+1):length(ff)])
  environment(res) <- parent.frame()
  return(res)
}

## need this implicitly ...
#' @importFrom lme4 GHrule
glFormula <- hack_function(lme4::glFormula,
                           string_target_start = "mkReTrms(findbars(RHSForm(formula)",
                           string_replace = "  reTrms <- mkReTrms(findbars(lme4:::RHSForm(formula)), fr, phylonm, phyloZ)",
                           add_args = list(phylonm = NULL, phyloZ = NULL),
                           after_arg = "control",
                           subs_fns = list(lme4 = c("checkArgs", "checkCtrlLevels", "checkFormulaData", "checkNlevels", "checkZdims",
                                                    "checkZrank", "RHSForm", "reOnly", "chkRank.drop.cols",
                                                    "checkScaleX")))


## @rdname lme4_utils
## @inheritParams lme4_utils
## glFormula

## glFormula0 <- function(formula, data = NULL, family = gaussian, subset, weights,
##                       na.action, offset, contrasts = NULL, start, mustart, etastart,
##                       control = glmerControl(), phylonm, phyloZ, ...) {
##   control <- control$checkControl
##   mf <- mc <- match.call()
##   if (is.character(family)) {
##     family <- get(family, mode = "function", envir = parent.frame(2))
##   }
##   if (is.function(family)) {
##     family <- family()
##   }
##   if (isTRUE(all.equal(family, gaussian()))) {
##     mc[[1]] <- quote(lme4::lFormula)
##     mc["family"] <- NULL
##     return(eval(mc, parent.frame()))
##   }
##   if (family$family %in% c(
##     "quasibinomial", "quasipoisson",
##     "quasi"
##   )) {
##     stop("\"quasi\" families cannot be used in glmer")
##   }
##   ignoreArgs <- c(
##     "start", "verbose", "devFunOnly", "optimizer",
##     "control", "nAGQ"
##   )
##   l... <- list(...)
##   l... <- l...[!names(l...) %in% ignoreArgs]
##   do.call(lme4:::checkArgs, c(list("glmer"), l...))
##   cstr <- "check.formula.LHS"
##   lme4:::checkCtrlLevels(cstr, control[[cstr]])
##   denv <- lme4:::checkFormulaData(formula, data, checkLHS = control$check.formula.LHS ==
##     "stop")
##   mc$formula <- formula <- as.formula(formula, env = denv)
##   m <- match(c(
##     "data", "subset", "weights", "na.action", "offset",
##     "mustart", "etastart"
##   ), names(mf), 0L)
##   mf <- mf[c(1L, m)]
##   mf$drop.unused.levels <- TRUE
##   mf[[1L]] <- quote(stats::model.frame)
##   fr.form <- subbars(formula)
##   environment(fr.form) <- environment(formula)
##   for (i in c("weights", "offset")) {
##     if (!eval(bquote(missing(x = .(i))))) {
##       assign(i, get(i, parent.frame()), environment(fr.form))
##     }
##   }
##   mf$formula <- fr.form
##   fr <- eval(mf, parent.frame())
##   fr <- factorize(fr.form, fr, char.only = TRUE)
##   attr(fr, "formula") <- formula
##   attr(fr, "offset") <- mf$offset
##   if (!missing(start) && is.list(start)) {
##     attr(fr, "start") <- start$fixef
##   }
##   n <- nrow(fr)
##   reTrms <- mkReTrms(findbars(lme4:::RHSForm(formula)), fr, phylonm, phyloZ)
##   wmsgNlev <- lme4:::checkNlevels(reTrms$flist, n = n, control, allow.n = TRUE)
##   wmsgZdims <- lme4:::checkZdims(reTrms$Ztlist, n = n, control, allow.n = TRUE)
##   wmsgZrank <- lme4:::checkZrank(reTrms$Zt,
##     n = n, control, nonSmall = 1e+06,
##     allow.n = TRUE
##   )
##   fixedform <- formula
##   lme4:::RHSForm(fixedform) <- nobars(lme4:::RHSForm(fixedform))
##   mf$formula <- fixedform
##   fixedfr <- eval(mf, parent.frame())
##   attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(
##     fixedfr,
##     "terms"
##   ), "predvars")
##   ranform <- formula
##   lme4:::RHSForm(ranform) <- subbars(lme4:::RHSForm(lme4:::reOnly(formula)))
##   mf$formula <- ranform
##   ranfr <- eval(mf, parent.frame())
##   attr(attr(fr, "terms"), "predvars.random") <- attr(
##     terms(ranfr),
##     "predvars"
##   )
##   X <- model.matrix(fixedform, fr, contrasts)
##   if (is.null(rankX.chk <- control[["check.rankX"]])) {
##     rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
##   }
##   X <- lme4:::chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
##   if (is.null(scaleX.chk <- control[["check.scaleX"]])) {
##     scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
##   }
##   X <- lme4:::checkScaleX(X, kind = scaleX.chk)
##   list(
##     fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
##     wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)
##   )
## }

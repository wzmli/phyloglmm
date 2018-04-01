#### Helper functions for phyloglmm

phylo.to.Z <- function(r,stand=FALSE){
	ntip <- length(r$tip.label)
	Zid <- Matrix(0.0,ncol=length(r$edge.length),nrow=ntip)
	nodes <- (ntip+1):max(r$edge)
	root <- nodes[!(nodes %in% r$edge[,2])]
	for (i in 1:ntip){
		cn <- i  ## current node
		while (cn != root){
			ce <- which(r$edge[,2]==cn)   ## find current edge
			Zid[i,ce] <- 1   ## set Zid to 1
			cn <- r$edge[ce,1]            ## find previous node
		}
	}
	sig <- det(vcv(r))^(1/ntip)
	Z <- t(sqrt(r$edge.length) * t(Zid))
	if(stand){Z <- t(sqrt(r$edge.length/sig) * t(Zid))}
	rownames(Z) <- r$tip.label
	colnames(Z) <- 1:length(r$edge.length)
	return(Z)                                  
}

## Can me make this simpler?
split_blkMat <- function(M,ind){
	res <- list()
	if (length(ind)==1){
		return(list(M))
		}
	for (i in 1:(length(ind)-1)) {
		v <- (ind[i]+1):ind[i+1]
		res[[i]] <- M[v,v]
	}
	return(res)
}

modify_phylo_retrms <- function(rt,phylo,phylonm,phyloZ,nsp){
	## FIXME: better way to specify phylonm
	## need to replace Zt, Lind, Gp, flist, Ztlist
	n.edge <- nrow(phylo$edge)
	
	## Find the location of phylo random effect
	phylo.pos <- c()
	for(i in 1:length(phylonm)){
	  phylo.pos <- c(phylo.pos,which(names(rt$cnms)==phylonm[[i]]))
	}
	
	## need to know number of number of speices to split index
	## Fixme, the current index spliting is broken
	if(is.null(nsp)){
	nsp <- nrow(rt[["Lambdat"]])/length(rt[["cnms"]][[phylonm]]) 
	}
	inds <- c(0,cumsum(sapply(rt$Ztlist,nrow)))
	## Zt: substitute phylo Z for previous dummy (scalar-intercept) Z
	## Gp: substitute new # random effects (n.edge) for old # (n.phylo)
	Gpdiff <- diff(rt$Gp)  ## old numbers
	Gpdiff_new <- Gpdiff
	rt[["Gp"]] <- as.integer(c(0,cumsum(Gpdiff_new))) ## reconstitute

	## Split index length according to RE complexity w.r.t cov-triangle for each RE
	Lind_split_length <- sapply(rt[["cnms"]]
    , function(i){
      (length(i)*(length(i)+1))/2
    })
  Lind_list <- list()
  ##Fixme, manually hacked Lind_list with a for loop. This will break if we have 
  ## terms on the left side of the bar/pipe
  
	 ### lFormula is creating the all RE index w.r.t nsp lengths into a simple vector, we have to use the function above to split properly
	Lind_list <- split(rt[["Lind"]],rep(seq_along(Lind_split_length),Lind_split_length*nsp))
  if(names(rt[["cnms"]][1]) == "sp:site"){ ### hack
    for(i in 1:length(rt$cnms)){
	    Lind_list[[i]] <- rep(i,sum(rt$Lind==i))
    }
    if(length(rt[["cnms"]][2]$sp) == 2){
      Lind_list[[2]] <- rep(c(2,3,4),sum(rt$Lind==2)) 
    }
  }
  
	## Lambdat: replace block-diagonal element in Lambdat with a
	## larger diagonal matrix
	Lambdat_list <- split_blkMat(rt[["Lambdat"]],inds)
	
	for(i in phylo.pos){
		## each sp is being rep w.r.t the complexity of RE in their respective Zt
		repterms <- length(rt[["cnms"]][[i]])
		if(names(rt[["cnms"]][i]) == "sp:site"){
		  repterms <- nsite ### If this is a special case, it will always be number of sites
		  n.edge <- n.edge*repterms ## This is simply a term to create correct dim for Lind and Lambdat
		}
		## reconstitute Zt from new Ztlist
		## We have to rep the same number of sp terms and edges in phyloZ to match Zt 
		rt[["Ztlist"]][[i]] <- t(kronecker(diag(repterms),phyloZ))%*% rt[["Ztlist"]][[i]]
		## switch places inside kronecker
		Gpdiff_new[i] <- n.edge  ## replace

		## We have to create the Lind to match the theta field with larger reps w.r.t n.edges
		Lind_num <- unique(Lind_list[[i]])
		Lind_list[[i]] <- rep(Lind_list[[i]][seq_along(1:length(Lind_num))],n.edge)
		Lambdat_list[[i]] <- Diagonal(n.edge)

		left_RE_pipe <- length(rt[["cnms"]][[i]])
		
		## function to create sparse matrix templete w.r.t RE complexity
		SM_template <- function(n){
		  rr <- rep(1:n,n:1)
		  cc <- unlist(lapply(seq(n),function(ll){seq(ll,n)}))
		  xx <- as.numeric(rr==cc)
		  return(sparseMatrix(i=rr,j=cc,x=xx))
		}
    temp_lambda <- SM_template(left_RE_pipe)
		Lambdat_list[[i]] <- bdiag(replicate(n.edge,temp_lambda))
	}
	rt[["Zt"]] <- do.call(rbind,rt[["Ztlist"]])
	rt[["Lind"]] <- unlist(Lind_list)
	rt[["Lambdat"]] <- Matrix::.bdiag(Lambdat_list)
	## flist: Not sure how this part is being used.
	rt[["flist"]] <- as.list(rt[["flist"]])
	## Todo: fix flist
# 	for(i in 1:length(rt[["flist"]])){
# 	  rt[["flist"]][i] <- factor(paste0("edge_",seq(n.edge)))
# 	}
	return(rt)
}


phylo_lmm <- function(formula,data,phylo,phylonm,phyloZ,nsp=NULL,control,REML){
	lmod <- lFormula(formula=formula,data = data,control=control, REML=REML)
	lmod$reTrms <- modify_phylo_retrms(lmod$reTrms,phylo,phylonm,phyloZ,nsp)
	devfun <- do.call(mkLmerDevfun, lmod)
	opt <- optimizeLmer(devfun,control=control$optCtrl)
	mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
}

phylo_glmm <- function(formula,data,phylo,phylonm,phyloZ,control,family){
  glmod <- glFormula(formula=formula,data = data,control=control,family)
  glmod$reTrms <- modify_phylo_retrms(glmod$reTrms,phylo,phylonm,phyloZ)
  devfun <- do.call(mkGlmerDevfun, glmod)
  opt <- optimizeGlmer(devfun)
  devfun <- updateGlmerDevfun(devfun,glmod$reTrms)
  opt <- optimizeGlmer(devfun,stage=2)
  mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr)
}
# 
# .simulateFun <- function(object, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA, 
#           ReForm, REForm, REform, newdata = NULL, newparams = NULL, 
#           formula = NULL, family = NULL, weights = NULL, offset = NULL, 
#           allow.new.levels = FALSE, na.action = na.pass, cond.sim = TRUE, 
#           ...){
#   nullWts <- FALSE
#   if (is.null(weights)) {
#     if (is.null(newdata)) 
#       weights <- weights(object)
#     else {
#       nullWts <- TRUE
#       weights <- rep(1, nrow(newdata))
#     }
#   }
#   if (missing(object)) {
#     if (is.null(formula) || is.null(newdata) || is.null(newparams)) {
#       stop("if ", sQuote("object"), " is missing, must specify all of ", 
#            sQuote("formula"), ", ", sQuote("newdata"), ", and ", 
#            sQuote("newparams"))
#     }
#     if (is.character(family)) 
#       family <- get(family, mode = "function", envir = parent.frame())
#     if (is.function(family)) 
#       family <- family()
#     if (is.null(family) || (family$family == "gaussian" && 
#                             family$link == "identity")) {
#       lmod <- lFormula(formula, newdata, weights = weights, 
#                        offset = offset, control = lmerControl(check.formula.LHS = "ignore"))
#       devfun <- do.call(mkLmerDevfun, lmod)
#       object <- mkMerMod(environment(devfun), opt = list(par = NA, 
#                                                          fval = NA, conv = NA), lmod$reTrms, fr = lmod$fr)
#     }
#     else {
#       glmod <- glFormula(formula, newdata, family = family, 
#                          weights = weights, offset = offset, control = glmerControl(check.formula.LHS = "ignore"))
#       devfun <- do.call(mkGlmerDevfun, glmod)
#       object <- mkMerMod(environment(devfun), opt = list(par = NA, 
#                                                          fval = NA, conv = NA), glmod$reTrms, fr = glmod$fr)
#     }
#   }
#   stopifnot((nsim <- as.integer(nsim[1])) > 0, is(object, "merMod"))
#   if (!is.null(newparams)) {
#     object <- setParams(object, newparams)
#   }
#   re.form.miss <- missing(re.form)
#   re.form <- reFormHack(re.form, ReForm, REForm, REform)
#   if (!missing(use.u)) {
#     if (!re.form.miss) {
#       stop("should specify only one of ", sQuote("use.u"), 
#            " and ", sQuote("re.form"))
#     }
#     re.form <- if (use.u) 
#       NULL
#     else ~0
#   }
#   if (is.null(re.form)) {
#     re.form <- noLHSform(formula(object))
#   }
#   if (!is.null(seed)) 
#     set.seed(seed)
#   if (!exists(".Random.seed", envir = .GlobalEnv)) 
#     runif(1)
#   RNGstate <- .Random.seed
#   sigma <- sigma(object)
#   etapred <- predict(object, newdata = newdata, re.form = re.form, 
#                      type = "link", na.action = na.omit, allow.new.levels = allow.new.levels)
#   n <- length(etapred)
#   makeOp <- function(x, y, op = NULL) {
#     if (is.null(op)) {
#       substitute(OP(X), list(X = x, OP = y))
#     }
#     else substitute(OP(X, Y), list(X = x, OP = op, Y = y))
#   }
#   compReForm <- reOnly(formula(object))
#   if (!noReForm(re.form)) {
#     rr <- reOnly(re.form)[[2]]
#     ftemplate <- substitute(. ~ . - XX, list(XX = rr))
#     compReForm <- update.formula(compReForm, ftemplate)[-2]
#   }
#   sim.reff <- if (!is.null(findbars(compReForm))) {
#     if (is.null(newdata) && is.null(re.form)) {
#       newRE <- mkNewReTrms(object, newdata, compReForm, 
#                            na.action = na.action, allow.new.levels = allow.new.levels)
#     }
#     else {
#       cat("using old REs\n")
#       newRE <- getME(object, c("Lambdat", "Zt"))
#     }
#     U <- t(newRE$Lambdat %*% newRE$Zt)
#     u <- rnorm(ncol(U) * nsim)
#     as(U %*% matrix(u, ncol = nsim), "matrix")
#   }
#   else 0
#   val <- if (isLMM(object)) {
#     etapred + sigma * (sim.reff + if (cond.sim) 
#       matrix(rnorm(n * nsim), ncol = nsim)
#       else 0)
#   }
#   else if (isGLMM(object)) {
#     etasim <- etapred + sim.reff
#     family <- normalizeFamilyName(object@resp$family)
#     musim <- family$linkinv(etasim)
#     if (family$family == "binomial" && is.matrix(r <- model.response(object@frame))) {
#       if (nullWts) 
#         weights <- rowSums(r)
#     }
#     if (is.null(sfun <- simfunList[[family$family]]) && is.null(family$simulate)) 
#       stop("simulation not implemented for family", family$family)
#     if (cond.sim) {
#       val <- sfun(object, nsim = 1, ftd = rep_len(musim, 
#                                                   n * nsim), wts = weights)
#     }
#     else {
#       val <- rep_len(musim, n * nsim)
#     }
#     if (family$family == "binomial" && is.matrix(r <- model.response(object@frame))) {
#       lapply(split(val[[1]], gl(nsim, n, 2 * nsim * n)), 
#              matrix, ncol = 2, dimnames = list(NULL, colnames(r)))
#     }
#     else if (family$family == "binomial" && is.factor(val[[1]])) {
#       split(val[[1]], gl(nsim, n))
#     }
#     else split(val, gl(nsim, n))
#   }
#   else stop("simulate method for NLMMs not yet implemented")
#   if (!is.list(val)) {
#     dim(val) <- c(n, nsim)
#     val <- as.data.frame(val)
#   }
#   else class(val) <- "data.frame"
#   names(val) <- paste("sim", seq_len(nsim), sep = "_")
#   f <- fitted(object)
#   nm <- names(f)[!is.na(f)]
#   if (length(nm) == 0) {
#     nm <- as.character(seq(n))
#   }
#   else if (!is.null(newdata)) {
#     nm <- rownames(newdata)
#   }
#   row.names(val) <- nm
#   fit.na.action <- attr(model.frame(object), "na.action")
#   if (!missing(na.action) && !is.null(fit.na.action)) {
#     class.na.action <- class(attr(na.action(NA), "na.action"))
#     if (class.na.action != class(fit.na.action)) {
#       class(fit.na.action) <- class.na.action
#     }
#   }
#   nafun <- function(x) {
#     x[] <- apply(x, 2L, napredict, omit = fit.na.action)
#     x
#   }
#   val <- if (is.matrix(val[[1]])) {
#     structure(lapply(val, nafun), class = "data.frame")
#   }
#   else {
#     as.data.frame(lapply(val, napredict, omit = fit.na.action))
#   }
#   nm2 <- if (is.null(newdata)) 
#     names(napredict(na.omit(f), omit = fit.na.action))
#   else rownames(napredict(newdata, omit = fit.na.action))
#   if (length(nm2) > 0) 
#     row.names(val) <- nm2
#   structure(val, na.action = fit.na.action, seed = RNGstate)
# }

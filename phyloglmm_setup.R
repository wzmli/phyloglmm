#### Helper functions for phyloglmm

phylo.to.Z <- function(r){
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
	Z <- t(sqrt(r$edge.length) * t(Zid))     
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
	# Lind_list <- split(rt[["Lind"]],rep(seq_along(Lind_split_length),Lind_split_length*nsp))
  # if(names(rt[["cnms"]][1]) == "sp:site_name"){ ### hack
  
  for(i in 1:length(rt$cnms)){
	  Lind_list[[i]] <- rep(i,sum(rt$Lind==i))
   }
  # }
	## Lambdat: replace block-diagonal element in Lambdat with a
	## larger diagonal matrix
	Lambdat_list <- split_blkMat(rt[["Lambdat"]],inds)
	
	for(i in phylo.pos){
		## each sp is being rep w.r.t the complexity of RE in their respective Zt
		repterms <- length(rt[["cnms"]][[i]])
		if(names(rt[["cnms"]][i]) == "sp:site"){
		  repterms <- 20 ### Hacked number sites, need to think about how to do this
		  n.edge <- n.edge*repterms ## This is simply a term to create correct dim for Lind and Lambdat
		}
		## reconstitute Zt from new Ztlist
		## We have to rep the same number of sp terms and edges in phyloZ to match Zt 
		rt[["Ztlist"]][[i]] <- t(kronecker(phyloZ,diag(repterms)))%*% rt[["Ztlist"]][[i]]
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


phylo_lmm <- function(formula,data,phylo,phylonm,phyloZ,nsp=NULL,control){
	lmod <- lFormula(formula=formula,data = data,control=control)
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

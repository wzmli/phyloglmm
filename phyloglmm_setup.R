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

modify_phylo_retrms <- function(rt,phylo,phylonm,phyloZ=phylo.to.Z(phylo),sp){
	## FIXME: better way to specify phylonm
	## need to replace Zt, Lind, Gp, flist, Ztlist
	## we have the same number of parameters (theta, lower),
	##  same number of obs
	n.edge <- nrow(phylo$edge)
	phylo.pos <- which(names(rt$cnms)==phylonm)
	inds <- c(0,cumsum(sapply(rt$Ztlist,nrow)))
	## Zt: substitute phylo Z for previous dummy (scalar-intercept) Z
	## Gp: substitute new # random effects (n.edge) for old # (n.phylo)
	Gpdiff <- diff(rt$Gp)  ## old numbers
	Gpdiff_new <- Gpdiff
	rt[["Gp"]] <- as.integer(c(0,cumsum(Gpdiff_new))) ## reconstitute
	## Lind: replace phylo block with the same element, just more values
	Lind_list <- split(rt[["Lind"]],rep(seq_along(Gpdiff),Gpdiff))
	## Lambdat: replace block-diagonal element in Lambdat with a
	## larger diagonal matrix
	Lambdat_list <- split_blkMat(rt[["Lambdat"]],inds)
	for(i in phylo.pos){
		repterms <- nrow(rt[["Ztlist"]][[i]])/length(unique(sp))
		## reconstitute Zt from new Ztlist
		rt[["Ztlist"]][[i]] <- (t(KhatriRao(phyloZ
			, matrix(1,ncol=ncol(phyloZ),nrow=repterms))) 
			%*% rt[["Ztlist"]][[i]]
			)
		Gpdiff_new[i] <- n.edge  ## replace
		Lind_list[[i]] <- rep(Lind_list[[i]][1],n.edge)
		Lambdat_list[[i]] <- (KhatriRao(diag(n.edge)
			, Matrix(1, ncol=n.edge, nrow=repterms))
			)
	}
	rt[["Zt"]] <- do.call(rbind,rt[["Ztlist"]])
	rt[["Lind"]] <- unlist(Lind_list)
	rt[["Lambdat"]] <- Matrix::.bdiag(Lambdat_list)
	## flist: 
	rt[["flist"]] <- as.list(rt[["flist"]])
	rt[["flist"]][[phylonm]] <- factor(paste0("edge_",seq(n.edge)))
	return(rt)
}


phylo_lmm <- function(formula,data,phylo,phylonm,phyloZ,control,sp){
	lmod <- lFormula(formula=formula,data = data,control=control)
	lmod$reTrms <- modify_phylo_retrms(lmod$reTrms,phylo,phylonm,phyloZ,sp)
	devfun <- do.call(mkLmerDevfun, lmod)
	opt <- optimizeLmer(devfun)
	mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
}



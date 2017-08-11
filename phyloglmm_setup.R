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

modify_phylo_retrms <- function(rt,phylo,phylonm,phyloZ,sp,correlated){
	## FIXME: better way to specify phylonm
	## need to replace Zt, Lind, Gp, flist, Ztlist
	n.edge <- nrow(phylo$edge)
	
	## Find the location of phylo random effect
	phylo.pos <- which(names(rt$cnms)==phylonm) 
	
	## need to know number of number of speices to split index
	nsp <- length(unique(sp))
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

	 ### lFormula is creating the all RE index w.r.t nsp lengths into a simple vector, we have to use the function above to split properly
	Lind_list <- split(rt[["Lind"]],rep(seq_along(Lind_split_length),Lind_split_length*nsp))

	## Lambdat: replace block-diagonal element in Lambdat with a
	## larger diagonal matrix
	Lambdat_list <- split_blkMat(rt[["Lambdat"]],inds)
	
	for(i in phylo.pos){
		## each sp is being rep w.r.t the complexity of RE in their respective Zt
		repterms <- nrow(rt[["Ztlist"]][[i]])/length(unique(sp))
		## reconstitute Zt from new Ztlist
		## We have to rep the same number of sp terms and edges in phyloZ to match Zt 
		rt[["Ztlist"]][[i]] <- (t(kronecker(phyloZ,matrix(1,nrow=repterms,ncol=repterms)))
			%*% rt[["Ztlist"]][[i]]
			)
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
		if(left_RE_pipe>1){
			temp_lambda <- SM_template(left_RE_pipe)
			Lambdat_list[[i]] <- bdiag(replicate(n.edge,temp_lambda))
		}
	}
	rt[["Zt"]] <- do.call(rbind,rt[["Ztlist"]])
	rt[["Lind"]] <- unlist(Lind_list)
	rt[["Lambdat"]] <- Matrix::.bdiag(Lambdat_list)
	## flist: 
	rt[["flist"]] <- as.list(rt[["flist"]])
	rt[["flist"]][[phylonm]] <- factor(paste0("edge_",seq(n.edge)))
	return(rt)
}


phylo_lmm <- function(formula,data,phylo,phylonm,phyloZ,control,sp,correlated){
	lmod <- lFormula(formula=formula,data = data,control=control)
	lmod$reTrms <- modify_phylo_retrms(lmod$reTrms,phylo,phylonm,phyloZ,sp,correlated)
	devfun <- do.call(mkLmerDevfun, lmod)
	opt <- optimizeLmer(devfun)
	mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
}

# 
# nb <- ncol(lmod$reTrms$Lambdat)
# lmod$reTrms$Lind <- rep(1:3,nb)
# temp_lambda <- sparseMatrix(i=c(1,1,2),j=c(1,2,2),x=c(1,0,1))
# lmod$reTrms$Lambdat <- bdiag(replicate(nb,temp_lambda))


### simulate phylogenetic tree

library(ape)
library(MASS)
library(sparseMVN)
library(Matrix)
library(dplyr)
## library(MASS)  ## for mvrnorm() ## don't load so we don't screw up dplyr::select()

# tree_seed <- 101
set.seed(tree_seed)

# source("parameters.R")
# nspp <- 10
# nsite <- 5
# nrep <- 1
phy <- rtree(n = nspp)

Vphy <- vcv(phy)

# Generate environmental site variable

X <- rnorm(n=nspp*nsite*nrep, sd=Xsd)

# Generate phylogenetic intercept and slope

phycormat <- matrix(c(1,phyrho.B01,phyrho.B01,1),2,2)
physdvec <- c(physd.B0,physd.B1)
phyvarmat <- physdvec %*% t(physdvec)
phycovmat <- phyvarmat * phycormat


phySigma <- kronecker(phycovmat,Vphy)  ## phylo blocks
# phySigma <- kronecker(Vphy,phycovmat)


b_phy <- MASS::mvrnorm(n=1
	, mu=rep(c(beta0,beta1),each=nspp)
	, Sigma=phySigma
)

Y.phy <- rep(head(b_phy, nspp), each = nrep*nsite) 
## phy_index <- 10000*(as.numeric(factor(Y.phy)) - 1)
phy_index <- 1000000*(as.numeric(factor(rep(seq(nspp), each=nrep*nsite))) - 1)

X.phy <- rep(tail(b_phy, nspp), each = nrep*nsite)*X


# Generate random intercept and slope

cormat <- matrix(c(1, rho.B01, rho.B01, 1), 2, 2)

sdvec <- c(sd.B0, sd.B1)
varmat <- sdvec %*% t(sdvec)
covmat <- varmat * cormat

bSigma <- kronecker(covmat,diag(nspp))

# b <- MASS::mvrnorm(n=1
# 	, mu=rep(c(beta0, beta1), each = nspp)
# 	, Sigma = bSigma)
# 
# Y.re <- rep(head(b, nspp), each = nsite*nrep)
# X.re <- rep(tail(b, nspp), each = nsite*nrep)*X

b <- MASS::mvrnorm(n=nspp
                   , mu=c(beta0, beta1)
                   , Sigma = covmat
                   , empirical = TRUE)

Y.re <- rep(b[,1],each=nsite*nrep)
beta_index <- 1000*(as.numeric(factor(Y.re)) - 1)
X.re <- rep(b[,2],each= nsite*nrep)*X

# Generate random sites and phylogenetic species-site interaction

site <- rep(1:nsite, nspp*nrep)
site_index <- 10000*(as.numeric(factor(site)) - 1)
b_site <- rnorm(n=nsite, mean=beta0, sd = sd.site)
Y.site <- rep(b_site, nspp*nrep)

interaction_sdvec <- rep(sd.interaction, nsite)
interaction_varmat <- interaction_sdvec %*% t(interaction_sdvec)
interaction_covmat <- interaction_varmat * Diagonal(nsite)

interactionSigma <- kronecker(interaction_covmat, Vphy)
# interactionSigma <- kronecker(Vphy,interaction_covmat)

image(Matrix(interactionSigma))

interaction_Cholesky <- Cholesky(interactionSigma)

b_interaction <- sparseMVN::rmvn.sparse(n=1
	, mu=rep(beta0, each = nsite*nspp)
	, CH = interaction_Cholesky
	, prec = FALSE # use covariance_Cholesky when F
)

## interaction_index <- as.numeric(factor(b_interaction)) - 1
interaction_index <- as.numeric(interaction(phy_index,site_index))-1

# Generate observation error

Y.e <- rnorm(nspp*nsite*nrep, sd = sd.resid)

Y <- Y.e #+ Y.site #+ Y.phy + Y.re + X.phy + X.re

dat_nointeraction <- data.frame(sp = rep(rownames(Vphy), each = nrep*nsite)
	, phy_index
	, Y
	, X
	, site
	, site_index)

dat <- (dat_nointeraction
	%>% mutate(sp = factor(sp, levels = rownames(Vphy))
		, obs = sp
		)
	# %>% group_by(sp,site)
	%>% group_by(site,sp)
	%>% mutate(rep=seq(nrep))
	%>% ungroup()
	## this is the order we want if we kronecker(Vphy,interaction_covmat) above:
	##  if we kronecker(interaction_covmat,Vphy) then it should be (site,rep,species) (???)
	## Pez is wrong! 
	## rep pez mistake by screwing up the order (sp,rep,site) != kron(diag,Vphy)
	%>% arrange(sp,rep,site)
	# %>% arrange(site,sp,rep)
	%>% mutate(sp_site = rep(b_interaction, each=nrep)
		, new_y = Y + sp_site
		, interactionindex = rep(interaction_index, each = nrep)
		)
	%>% mutate(index = interactionindex+site_index+phy_index)
)









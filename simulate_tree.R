### simulate phylogenetic tree

library(ape)
library(MASS)
library(sparseMVN)
library(Matrix)
library(dplyr)
## library(MASS)  ## for mvrnorm() ## don't load so we don't screw up dplyr::select()
set.seed(tree_seed)

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
                   , Sigma = covmat)

Y.re <- rep(b[,1],each=nsite*nrep)
X.re <- rep(b[,2],each= nsite*nrep)*X

# Generate random sites and phylogenetic species-site interaction

site <- rep(1:nsite, nspp*nrep)
b_site <- rnorm(n=nsite, mean=beta0, sd = sd.site)
Y.site <- rep(b_site, nspp*nrep)

interaction_sdvec <- rep(sd.interaction, nsite)
interaction_varmat <- interaction_sdvec %*% t(interaction_sdvec)
interaction_covmat <- interaction_varmat * Diagonal(nsite)

# interactionSigma <- kronecker(interaction_covmat, Vphy)
interactionSigma <- kronecker(Vphy,interaction_covmat)


interaction_Cholesky <- Cholesky(interactionSigma)

b_interaction <- sparseMVN::rmvn.sparse(n=1
	, mu=rep(beta0, each = nsite*nspp)
	, CH = interaction_Cholesky
	, prec = FALSE # use covariance_Cholesky when F
)

# Generate observation error

Y.e <- rnorm(nspp*nsite*nrep, sd = sd.resid)

Y <- Y.phy + Y.re + X.phy + X.re + Y.site + Y.e

dat_nointeraction <- data.frame(sp = rep(rownames(Vphy), each = nrep*nsite)
	,Y, X, site)

dat <- (dat_nointeraction
	%>% mutate(sp = factor(sp, levels = rownames(Vphy))
		, obs = sp
		)
	%>% group_by(sp,site)
	%>% mutate(rep=seq(n()))
	%>% ungroup()
	## this is the order we want if we kronecker(Vphy,interaction_covmat) above:
	##  if we kronecker(interaction_covmat,Vphy) then it should be (site,rep,species) (???)
	%>% arrange(sp,rep,site)
	%>% mutate(sp_site = rep(b_interaction, each = nrep)
		, new_y = Y + sp_site
		)
)









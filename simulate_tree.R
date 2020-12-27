### simulate phylogenetic tree

library(ape)
library(MASS)
library(sparseMVN)
library(Matrix)
library(dplyr)
## library(MASS)  ## for mvrnorm() ## don't load so we don't screw up dplyr::select()

quadform <- function(sd_mat,cormat){
  sd_mat %*% cormat %*% sd_mat
}

# tree_seed <- 101
set.seed(tree_seed)

# source("parameters.R")
# nspp <- 10
# nsite <- 5
nrep <- 1
# nspp <- 100
# nsite <- 50
phy <- rtree(n = nspp)

Vphy <- vcv(phy)

# Generate data frame 

interactions <- interaction(1:nsite,1:nspp)
interactions <- interaction(1:nspp,1:nsite)


indexdat <- (data.frame(ints = levels(interactions))
       %>% rowwise()
       %>% mutate(ints = as.character(ints)
                  # , site = unlist(strsplit(ints,"[.]"))[1]
                  # , sp = unlist(strsplit(ints,"[.]"))[2]
                  , site = unlist(strsplit(ints,"[.]"))[2]
                  , sp = unlist(strsplit(ints,"[.]"))[1]
       )  
)

# Species frame: 
## Generate phylogenetic intercept and slope

phycormat <- matrix(c(1,phyrho.B01,phyrho.B01,1),2,2)
physdmat <- diag(c(physd.B0,physd.B1))
phycovmat <- quadform(physdmat,phycormat)

phySigma <- kronecker(phycovmat,Vphy)  ## phylo blocks
# phySigma <- kronecker(Vphy,phycovmat)

b_phy <- MASS::mvrnorm(n=1
	, mu=rep(c(beta0,beta1),each=nspp)
	, Sigma=phySigma
)


# Generate random intercept and slope

cormat <- matrix(c(1, rho.B01, rho.B01, 1), 2, 2)

sdmat <- diag(c(sd.B0, sd.B1))
covmat <- quadform(sdmat,cormat)

bSigma <- kronecker(covmat,diag(nspp))

b <- MASS::mvrnorm(n=nspp
                   , mu=c(beta0, beta1)
                   , Sigma = covmat
                   , empirical = TRUE)

spdat <- data.frame(sp = as.character(1:nspp)
  , b0_phy = head(b_phy,nspp) 
  , b1_phy = tail(b_phy,nspp)
  , b0 = b[,1]
  , b1 = b[,2]
  , tipname = rownames(Vphy)
)


# Site data frame

sitedat <- data.frame(site=as.character(1:nsite)
  , b_site = rnorm(n=nsite,mean=beta0,sd=sd.site)
)

# Generate random sites and phylogenetic species-site interaction

interaction_sdmat <- diag(rep(sd.interaction, nsite))
if(nsite == 1){interaction_sdmat <- matrix(sd.interaction,nrow=1)}
interaction_covmat <- quadform(interaction_sdmat,Diagonal(nsite))

interactionSigma <- kronecker(interaction_covmat, diag(nspp))

interaction_Cholesky <- Cholesky(interactionSigma)

b_interaction <- sparseMVN::rmvn.sparse(n=1
  , mu=rep(beta0, each = nsite*nspp)
  , CH = interaction_Cholesky
  , prec = FALSE # use covariance_Cholesky when F
)


interaction_sdmat <- diag(rep(sd.interaction, nsite))
if(nsite == 1){interaction_sdmat <- matrix(sd.interaction,nrow=1)}
interactionphy_covmat <- quadform(interaction_sdmat,Diagonal(nsite))

interactionphySigma <- Matrix(kronecker(interactionphy_covmat, Vphy),sparse = TRUE)

# interactionSigma <- kronecker(Vphy,interaction_covmat)


interactionphy_Cholesky <- Cholesky(interactionphySigma)

b_interactionphy <- sparseMVN::rmvn.sparse(n=1
	, mu=rep(beta0, each = nsite*nspp)
	, CH = interactionphy_Cholesky
	, prec = FALSE # use covariance_Cholesky when F
)

interactiondat <- data.frame(ints = levels(interactions) ## Check order above, if things don't add up, this is the step to check
  , b_intphy = c(b_interactionphy)    
  , b_int = c(b_interaction)
)

# Generate observation error and environment covariate

dat <- (indexdat[rep(1:nrow(indexdat),each=nrep),]
  %>% left_join(.,spdat)
  %>% left_join(.,sitedat)
  %>% left_join(.,interactiondat)
  %>% rowwise()
  %>% mutate(noise = rnorm(1,sd=sd.resid)
   , X = rnorm(1,sd=Xsd)
   , y_all = b0_phy + b1_phy*X + b0 + b1*X + b_site + b_intphy + noise
   , y = b0_phy + b1_phy*X + b0 + b1*X + b_site + b_intphy + noise
   , y_phy = b0_phy + b1_phy*X + noise
   , y_nophy = b0 + b1*X + noise
   , y_main = b0_phy + b1_phy*X + noise
   , y_interaction = b_int + noise
   , y_interactionphyonly = b_intphy + noise
   , y_interactionphy = b_int + b_intphy + noise
   , sp = tipname)
)

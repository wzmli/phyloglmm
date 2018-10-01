### simulate phylogenetic tree

library(ape)
library(mvtnorm)
library(Matrix)
library(dplyr)
library(lme4)

## library(MASS)  ## for mvrnorm() ## don't load so we don't screw up dplyr::select()
seed <- 104
nsite <- 1
nspp <- 3
source("parameters.R")

phy <- rtree(n = nspp)

## Compute branch lengths using other methods (optional)
#phy <- compute.brlen(phy, method = "Grafen", power = 0.5)

# standardize the phylogenetic covariance matrix to have determinant 1 (optional)
Vphy <- vcv(phy)
#Vphy <- Vphy/max(Vphy)
#Vphy <- Vphy/(det(Vphy)^(1/nspp)) ## for MCMCglmm?

Xsd <- 1

# Generate environmental site variable and standardize it

if(nsite == 1){
  X <- rnorm(nspp,sd=Xsd)
  #	X <- iD %*% rnorm(n=nspp,sd=Xsd)
}
if(nsite > 1){
  X <- rnorm(n = nspp*nsite, sd=Xsd)
  B_site <- matrix(1:nsite, nrow = 1, ncol = nsite)
}
site_name <- sapply(1:nsite, function(x){paste("site","_",x,sep="")})

cormat_phy <- matrix(c(1,rho.B01,rho.B01,1),2,2)
sdvec_phy <- c(sd.B0,sd.B1)
varmat_phy <- sdvec_phy %*% t(sdvec_phy)
covmat_phy <- varmat_phy * cormat_phy


cormat_sp <- matrix(c(1,rho.slopetip,rho.slopetip,1),2,2)
sdvec_sp <- c(sd.tip,sd.slope)
varmat_sp <- sdvec_sp %*% t(sdvec_sp)
covmat_sp <- varmat_sp * cormat_sp


#print(covmat)

Sigma_phy <- kronecker(covmat_phy, Vphy)
Sigma_sp <- kronecker(covmat_sp, diag(nspp))


#if((signal.B0==signal.B1) & (signal.B0 == FALSE)){
#	Sigma <- kronecker(covmat,diag(nspp))
#}

b_phy <- MASS::mvrnorm(n=1
  , mu = rep(c(beta0,beta1), each=nspp)
  , Sigma = Sigma_phy
)

b0_phy <- b_phy[1:nspp] 
b1_phy <- b_phy[(nspp+1):(2*nspp)]

b_sp <- MASS::mvrnorm(n=1
  , mu = rep(c(beta0,beta1), each=nspp)
  , Sigma = Sigma_sp
)

b0_sp <- b_sp[1:nspp] 
b1_sp <- b_sp[(nspp+1):(2*nspp)]

#y <- (matrix(outer(b0_phy, array(1, dim = c(1, nsite))), nrow = nspp, ncol = nsite) 
#  + matrix(outer(b1_phy, X), nrow = nspp, ncol = nsite)
#  + (matrix(outer(b0_sp, array(1, dim = c(1, nsite))), nrow = nspp, ncol = nsite))
#  + matrix(outer(b1_sp, X), nrow = nspp, ncol = nsite)
#)

Y <- (rep(b0_phy,nsite) 
+	rep(b0_sp,nsite)
#	+ rep(b1_sp,nsite)*X
#	+ rep(b1_phy,nsite)*X 
#	+ rnorm(nspp*nsite, sd=sd.resid) 
)

dat <- data.frame(X, Y, sp=rep(phy$tip.label,nsite))

#if(nsite == 1){
#  y <- matrix(outer(b0_phy, array(1, dim = c(1, nsite))), nrow = nspp,
#              ncol = nsite) + matrix(b1_phy*X, nrow = nspp, ncol = nsite)
#}
#e <- rnorm(nspp * nsite, sd = sd.resid) # add residual variance 
#y <- y + matrix(e, nrow = nspp, ncol = nsite)
#y <- matrix(y, nrow = nspp * nsite, ncol = 1)

#V <- kronecker(diag(nrow=nsite, ncol=nsite), Vphy)
#Y.ss <- rmvnorm(n=1, sigma=ss^2*V)

#Y <- y + t(Y.ss)
#Y <- matrix(y, nrow = nspp, ncol = nsite)

# name the simulated species 1:nspp and sites 1:nsites
#rownames(Y) <- 1:nspp
#colnames(Y) <- 1:nsite


# Transform data matrices into "long" form, and generate a data frame
#YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)

#XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow =
#               nspp * nsite, ncol = 1)

# if(nsite == 1){
#XX <- matrix(X,nrow=nspp,ncol=1)
# }

#site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol =
#                                           1)), nrow = nspp * nsite, ncol = 1)

#sp <- rep(phy$tip.label,nsite)

#dat <- data.frame(Y = YY, X=XX, site = as.factor(site), sp = as.factor(sp),site_name = rep(site_name,each=nspp))

#print(dim(dat))
#print(head(dat))

source("new_phylo_setup.R")

phyZ <- phylo.to.Z(phy,stand=FALSE)

dat2 <- (dat
        %>% mutate(obs = sp
#                   , site = site_name
                   , y_na = NA)
)	


tempmod <- phylo_lmm(Y ~ X
#                     + (1 | sp:site)
#+ (1 | sp)
+ (1 + X | sp)
#	+ (1 | obs)
  + (1 + X | obs)
# + (1 | site)
, data=dat2
, phylonm = c("sp","sp:site")
, phylo = phy
, phyloZ=phyZ
# , nsp = nsite
, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
, REML = TRUE
)

print(dat)
print(plot(phy))

print(summary(tempmod))

#rdnosave()

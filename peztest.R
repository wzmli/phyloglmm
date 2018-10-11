library(ape)
library(Matrix)
library(lme4)
library(dplyr)
library(pez)

seed <- 101
source("parameters.R")
nspp <- 30
nsite <- 20
rho.B01 <- 0.7
source("simulate_tree.R")
source("new_phylo_setup.R")

## pez and lme4 setup ----

sp.int <- list(1, sp = dat$sp, covar = diag(nspp))

# random intercept with species showing phylogenetic covariances
phy.int <- list(1, sp = dat$sp, covar = Vphy)

# random slope with species independent
sp.X <- list(dat$X, sp = dat$sp, covar = diag(nspp))

# random slope with species showing phylogenetic covariances
phy.X <- list(dat$X, sp = dat$sp, covar = Vphy)

phyZ <- phylo.to.Z(phy,stand=FALSE)

dat <- (dat %>% mutate(obs = sp))

## test 1 random int (good, all match up)----

pez1 <- communityPGLMM(Y~X
  , data = dat
  , sp = dat$sp
  , site = dat$site
  , random.effects = list(sp.int)
  , REML = TRUE
  , verbose = FALSE
)

lme41 <- lmer(Y ~ X + (1|sp), data=dat)

phylolmm1 <- phylo_lmm(Y ~ X + (1|obs)
  , data=dat
  , phylonm = c("sp","site:sp")
  , phylo = phy
  , phyloZ=phyZ
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
  )
  , REML = TRUE
)

print(summary(pez1))
print(summary(lme41))
print(summary(phylolmm1))


### test 2 random slope (pez slope optim control is different?) ----

pez2 <- communityPGLMM(Y~X
                       , data = dat
                       , sp = dat$sp
                       , site = dat$site
                       , random.effects = list(sp.int, sp.X)
                       , REML = TRUE
                       , verbose = FALSE
)

lme42 <- lmer(Y ~ X + (1|sp) + (0 + X | sp), data=dat)

phylolmm2 <- phylo_lmm(Y ~ X + (1|obs) + (0 + X|obs)
                       , data=dat
                       , phylonm = c("sp","site:sp")
                       , phylo = phy
                       , phyloZ=phyZ
                       , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
                       )
                       , REML = TRUE
)

print(summary(pez2))
print(summary(lme42))
print(summary(phylolmm2))


### test 3 random phylo intercept (different results) ----


pez3 <- communityPGLMM(Y~X
                       , data = dat
                       , sp = dat$sp
                       , site = dat$site
                       , random.effects = list(phy.int)
                       , REML = TRUE
                       , verbose = FALSE
)

phylolmm3 <- phylo_lmm(Y ~ X + (1|sp)
                       , data=dat
                       , phylonm = c("sp","site:sp")
                       , phylo = phy
                       , phyloZ=phyZ
                       , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
                       )
                       , REML = TRUE
)

print(summary(pez3))
print(summary(phylolmm3))

### use ives tree ----

set.seed(seed)

phy <- rcoal(n = nspp)

## Compute branch lengths using other methods (optional)
# phy <- compute.brlen(phy, method = "Grafen", power = 0.5)

# standardize the phylogenetic covariance matrix to have determinant 1 (optional)
Vphy <- vcv(phy)
Vphy <- Vphy/(det(Vphy)^(nspp)) ## for MCMCglmm?

iD <- t(chol(Vphy))
Xsd <- 1

# Generate environmental site variable and standardize it

if(nsite == 1){
  X <- rnorm(n=nspp,sd=Xsd)
  #	X <- iD %*% rnorm(n=nspp,sd=Xsd)
}
if(nsite > 1){
  X <- rnorm(n = nspp*nsite, sd=Xsd)
  B_site <- matrix(1:nsite, nrow = 1, ncol = nsite)
  # X <- (X - mean(X))/sd(X)
}
site_name <- sapply(1:nsite, function(x){paste("site","_",x,sep="")})

cormat <- matrix(c(1,rho.B01,rho.B01,1),2,2)
sdvec <- c(sd.B0,sd.B1)
varmat <- sdvec %*% t(sdvec)
covmat <- varmat * cormat

#print(covmat)

Sigma <- kronecker(covmat,Vphy)

#if((signal.B0==signal.B1) & (signal.B0 == FALSE)){
#	Sigma <- kronecker(covmat,diag(nspp))
#}

b.all <- MASS::mvrnorm(n=1
                       , mu=rep(c(beta0,beta1),each=nspp)
                       , Sigma=Sigma
)

b0 <- b.all[1:nspp] 
b1 <- b.all[(nspp+1):(2*nspp)]

y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp,
            ncol = nsite) + matrix(outer(b1, X), nrow = nspp, ncol = nsite)

if(nsite == 1){
  y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp,
              ncol = nsite) + matrix(b1*X, nrow = nspp, ncol = nsite)
}
e <- rnorm(nspp * nsite, sd = sd.resid) # add residual variance 
y <- y + matrix(e, nrow = nspp, ncol = nsite)
y <- matrix(y, nrow = nspp * nsite, ncol = 1)

Y <- y
Y <- matrix(Y, nrow = nspp, ncol = nsite)

# name the simulated species 1:nspp and sites 1:nsites
rownames(Y) <- 1:nspp
colnames(Y) <- 1:nsite


# Transform data matrices into "long" form, and generate a data frame
YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)

XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow =
               nspp * nsite, ncol = 1)

# if(nsite == 1){
XX <- matrix(X,nrow=nspp,ncol=1)
# }

site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol =
                                           1)), nrow = nspp * nsite, ncol = 1)

sp <- rep(phy$tip.label,nsite)

dat <- data.frame(Y = YY, X=XX, site = as.factor(site), sp = as.factor(sp),site_name = rep(site_name,each=nspp))

### setup for new tree ----
phy.int <- list(1, sp = dat$sp, covar = Vphy)


### test 4 testing ives tree (with no standardizing phyZ) still doesn't match

pez4 <- communityPGLMM(Y~X
                       , data = dat
                       , sp = dat$sp
                       , site = dat$site
                       , random.effects = list(phy.int)
                       , REML = TRUE
                       , verbose = FALSE
)

phyZ <- phylo.to.Z(phy,stand=FALSE)

phylolmm4 <- phylo_lmm(Y ~ X + (1|sp)
                       , data=dat
                       , phylonm = c("sp","site:sp")
                       , phylo = phy
                       , phyloZ=phyZ
                       , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
                       )
                       , REML = TRUE
)

print(summary(pez4))
print(summary(phylolmm4))


### test 5 testing ives tree (with standardizing phyZ)

phyZ <- phylo.to.Z(phy,stand=TRUE)

phylolmm4stand <- phylo_lmm(Y ~ X + (1|sp)
                       , data=dat
                       , phylonm = c("sp","site:sp")
                       , phylo = phy
                       , phyloZ=phyZ
                       , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
                       )
                       , REML = TRUE
)

print(summary(pez4))
print(summary(phylolmm4stand))


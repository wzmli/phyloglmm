library(pez)
library(ape)
library(MASS)
library(mvtnorm)
library(dplyr)

# Simulate a phylogeny that has a lot of phylogenetic signal (power = 1.3)

seed <- 101

source("parameters.R")
nspp <- 25
nsite <- 20

phy <- rtree(n = nspp)

#phy <- compute.brlen(rtree(n = nspp), method = "Grafen", power = 1.3)

# Simulate species means
mean.sp <- rTraitCont(phy, model = "BM", sigma=sd.B0^2)
# Replicate values of mean.sp over sites
Y.sp <- rep(mean.sp, times=nsite)

# Simulate site means
mean.site <- rnorm(nsite, sd=sd.site)
# Replicate values of mean.site over sp
Y.site <- rep(mean.site, each=nspp)

# Compute the covariance matrix for phylogenetic attraction.
Vphy <- vcv(phy)

# Standardize Vphy to have determinant = 1.
#Vphy <- Vphy/(det(Vphy)^(1/nspp))

# Construct the overall covariance matrix for phylogenetic attraction.
V <- kronecker(diag(nrow=nsite, ncol=nsite), Vphy)
Y.ss <- array(t(rmvnorm(n=1, sigma=ss^2*V)))

# Simulate residual errors
Y.resid <- rnorm(nspp*nsite, sd=sd.resid)

# Construct the dataset
d <- data.frame(sp=rep(phy$tip.label, times=nsite), site=rep(1:nsite, each=nspp))

# Simulate abundance data
d$Y <- Y.sp + Y.site + Y.ss + Y.resid

dat <- d %>% mutate(sp = factor(sp), site = factor(site))

print(head(d))

sp.int <- list(1, sp = d$sp, covar = diag(nspp))

# random intercept with species showing phylogenetic covariances
phy.int <- list(1, sp = dat$sp, covar = Vphy)
# sp:site
phy.interaction <- list(1, sp = dat$sp, covar = Vphy, site = dat$site)

site.int <- list(1, site=dat$site, covar = diag(nsite))

# Analyze the model
pezfit <- communityPGLMM(Y ~ 1
                         , family = "gaussian"
                         , sp = dat$sp
                         , site = dat$site
                         , random.effects = list(phy.interaction
                                                 , phy.int
                                                 , sp.int
                                                 , site.int
                         ),	data=dat)

print(summary(pezfit))

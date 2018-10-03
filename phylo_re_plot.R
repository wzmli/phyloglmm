## What are phylogenetic signals?
library(ape)
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtree)
library(cowplot)
library(Hmisc)

tree_seed <- 101

set.seed(tree_seed)
nspp <- 200
nrep <- 50
phy <- rtree(n = nspp)

x <- rnorm(nspp*nrep, sd=1)

gg_tree <- (ggtree(phy)
  + geom_tiplab()
)

# print(plot(phy))

# print(phy$tip.label)

Vphy <- vcv(phy)

physd.y <- 10
physd.x <- 0

sd.y <- 5
sd.x <- 0 

sd.resid <- 2

phycormat <- matrix(c(1,0,0,1),nrow=2)
physdvec <- c(physd.y, physd.x)
phyvarmat <- physdvec %*% t(physdvec)
phycovmat <- phyvarmat * phycormat

physigma <- kronecker(phycovmat, Vphy) 

betas <- c(0,0)

seed <- 1011

set.seed(seed)

b_phy <- MASS::mvrnorm(n=1, mu=rep(betas,each=nspp), Sigma=physigma)

Y.phy <- rep(head(b_phy,nspp), each = nrep)
X.phy <- rep(tail(b_phy,nspp), each = nrep) 


cormat <- matrix(c(1,0,0,1),nrow=2)

print(cormat)
sdvec <- c(sd.y, sd.x)
varmat <- sdvec %*% t(sdvec)
covmat <- varmat * cormat


sigma_b <- kronecker(covmat, diag(nspp)) 

b <- MASS::mvrnorm(n=1, mu=rep(betas,each=nspp), Sigma=sigma_b)

Y.int <- rep(head(b,nspp), each = nrep)

X.int <- rep(tail(b,nspp), each = nrep)


Y.e <- rnorm(nspp*nrep, sd=sd.resid)
X.e <- rnorm(nspp*nrep, sd=1)

Y.phyint <- Y.phy + Y.int
X.phyint <- X.phy + X.int
X.phyint.e <- X.phyint + X.e
Y.phyint.e <- Y.phyint + Y.e
Y <- Y.phyint.e + X.phyint.e

dat <- data.frame(sp = rep(rownames(Vphy), each = nrep)
  , Y.phy
  , Y.phyint
  , Y.phyint.e
  , Y
  , X.phy
  , X.phyint.e
)


dat <- dat %>% mutate(sp = factor(sp,levels=rownames(Vphy)), obs=sp, site=rep(1:nrep,nspp))

source("new_phylo_setup.R")

phyZ <- phylo.to.Z(phy, stand=FALSE)

mod1 <- phylo_lmm(Y.phyint.e~1
  + (1 | sp)
  + (1 | obs)
  , data=dat
  , phylonm = c("sp","sp:site")
  , phylo = phy
  , phyloZ=phyZ
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  , REML = FALSE                  
)

source("phyr_hacked.R")

mod2 <- phyr_hacked(Y.phyint.e ~ 1 + (1|sp__)
  , data = dat
  , family = "gaussian"
  , tree = phy
  , REML = FALSE
  , cpp = TRUE
)

mod3 <- phyr::communityPGLMM(Y.phyint.e ~ 1 + (1|sp__)
  , data = dat
  , family = "gaussian"
  , tree = phy
  , REML = FALSE
  , cpp = TRUE
)


print(summary(mod1))
print(mod2)
print(mod3)

## What are phylogenetic signals?
library(ape)
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtree)
library(cowplot)
library(Hmisc)

tree_seed <- 121

set.seed(tree_seed)
nspp <- 1200
nsite <- 10
nrep <- 10

phy <- rtree(n = nspp)

X <- rnorm(nspp*nrep*nsite, sd=1)

gg_tree <- (ggtree(phy)
  + geom_tiplab()
)

# print(plot(phy))

# print(phy$tip.label)

Vphy <- vcv(phy)

physd.y <- 10
physd.x <- 10

sd.y <- 5
sd.x <- 2 

sd.site <- 2
sd.resid <- 1

phycormat <- matrix(c(1,0.7,0.7,1),nrow=2)
physdvec <- c(physd.y, physd.x)
phyvarmat <- physdvec %*% t(physdvec)
phycovmat <- phyvarmat * phycormat

physigma <- kronecker(phycovmat, Vphy) 

betas <- c(0,0)

seed <- 192409

set.seed(seed)

b_phy <- MASS::mvrnorm(n=1, mu=rep(betas,each=nspp), Sigma=physigma)

Y.phy <- rep(head(b_phy,nspp), each = nrep*nsite)
X.phy <- rep(tail(b_phy,nspp), each = nrep*nsite)*X 

## each trait will have nrep*nsite repeats

cormat <- matrix(c(1,0.4,0.4,1),nrow=2)


print(cormat)
sdvec <- c(sd.y, sd.x)
varmat <- sdvec %*% t(sdvec)
covmat <- varmat * cormat


sigma_b <- kronecker(covmat, diag(nspp)) 

b <- MASS::mvrnorm(n=1, mu=rep(betas,each=nspp), Sigma=sigma_b)


Y.int <- rep(head(b,nspp), each = nrep*nsite)

X.int <- rep(tail(b,nspp), each = nrep*nsite)*X

b_site <- rnorm(n=nsite, mean=0, sd=sd.site)

Y.site <- rep(rep(b_site,each=1),nspp*nrep)

Y.e <- rnorm(nspp*nrep*nsite, sd=sd.resid)

Y.phyint <- Y.phy + Y.int
X.phyint <- X.phy + X.int
Y.phyint.site <- Y.phyint + Y.site
Y.phyint.site.e <- Y.phyint + Y.site + Y.e
Y <- Y.phyint.site.e + X.phyint

dat <- data.frame(sp = rep(rownames(Vphy), each = nrep*nsite)
  , X
  , Y.phy
  , Y.phyint
  , Y.phyint.site
  , Y.phyint.site.e
  , Y
  , site = rep(rep(1:nsite,each=1),nspp*nrep)
  # , Y
  # , X.phy
  # , X.phyint.e
)


dat <- dat %>% mutate(sp = factor(sp,levels=rownames(Vphy)), obs=sp)

source("new_phylo_setup.R")

phyZ <- phylo.to.Z(phy, stand=FALSE)

mod1 <- phylo_lmm(Y~1+X
  + (1 + X| obs)
  + (1 + X| sp)
  + (1 | site)
  , data=dat
  , phylonm = c("sp","sp:site")
  , phylo = phy
  , phyloZ=phyZ
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  , REML = FALSE                  
)

source("phyr_hacked.R")

# mod2 <- phyr_hacked(Y.phyint.site.e ~ 1 + (1|sp__) + (1|site)
#   , data = dat
#   , family = "gaussian"
#   , tree = phy
#   , REML = FALSE
#   , cpp = TRUE
# )

# mod3 <- phyr::communityPGLMM(Y.phyint.e ~ 1 + (1|sp__)
#   , data = dat
#   , family = "gaussian"
#   , tree = phy
#   , REML = FALSE
#   , cpp = TRUE
# )


print(summary(mod1))
# print(mod2)
# print(mod3)
# View(dat %>% filter(sp == "t1") )

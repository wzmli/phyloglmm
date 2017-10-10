## simulate and fit compound symmetric case
library(ape)
library(Matrix)
library(lme4)
library(dplyr)
library(pez)


phyZ <- phylo.to.Z(phy)
# 
# debug(phylo_lmm)
# debug(modify_phylo_retrms)


tempmod <- phylo_lmm(Y ~ 1
  + (1|sp)
  + (1|sp:site)
	, data=dat
	, phylonm = c("sp","sp:site")
	, phylo = phy
	, phyloZ=phyZ
	, nsp = nsite
	, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
)

new_y <- simulate(tempmod, newparams=list(theta=c(6/sigma(tempmod),2/sigma(tempmod)),beta=0))
dat$new_y <- new_y[[1]]

fitmod <- phylo_lmm(new_y ~ 1
  + (1|sp)
  + (1|sp:site)
  , data=dat
  , phylonm = c("sp","sp:site")
  , phylo = phy
  , phyloZ=phyZ
  , nsp = nsite
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
)

print(summary(fitmod))
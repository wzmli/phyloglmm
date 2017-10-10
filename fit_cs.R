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

# random intercept with species independent
re.1 <- list(1, sp = dat$sp, covar = diag(nspp))

# random intercept with species showing phylogenetic covariances
re.2 <- list(1, sp = dat$sp, covar = Vphy)

# random slope with species independent
re.3 <- list(dat$X, sp = dat$sp, covar = diag(nspp))

# random slope with species showing phylogenetic covariances
re.4 <- list(dat$X, sp = dat$sp, covar = Vphy)

# sp:site
re.nested.phy <- list(1, sp = dat$sp, covar = Vphy, site = dat$site)



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

if(platform == "pez"){
  fittime <- system.time(
    fitmod <- communityPGLMM(new_y ~ 1
        , data = dat
        , family = "gaussian"
        , sp = dat$sp
        , site = dat$site
        , random.effects = list(re.2, re.nested.phy)
        , REML = FALSE
        , verbose = FALSE
    )
  )
}


if(platform == "lme4"){
  fittime <- system.time(
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
  )
}

print(fittime)
print(summary(fitmod))

saveRDS(list(fittime,fitmod),file=paste(datadir,paste("cs",platform,seed,"Rds",sep="."),sep=""))

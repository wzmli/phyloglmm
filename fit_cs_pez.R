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


tempmod <- phylo_lmm(Y ~ X
                     #  + (1|sp)
                     + (1 + X|sp)
                     + (1|sp:site)
                     , data=dat
                     , phylonm = c("sp","sp:site")
                     , phylo = phy
                     , phyloZ=phyZ
                     , nsp = nsite
                     , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
)

t1 <- sd.B0/sigma(tempmod)
t2 <- rho.B01*sd.B1/sigma(tempmod)
t3 <- sqrt((sd.B1/sigma(tempmod))^2 - t2^2)

new_y <- simulate(tempmod, newparams=list(theta=c(ss/sigma(tempmod)
                                                  , t1
                                                  , t2
                                                  , t3)
                                          , beta = c(beta0,beta1)))

dat$new_y <- new_y[[1]]

fittime <- system.time(
  fitpez <- communityPGLMM(new_y ~ X
        , data = dat
        , family = "gaussian"
        , sp = dat$sp
        , site = dat$site
        , random.effects = list(re.2, re.4, re.nested.phy)
        , REML = FALSE
        , verbose = FALSE
    )
  )

print(fittime)
print(summary(fitpez))

saveRDS(list(fittime,fitpez),file=paste(datadir,paste("pez_cs",size,seed,"rds",sep="."),sep=""))

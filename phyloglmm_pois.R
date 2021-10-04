## poisson fit
#### Fitting phyloglmm via lme4

load('.simulate_poistree.RData')
library(phyloglmm)

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

print(covmat)
print(cov2cor(covmat))

phyZ <- phylo.to.Z(phy)

dat <- (dat
  %>% mutate(obs = sp)
)

# lme4timecor <- system.time(
  lme4fitcor <- phylo_glmm(Y ~ X + (1|sp) + (1|obs)
    , data=dat
    , phylonm = "sp"
    , family = poisson
    , phylo = phy
    , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  )
# )

# print(lme4timecor)
print(summary(lme4fitcor))

library(MCMCglmm)

nitt <- 5e3 ## was 5e6
inv.phylo <- inverseA(phy,nodes="TIPS",scale=FALSE)
prior <- list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
MCMC_time <- system.time(
  MCMCglmm_fit <- MCMCglmm(Y~X,random=~sp,
                           family="poisson",ginverse=list(sp=inv.phylo$Ainv),
                           prior=prior,data=dat,nitt=nitt,burnin=100,
                           thin=nitt/100,verbose=TRUE))

ss <- summary(MCMCglmm_fit)



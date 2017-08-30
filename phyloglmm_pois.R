## poisson fit
#### Fitting phyloglmm via lme4

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

lme4timecor <- system.time(
  lme4fitcor <- phylo_glmm(Y ~ X + (1+X|sp)
    , data=dat
    , phylonm = "sp"
    , family = poisson
    , phylo = phy
    , phyloZ=phyZ
    , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  )
)

print(lme4timecor)
print(summary(lme4fitcor))


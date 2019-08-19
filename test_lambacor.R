library(ape)
library(Matrix)
library(lme4)
library(dplyr)

phyZ <- phylo.to.Z(phy,stand=FALSE)

dat <- (dat
        %>% mutate(obs = sp)
)	


lme4fit <- phylo_lmm(Y ~ X + (1|sp)
                     , data=dat
                     , phylonm = c("sp","site:sp")
                     , phylo = phy
                     , phyloZ=phyZ
                     , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
                     )
                     , REML = FALSE
                     , lambhack = FALSE
)
lme4fit
aa <- lme4fit
aa

# debug(mkReTrms)
debug(mkLmerDevfun)
debug(phylo_lmm)
lme4fit2 <- phylo_lmm(Y ~ X + (1|obs)
                     , data=dat
                     , phylonm = c("sp","site:sp")
                     , phylo = phy
                     , phyloZ=phyZ
                     , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
                     )
                     , REML = FALSE
                     , lambhack = TRUE
)
lme4fit2
bb <- lme4fit2


### Fitting phyloglmm with TMB

library(dplyr)
library(glmmTMB)
library(lme4)
library(Matrix)

phyZ <- phylo.to.Z(phy,stand=TRUE)

# dat <- (dat
#         %>% rowwise()
#         %>% mutate(phylo=paste("t",sp,sep="")
#                    , obs=phylo
#         )
# )
dat <- (dat
        %>% mutate(obs = sp)
)	

dat <- data.frame(dat)

# source("glmmTMBhacked.R")
# debug(glmmTMBhacked)
# debug(mktempmodhacked)
# debug(mkTMBStruchacked)
# debug(getXReTrmshacked)
# debug(mkReTrms)
hackedmod <- glmmTMBhacked(Y ~ X  + (1 + X|sp)
  , data=dat
  , phyloZ = phyZ
  , phylonm = "sp"
  , doFit=TRUE
  , dispformula = ~1
  , control=glmmTMBControl(optCtrl=list(trace=1))
  , REML = TRUE
  ) 


print(summary(hackedmod))

#rdnosave()

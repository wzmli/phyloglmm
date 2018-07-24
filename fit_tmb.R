### Fitting phyloglmm with TMB

library(dplyr)
library(glmmTMB)
library(lme4)
library(Matrix)

phyZ <- phylo.to.Z(phy)

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

source("glmmTMBhacked.R")
source("new_phylo_setup.R")
# debug(glmmTMBhacked)
# debug(mktempmodhacked)
# debug(getXReTrmshacked)
# debug(mkReTrms)
hackedmod <- glmmTMBhacked(Y ~ X  + (1|sp)
  , data=dat
  , phyloZ = phyZ
  , phylonm = "sp"
  , doFit=FALSE
  ) # doFit=FALSE) in BB's update

tempmod <- glmmTMB(Y ~ X  + (1|sp)
  , data=dat
  , doFit=FALSE
)

n.edge <- ncol(phyZ)

tempmod$condReStruc$`1 | sp`$blockReps <- n.edge
tempmod$condList$Z <- t(hackedmod$condList$reTrms$Zt)
## data *inside* data.tmb is actually the most critical to allow correct fit
tempmod$data.tmb$terms$`1 | sp`$blockReps <- n.edge
tempmod$data.tmb$Z <- t(hackedmod$condList$reTrms$Zt)
tempmod$parameters$b <- rep(0,ncol(hackedmod$data.tmb$Z))

ff <- glmmTMB:::fitTMB(tempmod)

print(ff)
print(tempmod)

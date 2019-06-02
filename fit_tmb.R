### Fitting phyloglmm with TMB

library(dplyr)
library(glmmTMB)
library(lme4)
library(Matrix)
t1 <- proc.time()

phyZ <- phylo.to.Z(phy,stand=FALSE)

dat <- (dat
  %>% mutate(obs = sp
      )
)	

if(numsite == "ss"){
  hackedmod <- glmmTMBhacked(Y ~ X  
    # + (1|sp) + (0 + X| sp) 
    + (1 + X|sp)
    , data=dat
    , phyloZ = phyZ
    , phylonm = c("sp", "sp:site")
    , doFit=TRUE
    # , dispformula = ~1
    , control=glmmTMBControl(optCtrl=list(trace=1,iter.max=1e5,eval.max=1e5))
    , REML = FALSE
  ) 
}

# source("glmmTMBhacked.R")
# debug(glmmTMBhacked)
# debug(mkTMBStruchacked)
# debug(getXReTrmshacked)
# debug(mkReTrms)

if(numsite == "ms"){
	hackedmod <- glmmTMBhacked(new_y ~ X  
+ (1 | sp:site)
#	+ (1 | sp)
	+ (1 + X | sp)
#	+ (1 | obs)
	+ (1 + X | obs)
	+ (1 | site)
  , data=dat
  , phyloZ = phyZ
  , phylonm = c("sp", "sp:site")
  , doFit=TRUE
  , dispformula = ~1
  , control=glmmTMBControl(optCtrl=list(trace=1,iter.max=1e4,eval.max=1e4))
  , REML = FALSE
  ) 
}
t2 <- proc.time()
glmmTMBtime <- t2-t1
print(summary(hackedmod))

glmmTMB_list <- list(hackedmod, glmmTMBtime)


print(summary(hackedmod))
print(hackedmod$fit$convergence)

saveRDS(glmmTMB_list, file=paste("datadir/glmmTMB/glmmTMB",numsite,size,tree_seed,"rds",sep="."))

#rdnosave()


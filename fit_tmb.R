### Fitting phyloglmm with TMB

library(dplyr)
library(glmmTMB)
library(lme4)
library(Matrix)
t1 <- proc.time()

phyZ <- phylo.to.Z(phy,stand=TRUE)

dat <- (dat
  %>% mutate(obs = sp
      , site = site_name
      )
)	

t2 <- proc.time()

if(numsite == "ss"){
  hackedmod <- glmmTMBhacked(Y ~ X  
    + (1| sp) 
    , data=dat
    , phyloZ = phyZ
    , phylonm = c("sp", "sp:site")
    , doFit=TRUE
    , dispformula = ~0
    # , control=glmmTMBControl(optCtrl=list(trace=1))
    , REML = TRUE
  ) 
}

t3 <- proc.time()
glmmTMBtime <- t3-t1

# source("glmmTMBhacked.R")
# debug(glmmTMBhacked)
# debug(mkTMBStruchacked)
# debug(getXReTrmshacked)
# debug(mkReTrms)

if(numsite == "ms"){
  tempmod <- phylo_lmm(Y ~ X
    + (1 | sp:site)
    + (1 + X |sp)
    + (1 + X | obs)
    + (1 | site)
    , data=dat
    , phylonm = c("sp","sp:site")
    , phylo = phy
    , phyloZ=phyZ
    , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
    , REML = TRUE
  )
  
  t1 <- sd.B0/sigma(tempmod)
  t2 <- rho.B01*sd.B1/sigma(tempmod)
  t3 <- sqrt((sd.B1/sigma(tempmod))^2 - t2^2)
  
  t4 <- sd.tip/sigma(tempmod)
  t5 <- rho.slopetip*sd.slope/sigma(tempmod)
  t6 <- sqrt((sd.slope/sigma(tempmod))^2 - t5^2)
  new_y <- simulate(tempmod
    , newparams=list(theta=c(ss/sigma(tempmod)
        , t1, t2, t3
        , t4, t5, t6
        , sd.site/sigma(tempmod)
        )
      , beta = c(beta0, beta1)
      )
    )
  dat$new_y <- new_y[[1]]
  t4 <- proc.time()
hackedmod <- glmmTMBhacked(new_y ~ X  
	+ (1 | sp:site)
	+ (1 + X | sp) 
	+ (1 + X | obs)
	+ (1 | site)
  , data=dat
  , phyloZ = phyZ
  , phylonm = c("sp", "sp:site")
  , doFit=TRUE
  , dispformula = ~1
  , control=glmmTMBControl(optCtrl=list(trace=1))
  , REML = TRUE
  ) 
t5 <- proc.time()
glmmTMBtime <- t2-t1 + (t5-t4)
}
print(summary(hackedmod))

glmmTMB_list <- list(hackedmod, glmmTMBtime)


saveRDS(glmmTMB_list, file=paste("datadir/glmmTMB",numsite,size,seed,"rds",sep="."))

#rdnosave()


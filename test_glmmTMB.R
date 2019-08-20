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

mod1 <- glmmTMBhacked(Y ~ X  
                             # + (1|sp) + (0 + X| sp) 
                             + (1 |sp)
                      # + (0 + X |sp)
                             , data=dat
                             , phyloZ = phyZ
                             , phylonm = c("sp", "sp:site")
                             , doFit=TRUE
                             # , dispformula = ~1
                             , control=glmmTMBControl(optCtrl=list(trace=1,iter.max=1e5,eval.max=1e5))
                             , REML = FALSE
                             , lambhack=FALSE
  ) 
mod1

mod2 <- glmmTMBhacked(Y ~ X  
                      # + (1|sp) + (0 + X| sp) 
                      + (1 |obs)
                      # + (0 + X | obs)
                      , data=dat
                      , phyloZ = phyZ
                      , phylonm = c("sp", "sp:site")
                      , doFit=TRUE
                      # , dispformula = ~1
                      , control=glmmTMBControl(optCtrl=list(trace=1,iter.max=1e5,eval.max=1e5))
                      , REML = FALSE
                      , lambhack=TRUE
) 
mod2

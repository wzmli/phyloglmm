### Fitting phyloglmm with TMB

library(dplyr)
library(glmmTMB)
library(lme4)
library(Matrix)
library(phyloglmm)

t1 <- proc.time()

dat <- (dat
  %>% mutate(sp = factor(sp),
             obs = sp
      )
)

print(dat)

quit()

if (numsite == "ss") {
  glmmTMBmod <- phylo_glmmTMB(y_main ~ X
    + (1 + X|sp)
    , data=dat
    , phylo = phy
    , phylonm = "sp"
    , doFit=TRUE
    , control=glmmTMBControl(optCtrl=list(trace=1,iter.max=1e4,eval.max=1e4))
    , REML = FALSE
  )
}

if (numsite %in% c("ms","mms")){
  glmmTMBmod <- phylo_glmmTMB(y_all ~ X
                              + (1 | sp:site)
                              + (1 + X | sp)
                              + (1 + X | obs)
                              + (1 | site)
                            , data=dat
                            , phylo = phy
                            , phylonm = c("sp", "sp:site")
                            , doFit=TRUE
                            # , dispformula = ~1
                            , control=glmmTMBControl(optCtrl=list(trace=1,
                                                                  iter.max=1e4,eval.max=1e4))
                            , REML = FALSE
                              )
}
t2 <- proc.time()
glmmTMBtime <- t2-t1
print(summary(glmmTMBmod))

glmmTMB_list <- list(glmmTMBmod, glmmTMBtime)


print(summary(glmmTMBmod))
print(glmmTMBmod$fit$convergence)

saveRDS(glmmTMB_list, file=paste("datadir/glmmTMB/glmmTMB",numsite,size,tree_seed,"rds",sep="."))

#rdnosave()


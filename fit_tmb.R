### Fitting phyloglmm with TMB

library(dplyr)
library(glmmTMB)
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

TMBstruc <- glmmTMB(Y ~ X  + (1+X|sp)
  , data=dat
  , doFit=FALSE) # doFit=FALSE) in BB's update

TMBstruc_new <- modify_TMBstruc(TMBstruc,phy,phylonm="sp",phyloZ=phyZ)

# glmmTMB_fit <- glmmTMB:::fitTMB(TMBstruc_new)
# tt <- tidy(glmmTMB_fit,scales=c(ran_pars="vcov",fixed=NA))
# glmmTMB_res <- tt[,c("term","estimate","std.error")]

glmmTMB_fit <- fit_TMBstruc(TMBstruc_new)
print(glmmTMB_fit$report)


glmmTMB_res <- matrix(c(glmmTMB_fit$fit$par[1:2]
 	, exp(glmmTMB_fit$fit$par[3:4])^2)
	, ncol=1, 
	dimnames=list(c("intercept","cofactor","var.phylo","var.obs"),NULL)
)

print(glmmTMB_res)

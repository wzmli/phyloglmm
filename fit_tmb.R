### Fitting phyloglmm with TMB

library(dplyr)
library(glmmTMB)
library(Matrix)

phyZ <- phylo.to.Z(phy)

dat <- (dat
        %>% rowwise()
        %>% mutate(phylo=paste("t",sp,sep="")
                   , obs=phylo
        )
)

dat <- data.frame(dat)

TMBstruc <- glmmTMB(Y ~ noise  + (1|phylo) + (1|obs)
  , data=dat
  , debug=TRUE) # doFit=FALSE) in BB's update

TMBstruc_new <- modify_TMBstruc(TMBstruc,phy,phylonm="phylo",phyloZ=phyZ,sp=dat$phylo)

# glmmTMB_fit <- glmmTMB:::fitTMB(TMBstruc_new)
# tt <- tidy(glmmTMB_fit,scales=c(ran_pars="vcov",fixed=NA))
# glmmTMB_res <- tt[,c("term","estimate","std.error")]

glmmTMB_fit <- fit_TMBstruc(TMBstruc_new)
print(summary(glmmTMB_fit))

glmmTMB_res <- matrix(c(glmmTMB_fit$fit$par[1:2]
 	, exp(glmmTMB_fit$fit$par[3:4])^2)
	, ncol=1, 
	dimnames=list(c("intercept","cofactor","var.phylo","var.obs"),NULL)
)

print(glmmTMB_res)

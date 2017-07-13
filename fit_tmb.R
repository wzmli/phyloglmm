### Fitting phyloglmm with TMB

library(glmmTMB)
library(Matrix)
phyZ <- phylo.to.Z(phy)

dat$obs <- dat$sp

TMBstruc <- glmmTMB(Y ~ X  + (1|obs) + (1|sp) + (0 + X|obs) + (0+X|sp)
  , data=dat
  , debug=TRUE) # doFit=FALSE) in BB's update

TMBstruc_new <- modify_TMBstruc(TMBstruc,phy,phylonm="sp",phyloZ=phyZ,sp)

# glmmTMB_fit <- glmmTMB:::fitTMB(TMBstruc_new)
# tt <- tidy(glmmTMB_fit,scales=c(ran_pars="vcov",fixed=NA))
# glmmTMB_res <- tt[,c("term","estimate","std.error")]

glmmTMB_fit <- fit_TMBstruc(TMBstruc_new)
print(glmmTMB_fit)

# glmmTMB_res <- matrix(c(glmmTMB_fit$fit$par[1:2],exp(glmmTMB_fit$fit$par[3:4])^2),
#                       ncol=1, dimnames=list(c("intercept","cofactor","var.phylo","var.obs"),
#                       NULL))
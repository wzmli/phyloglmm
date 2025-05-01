## an attempt to do phyloglmm models by post-processing rather than needing to copy/hack setups:

library(reformulas)
library(glmmTMB)
library(lme4)

## phyloglmm loads, not currently working; need it for Z construction etc.
library(phyloglmm)
phylo <- garamszegi_phy
phyloZ <- phylo.to.Z(phylo)

datG <- garamszegi_simple
## create observation-level grouping variable
datG$obs <- factor(seq(nrow(datG)))
datG$sp <- factor(datG$phylo)

m1 <- glmmTMB(phen ~ cofactor + (obs|sp), doFit = FALSE,
              control = glmmTMB::glmmTMBControl(rank_check = "skip"),
              data = datG)
## now modify
## re-make Ztlist
## modify the appropriate component
mkReTrms(...)

finalizeTMB(TMBStruc, obj, fit, h = NULL, data.tmb.old = NULL)



#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)


t1 <- proc.time()

phyZ <- phylo.to.Z(phy,stand=FALSE)
phyZ <- phyZ[order(rownames(phyZ)),]

dat <- (dat
  %>% mutate(obs = sp)
  %>% ungroup()
)	


# phyonly <- phylo_lmm(y_phy ~ X
#   + (1 + X | sp)
#   , data=dat
#   , phylonm = c("sp","sp:site")
#   , phylo = phy
#   , phyloZ=phyZ
#   , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#   , REML = FALSE
# )
# 
# nophy <- phylo_lmm(y_nophy ~ X
#   + (1 + X | obs)
#   , data=dat
#   , phylonm = c("sp","sp:site")
#   , phylo = phy
#   , phyloZ=phyZ
#   , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#   , REML = FALSE
# )
# 
# main <- phylo_lmm(y_main ~ X
#   + (1 + X | sp)
#   + (1 + X | obs)
#   , data=dat
#   , phylonm = c("sp","sp:site")
#   , phylo = phy
#   , phyloZ=phyZ
#   , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#   , REML = FALSE
# )
# 
# interaction_mod <- phylo_lmm(y_interaction ~ X
#      + (1 | site:obs)
#      , data=dat
#      , phylonm = c("sp","sp:site")
#      , phylo = phy
#      , phyloZ=phyZ
#      , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#      , REML = FALSE
# )
# 
# interaction_phyonly <- phylo_lmm(y_interactionphyonly ~ X
#    + (1 | sp:site)
#    , data=dat
#    , phylonm = c("sp","sp:site")
#    , phylo = phy
#    , phyloZ=phyZ
#    , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#    , REML = FALSE
# )
# 
# interactionphy <- phylo_lmm(y_interactionphy ~ X
#   + (1 | sp:site)
#   + (1 | obs:site)
#   , data=dat
#   , phylonm = c("sp","sp:site")
#   , phylo = phy
#   , phyloZ=phyZ
#   , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#   , REML = FALSE
# )

# fullmod <- phylo_lmm(y ~ X
#   + (1 + X | sp)
#   + (1 + X | obs)
#   + (1 | site)
#   + (1 | sp:site)
#   # + (1 | obs:site)
#   , data=dat
#   , phylonm = c("sp","sp:site")
#   , phylo = phy
#   , phyloZ=phyZ
#   , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#   , REML = FALSE
# )


fullmod <- phylo_lmm(y_main ~ X
                     + (1 + X | sp)
                     + (1 + X | obs)
                     # + (1 | site)
                     # + (1 | sp:site)
                     # + (1 | obs:site)
                     , data=dat
                     , phylonm = c("sp","sp:site")
                     , phylo = phy
                     , phyloZ=phyZ
                     , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                     , REML = FALSE
)

# ll <- list(phyonly, nophy, main
#   , interaction_mod
#   , interaction_phyonly
#   , interactionphy, fullmod
# )

saveRDS(fullmod, file=paste("lme4test/lme4",numsite,size,tree_seed,"rds",sep="."))
#rdnosave()

### rep li and ives example 

library(ape)
library(Matrix)
library(lme4)
library(dplyr)
library(pez)
library(glmmTMB)

source("new_phylo_setup.R")
source("glmmTMBhacked.R")
dd <- readRDS("dune_dat.RDS")

phy <- dd[[1]]
dat <- dd[[2]]
REs <- dd[[3]]

dd <- data.frame(dat)
print(dd %>% count(site))
## 28 sp in each site

phyZ <- phylo.to.Z(phy,stand = TRUE)
phyZ <- phyZ[order(rownames(phyZ)),]

dat <- (dat
        %>% rowwise()
        %>% mutate(obs = sp
        )
)

## If you hack mkBlist in phylo_lmm and switch the order in the kronecker product step (bad!)
lme4time_1 <- system.time(
  lme4fit_1 <- phylo_lmm(Y ~ 1 + log.sla + annual 
# + (1|obs) 
# + (1|sp)
+ (1 | sp:site)
# + (0 + log.sla | site)
# + (1|site) 
, data=dat
, phylonm = c("sp","sp:site")
, phyloZ=phyZ
, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
, REML = FALSE
  )
)

re.site <- REs[[1]]
re.sp <- REs[[2]]
re.sp.phy <- REs[[3]]
re.nested.phy <- REs[[4]]
nsite <- length(unique(dat[["site"]]))
re.sla = list(unname(unlist(dat["log.sla"])), site = dat$site, covar = diag(nsite))


peztime_1 <- system.time(
  pezfit_1 <-  communityPGLMM(formula = "Y ~ 1 + log.sla + annual"
     , data = dat
     , family = "gaussian"
     , sp = dat$sp
     , site = dat$site
     , random.effects = list(re.sp.phy
    , re.sp
    , re.nested.phy
    , re.sla
    , re.site
     )
     , REML = F
     , verbose = F
     # , s2.init = c(1.5, rep(0.01, 4))
  )
)

# 
print(peztime_1)
print(summary(pezfit_1))
print(lme4time_1)
print(summary(lme4fit_1))

print(peztime_1/lme4time_1)


peztime_1 <- system.time(
  pezfit_1 <-  communityPGLMM(formula = "Y ~ 1 + log.sla + annual"
     , data = dat
     , family = "gaussian"
     , sp = dat$sp
     , site = dat$site
     , random.effects = list(#re.sp.phy
    # , re.sp
    re.nested.phy
    # , re.sla
    # , re.site
     )
     , REML = F
     , verbose = F
     # , s2.init = c(1.5, rep(0.01, 4))
  )
)

chol_time <- system.time(cholphyt <- t(chol(phyZ %*% t(phyZ))))


hackedmod <- glmmTMBhacked(Y ~ 1 + log.sla + annual
  + (1|obs:site)
  # + (1|sp)
  # + (0 + log.sla | site)
  # + (1|site)
  # + (1 | sp:site)
  , data=dat
  , phylonm = c("sp","sp:site")
  , phyloZ = cholphyt
  , doFit=TRUE
  # , dispformula = ~1
  , control=glmmTMBControl(optCtrl=list(trace=1,iter.max=1e5,eval.max=1e5))
  , REML = FALSE
  , lambhack=TRUE
)

hackedmod

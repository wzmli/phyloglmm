### rep li and ives example 

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

dd <- data.frame(dat)
print(dd %>% count(site))
## 28 sp in each site

phy <- get_phylo(veg.long = dune.veg2
              , trait = dune.traits2[c(1, 2)]
              , trait.re = c("log.sla")
              , phylo = dune.phylo2
              , trans = "log")

phyZ <- phylo.to.Z(phy)

dat <- (dat
        %>% mutate(obs = sp)
)	

# print(phyZ)

lme4time <- system.time(
  lme4fit <- phylo_lmm(Y ~ 1 + log.sla + annual 
  	+ (1|obs) 
	  + (1|sp)
	  + (0 + site|sp)
  	+ (1 | log.sla)
  	+ (1|site) 
                          , data=dat
                          , phylonm = "sp"
                          , nsp = 28
                          , phylo = phy
                          , phyloZ=phyZ
                          , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  )
)

REs <- get_RE(veg.long = dune.veg2
              , trait = dune.traits2[c(1, 2)]
              , trait.re = c("log.sla")
              , phylo = dune.phylo2
              , trans = "log")

## get_RE returns a list of random effects in this order:
## re.site
## re.sp
## re.phy
## re.nested.phy

re.site <- REs[[1]]
re.sp <- REs[[2]]
re.sp.phy <- REs[[3]]
re.nested.phy <- REs[[4]]
re.sla = list(unname(unlist(dat["log.sla"])), site = dat$site, covar = diag(nsite))




wtf <-  communityPGLMM(formula = "Y ~ 1 + log.sla + annual",
                       data = dat, family = "gaussian", 
                       sp = dat$sp, site = dat$site, 
                       random.effects = list(re.sp
                                             , re.sp.phy
                                             , re.nested.phy
                                             , re.sla
                                             , re.site), 
                       REML = T, verbose = F, 
                       # s2.init = c(1.5, rep(0.01, 3)), 
                       reltol = 10e-20, maxit = 10000)

summary(wtf)

print(summary(lme4fit))

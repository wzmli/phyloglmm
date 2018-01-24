## Compare chol pez and lme4


### rep li and ives example 

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

# debug(phylo_lmm)
# debug(modify_phylo_retrms)

dd <- data.frame(dat)
print(dd %>% count(site))
## 28 sp in each site

phy <- get_phylo(veg.long = dune.veg2
                 , trait = dune.traits2[c(1, 2)]
                 , trait.re = c("log.sla")
                 , phylo = dune.phylo2
                 , trans = "log"
)

phyZ <- phylo.to.Z(phy)
phyZ <- phyZ[order(rownames(phyZ)),]

dat <- (dat
        %>% rowwise()
        %>% mutate(obs = sp
        )
)	

lme4_nonnested <- phylo_lmm(Y ~ 1 + log.sla + annual 
    + (1|obs) 
    + (1|sp)
    # + (1 | site:sp)
    + (0 + log.sla | site)
    + (1|site) 
    , data=dat
    , phylonm = c("sp","site:sp")
    , nsp = 28
    , phylo = phy
    , phyloZ=phyZ
    , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
)

lme4_nested <- phylo_lmm(Y ~ 1 + log.sla + annual 
  + (1 | site:sp)
  , data=dat
  , phylonm = c("sp","site:sp")
  , nsp = 28
  , phylo = phy
  , phyloZ=phyZ
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
)

REs <- get_RE(veg.long = dune.veg2
              , trait = dune.traits2[c(1, 2)]
              , trait.re = c("log.sla")
              , phylo = dune.phylo2
              , trans = "log"
)

# get_RE returns a list of random effects in this order:
# re.site
# re.sp
# re.phy
# re.nested.phy

re.site <- REs[[1]]
re.sp <- REs[[2]]
re.sp.phy <- REs[[3]]
re.nested.phy <- REs[[4]]
re.sla = list(unname(unlist(dat["log.sla"])), site = dat$site, covar = diag(nsite))

re.hacked <- re.sp.phy
re.hacked$covar <- kronecker(diag(20),re.sp.phy$covar)
dimnames(re.hacked$covar)[[1]] <- rep(dimnames(re.sp.phy$covar)[[1]],20)
dimnames(re.hacked$covar)[[2]] <- rep(dimnames(re.sp.phy$covar)[[2]],20)

pez_nonnested <-  communityPGLMM(formula = "Y ~ 1 + log.sla + annual"
  , data = dat
  , family = "gaussian"
  , sp = dat$sp
  , site = dat$site
  , random.effects = list(re.sp
    , re.sp.phy
    , re.sla
    , re.site
    )
  , REML = T
  , verbose = F
  , s2.init = c(1.5, rep(0.01, 4))
  , reltol = 10e-10
  , maxit = 1000
)

pez_nested <-  communityPGLMM(formula = "Y ~ 1 + log.sla + annual"
  , data = dat
	, family = "gaussian"
	, sp = dat$sp
	, site = dat$site
	, random.effects = list(re.nested.phy)
	, REML = T
	, verbose = F
	# , s2.init = c(1.5, rep(0.01, 3))
	, reltol = 10e-10
	, maxit = 1000
)

hacked_nested <-  hacked_pez(formula = "Y ~ 1 + log.sla + annual"
  , data = dat
  , family = "gaussian"
  , sp = dat$sp
  , site = dat$site
  , random.effects = list(#re.sp
  # , re.sp.phy
  re.hacked
  # , re.nested.phy
  # , re.sla
  # , re.site
  )
  , REML = T
  , verbose = F
  # , s2.init = c(1.5, rep(0.01, 4))
  , reltol = 10e-10
  , maxit = 1000
)



print(summary(pez_nonnested))
print(summary(lme4_nonnested))

print(summary(pez_nested))
print(summary(lme4_nested))
print(summary(hacked_nested))


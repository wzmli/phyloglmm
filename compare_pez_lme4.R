## Compare pez and lme4 nested
## Compare pez


### rep li and ives example 

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

maxit <- 10000000
reltol <- 0.00000001

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

lme4_nested <- phylo_lmm(Y ~ 1 + log.sla + annual
  + (1 | sp:site)
  , data=dat
  , phylonm = c("sp","sp:site")
  , nsp = 28
  , phylo = phy
  , phyloZ=phyZ
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  , REML = FALSE
)

REs <- get_RE(veg.long = dune.veg2
  , trait = dune.traits2[c(1, 2)]
  , trait.re = c("log.sla")
  , phylo = dune.phylo2
  , trans = "log"
  , stand = FALSE
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

hacked_nested <-  hacked_pez(formula = "Y ~ 1 + log.sla + annual"
  , data = dat
  , family = "gaussian"
  , sp = dat$sp
  , site = dat$site
  , random.effects = list(re.hacked)
  , REML = F
  , verbose = F
# , s2.init = c(1.5, rep(0.01, 4))
  , reltol = reltol
  , maxit = maxit
)

print(summary(lme4_nested))
print(summary(hacked_nested))


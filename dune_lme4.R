### rep li and ives example 

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

source('new_phylo_setup.R', echo=TRUE)


# maxit <- 100000000
# reltol <- 0.000000001

# debug(phylo_lmm)
# debug(phylo_lmm2)
# debug(lFormula2)
# debug(mkReTrms2)
# debug(mkBlist2)
# debug(modify_phylo_retrms2)
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
  # %>% arrange(site)
)	

lme4time_1 <- system.time(
  lme4fit_1 <- phylo_lmm2(Y ~ 1 + log.sla + annual 
		+ (1|obs) 
		+ (1|sp)
    + (1 | sp:site)
		 + (0 + log.sla | site)
		+ (1|site) 
		, data=dat
		, phylonm = c("sp","sp:site","site:sp")
		, nsp = 28
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
		, REML = FALSE
  )
)

lme4fit_12 <- phylo_lmm2(Y ~ 1 + log.sla + annual 
                        + (1|obs) 
                        + (1|sp)
                        + (1 | site:sp)
                        + (0 + log.sla | site)
                        + (1|site) 
                        , data=dat
                        , phylonm = c("sp","sp:site","site:sp")
                        , nsp = 54
                        , phylo = phy
                        , phyloZ=phyZ
                        , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                        , REML = FALSE
)


lme4time_2 <- system.time(
  lme4fit_2 <- phylo_lmm(Y ~ 1 + log.sla + annual
		+ (1|obs)
		+ (1|sp)
		+ (1 | sp:site)
		+ (0 + log.sla | site)
		+ (1|site)
		, data=dat
		, phylonm = c("sp","sp:site")
		, nsp = 28
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
		, REML = FALSE
	)
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


# peztime_1 <- system.time(
#   pezfit_1 <-  communityPGLMM(formula = "Y ~ 1 + log.sla + annual"
#     , data = dat
#     , family = "gaussian"
#     , sp = dat$sp
#     , site = dat$site
#     , random.effects = list(#re.sp
#     #  , re.sp.phy
#       re.nested.phy
#       #, re.sla
#       #, re.site
#       )
#     , REML = F
#     , verbose = F
#     , s2.init = c(1.5, rep(0.01, 4))
#   )
# )
# 
# peztime_2 <- system.time(
# 	pezfit_2 <-  communityPGLMM(formula = "Y ~ 1 + log.sla + annual"
# 		, data = dat
# 		, family = "gaussian"
# 		, sp = dat$sp
# 		, site = dat$site
# 		, random.effects = list(
# 		    re.sp
# 			, re.sp.phy
# 			, re.nested.phy
# 			, re.site
# 		)
# 	, REML = T
# 	, verbose = F
# 	, s2.init = c(1.5, rep(0.01, 3))
# 	, reltol = 10e-10
# 	, maxit = 1000
# 	)
# )
# 
# 
# print(peztime_1)
# print(summary(pezfit_1))
print(lme4time_1)
print(summary(lme4fit_1))

print(summary(lme4fit_12))

# print(peztime_2)
# print(summary(pezfit_2))
print(lme4time_2)
print(summary(lme4fit_2))


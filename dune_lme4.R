### rep li and ives example 

library(ape)
library(Matrix)
library(lme4)
library(dplyr)


# maxit <- 100000000
# reltol <- 0.000000001

# debug(phylo_lmm)
# debug(lFormula)
# debug(mkReTrms)
# debug(mkBlist)

dd <- data.frame(dat)
print(dd %>% count(site))
## 28 sp in each site

phy <- get_phylo(veg.long = dune.veg2
	, trait = dune.traits2[c(1, 2)]
	, trait.re = c("log.sla")
	, phylo = dune.phylo2
	, trans = "log"
)

phyZ <- phylo.to.Z(phy,stand = TRUE)
phyZ <- phyZ[order(rownames(phyZ)),]

dat <- (dat
#  %>% rowwise()
  %>% mutate(obs = factor(sp)
  )
)

# dat2 <- rbind(dat,dat[560,])
# dat2[561,"Y"] <- 0.5

lme4time_1 <- system.time(
  lme4fit_1 <- phylo_lmm(Y ~ 1 + log.sla + annual 
		+ (1|sp)
		+ (1|obs) 
      + (1 | site:sp)
		   + (0 + log.sla | site)
#		  + (1|site) 
		, data=dat
		, phylonm = c("sp","site:sp")
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
		, REML = FALSE
  )
)

print(summary(lme4fit_1))

REs <- get_RE(veg.long = dune.veg2
	, trait = dune.traits2[c(1, 2)]
	, trait.re = c("log.sla")
	, phylo = dune.phylo2
	, trans = "log"
	, stand = TRUE
)

re.site <- REs[[1]]
re.sp <- REs[[2]]
re.sp.phy <- REs[[3]]
re.nested.phy <- REs[[4]]
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
#      , re.site
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

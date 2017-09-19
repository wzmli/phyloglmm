#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

print(covmat)
print(cov2cor(covmat))

phyZ <- phylo.to.Z(phy)

dat <- (dat
	%>% mutate(obs = sp)
)	

debug(phylo_lmm)
debug(modify_phylo_retrms)


lme4timecor <- system.time(
	lme4_interactionRE <- phylo_lmm(Y ~ site_name 
		# + (1|obs) 
		+ (1|sp) 
		# + (0+site_name|sp) 
		+ (1|sp:site_name)
		, data=dat
		, phylonm = c("sp","sp:site_name")
		, phylo = phy
		, phyloZ=phyZ
		, nsp = nspp
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
	)
)

print(summary(lme4_interactionRE))
# 
# lme4_fullRE <- phylo_lmm(Y ~ site_name 
#                                   # + (1|obs) 
#                                   # + (1|sp) 
#                                   + (0+site_name|sp) 
#                                   # + (1|site_name:sp)
#                                   , data=dat
#                                   , phylonm = c("sp","site_name:sp")
#                                   , phylo = phy
#                                   , phyloZ=phyZ
#                                   , nsp = nspp
#                                   , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
# )
# 
# print(summary(lme4_fullRE))

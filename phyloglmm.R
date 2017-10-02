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

single_site_dat <- (dat 
	%>% filter(site == 1)
)

#debug(phylo_lmm)
#debug(modify_phylo_retrms)


lme4time <- system.time(
	lme4_single_site <- phylo_lmm(Y ~ 1
		# + (1|obs) 
		+ (1|sp) 
		# + (1|sp:site)
		, data=dat
		, phylonm = c("sp","sp:site")
		, phylo = phy
		, phyloZ=phyZ
		, nsp = nspp
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
	)
)

print(summary(lme4_single_site))


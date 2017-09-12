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

lme4timecor <- system.time(
	lme4fitcor <- phylo_lmm(Y ~ site_name 
		+ (1|obs) 
		+ (1|sp) 
		+ (0+site_name|obs) 
		+ (0+site_name|sp)
		, data=dat
		, phylonm = "sp"
		, phylo = phy
		, phyloZ=phyZ
		, nsp = nspp
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
	)
)

print(lme4timecor)
print(summary(lme4fitcor))

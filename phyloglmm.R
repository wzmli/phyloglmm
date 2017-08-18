#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

print(covmat)

phyZ <- phylo.to.Z(phy)

dat <- (dat
	%>% mutate(obs = sp)
)	

lme4time <- system.time(
	lme4fit <- phylo_lmm(Y ~ X
		#+ (1 | obs)
		+ (1 | sp)
		#+ (1+X|sp)
		#+ (1+X|obs)
		# + (0 + X | obs)
		+ (0 + X | sp)
		, data=dat
		, phylonm = "sp" 
		, sp = dat$sp
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
	)
)

print(summary(lme4fit))

lme4timecor <- system.time(
	lme4fitcor <- phylo_lmm(Y ~ X
		#+ (1 | obs)
		#+ (1 | sp)
		+ (1+X|sp)
		#+ (1+X|obs)
		# + (0 + X | obs)
		#+ (0 + X | sp)
		, data=dat
		, phylonm = "sp"
		, sp = dat$sp
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
	)
)


print(lme4timecor)
print(summary(lme4fitcor))

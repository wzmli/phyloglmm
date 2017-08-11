#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

phyZ <- phylo.to.Z(phy)

# debug(phylo_lmm)
# debug(modify_phylo_retrms)
# debug(mkLmerDevfun)


if(single.site==TRUE){
	dat <- (dat
		%>% mutate(obs = sp)
)
	lme4time <- system.time(
	lme4fit <- phylo_lmm(Y ~ X
	                     # + (1 | obs)
	                     # + (1 | sp)
	                     + (1+noise|sp)
	                     # + (1+noise|obs)
	                     # + (0 + noise | obs)
	                     # + (0 + noise | sp)
	   , data=dat
		, phylonm = "sp" 
		, sp = dat$sp
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
		, correlated=TRUE)
	)
}
# 
if(single.site==FALSE){
dat$obs <- dat$sp
lme4time <- system.time(
	lme4fit <- phylo_lmm(Y ~ X + (1|obs) + (1|sp)  + (0 + X|obs) + (0+X|sp)
	, data=dat
	, phylonm = "sp"
	, sp = dat$sp
	, phylo = phy
	, phyloZ=phyZ
	, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
	, correlated = FALSE)
)
}
print(lme4time)
print(summary(lme4fit))

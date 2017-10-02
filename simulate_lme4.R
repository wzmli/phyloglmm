### Simulating via lme4 object

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

# debug(phylo_lmm)
# debug(modify_phylo_retrms)

single_site_phylo <- phylo_lmm(Y ~ 1
	+ (1|sp) 
	, data=dat
	, phylonm = c("sp","sp:site")
	, phylo = phy
	, phyloZ=phyZ
	, nsp = nspp
	, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
)

print(summary(single_site_phylo))



new_response <- simulate(single_site_phylo,newparams=list(theta=c(1),beta=c(0)))
new_y <- new_response + rnorm(0,10)

dat$y <- new_y

fitmod <- update(single_site_phylo)

print(fitmod)

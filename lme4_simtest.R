## lazy simulation test 

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

sd0vec <- sd1vec <- corvec <- residvec <- numeric(simnum)
for(i in 1:simnum){
	seed <- i
	source("simulate_tree.R",echo=FALSE)
	phyZ <- phylo.to.Z(phy)
	dat <- dat %>% mutate(obs = sp)
	lme4fit <- phylo_lmm(Y ~ X + (1+X|sp)
		, data=dat
		, phylonm = "sp"
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
	)
	cov_df <- as.data.frame(lme4::VarCorr(lme4fit))
	sdcor <- cov_df$sdcor
	sd0vec[i] <- sdcor[1]
	sd1vec[i] <- sdcor[2]
	corvec[i] <- sdcor[3]
	residvec[i] <- sdcor[4]
}

df <- data.frame(sd0vec,sd1vec,corvec,residvec)
print(summary(df))

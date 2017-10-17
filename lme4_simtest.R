## lazy simulation test 

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

sd0vec_ss <- residvec_ss <- sd0vec_ms <- sd1vec_ms <- corvec_ms <- residvec_ms <- numeric(simnum)
for(i in 1:simnum){
	seed <- i
	source("simulate_tree.R",echo=FALSE)
	phyZ <- phylo.to.Z(phy)
	lme4fit_single_site <- phylo_lmm(Y ~ X + (1|sp)
		, data = dat
		, phylonm = "sp"
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
	)
	lme4fit_multi_site <- phylo_lmm(Y ~ X + (1 + X |sp)
		, data=dat
		, phylonm = "sp"
		, phylo = phy
		, phyloZ = phyZ
		, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
		)
	cov_ssdf <- as.data.frame(lme4::VarCorr(lme4fit_single_site))
	cov_msdf <- as.data.frame(lme4::VarCorr(lme4fit_multi_site))
	sdcor_ss <- cov_ssdf$sdcor
	sdcor_ms <- cov_msdf$sdcor
	sd0vec_ss[i] <- sdcor_ss[1]
	residvec_ss[i] <- sdcor_ss[2]
	sd0vec_ms[i] <- sdcor_ms[1]
	sd1vec_ms[i] <- sdcor_ms[2]
	corvec_ms[i] <- sdcor_ms[3]
	residvec_ms[i] <- sdcor_ms[4]
}

df <- data.frame(sd0vec_ss,residvec_ss,sd0vec_ms,sd1vec_ms,corvec_ms,residvec_ms)
print(summary(df))

# print(summary(sd0vec_ss))
# print(summary(residvec_ss))

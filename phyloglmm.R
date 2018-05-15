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

print(head(dat))

#debug(phylo_lmm)
#debug(modify_phylo_retrms)


lme4time <- system.time(
	lme4fit <- phylo_lmm(Y ~ X
		+ (1|sp)
		# + (0 + X|sp)
#		+ (1+X|sp)
		, data=dat
		, phylonm = c("sp","site:sp")
		, phylo = phy
		, phyloZ=phyZ
		, nsp = nspp
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
		)
		, REML = FALSE
	)
)

print(lme4time)

print(summary(lme4fit))


lme4_list <- list(lme4fit, lme4time)

# saveRDS(lme4_list, file=paste("datadir/lme4_test",size,seed,"rds",sep="."))

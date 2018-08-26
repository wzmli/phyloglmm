#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

print(covmat)
print(cov2cor(covmat))

t1 <- proc.time()

phyZ <- phylo.to.Z(phy,stand=TRUE)

dat <- (dat
	%>% mutate(obs = sp)
)	


#debug(phylo_lmm)
#debug(modify_phylo_retrms)


lme4fit <- phylo_lmm(Y ~ X + (1+X|sp)
		, data=dat
		, phylonm = c("sp","site:sp")
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
		)
		, REML = FALSE
	)

t2 <- proc.time()

lme4time <- t2-t1
print(lme4time)

print(summary(lme4fit))


lme4_list <- list(lme4fit, lme4time)

saveRDS(lme4_list, file=paste("datadir/lme4",ss,size,seed,"rds",sep="."))

#rdnosave()

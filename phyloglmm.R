#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)


t1 <- proc.time()

phyZ <- phylo.to.Z(phy,stand=FALSE)

dat <- (dat
	%>% mutate(obs = sp)
)	

t2 <- proc.time()
#debug(phylo_lmm)
#debug(modify_phylo_retrms)

if(numsite == "ss"){
	lme4fit <- phylo_lmm(Y ~ X + (1+X|sp)
		, data=dat
		, phylonm = c("sp","site:sp")
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
		)
		, REML = TRUE
	)
t3 <- proc.time()
lme4time <- t3-t1
}

if(numsite == "ms"){
  
t4 <- proc.time()
    lme4fit <- phylo_lmm(new_y ~ X
      + (1 + X | sp)
	 	+ (1 + X | obs)
      + (1 | site)
		+ (1 | sp:site)
      , data=dat
      , phylonm = c("sp","sp:site")
      , phylo = phy
      , phyloZ=phyZ
      , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
      , REML = TRUE
    )
t5 <- proc.time()
lme4time <- t2-t1 + (t5-t4)
}

print(lme4time)

print(summary(lme4fit))


lme4_list <- list(lme4fit, lme4time)

saveRDS(lme4_list, file=paste("datadir/lme4",numsite,size,tree_seed,"rds",sep="."))

#rdnosave()

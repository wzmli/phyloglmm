#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)
library(phyloglmm)

t1 <- proc.time()


dat <- (dat
	%>% mutate(obs = sp)
	%>% ungroup()
	# %>% arrange(sp)
)

if(numsite == "ss"){
  lme4fit <- phylo_lmm(y_main ~ X + (1+X|sp)
                     , data=dat
                     , phylonm = c("sp","site:sp")
                     , phylo = phy
                     , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                     , REML = FALSE
                       )
}

if(numsite == "ms"){
  t4 <- proc.time()
  lme4fit <- phylo_lmm(y_all ~ X
#	 	+ (1 | sp)
#		+ (1 | obs)
                       + (1 + X | sp)
                       + (1 + X | obs)
                       + (1 | site)
                       + (1 | sp:site)
                     , data=dat
                     , phylonm = c("sp","sp:site")
                     , phylo = phy
                     , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                     , REML = FALSE
                       )
}
t2 <- proc.time()

lme4time <- t2 - t1
print(lme4time)

print(summary(lme4fit))


lme4_list <- list(lme4fit, lme4time)

saveRDS(lme4_list, file=paste("datadir/lme4/lme4",numsite,size,tree_seed,"rds",sep="."))



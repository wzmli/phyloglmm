#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)
library(phyloglmm)
library(shellpipes)

loadEnvironments()

t1 <- proc.time()

phyZ <- phylo.to.Z(phy,stand=FALSE)
phyZ2 <- phyZ[order(rownames(phyZ)),]

dat <- (dat
  %>% mutate(sp = factor(sp), obs = sp)
#   %>% mutate(obs = sp)
  %>% ungroup()

	# %>% arrange(sp)
)

# print(dat %>% select(ints,y_main,y_all))

print(numsite)

if(numsite == "ss"){
	lme4fit <- phylo_lmm(y_main ~ X + (1+X|sp)
		, data=dat
		, phylonm = c("sp")
		, phylo = phy
		# , phyloZ = phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
		, REML = FALSE
	)
}


print(lme4fit)


if(numsite %in% c("ms","mms")){
	lme4fit <- phylo_lmm(y_all ~ X
	#	 	+ (1 | sp)
	#		+ (1 | obs)
		+ (1 + X | sp)
		+ (1 + X | obs)
		+ (1 | site)
		+ (1 | sp:site)
		, data=dat
		, phylonm = c("sp", "sp:site")
		, phylo = phy
		, phyloZ = phyZ2 
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
		, REML = FALSE
	)
}

t2 <- proc.time()

lme4time <- t2 - t1
print(lme4time)

print(summary(lme4fit))


lme4_list <- list(lme4fit, lme4time)

rdsSave(lme4_list,target=paste0("datadir/lme4/",targetname()))



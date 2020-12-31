## simulate and fit compound symmetric case
library(ape)
library(Matrix)
library(lme4)
library(dplyr)
#library(pez)
library(phyr)

dat <- (dat
	%>% mutate(obs = sp
		)
	%>% arrange(sp,site)
)

print(dat)
phyZ <- phylo.to.Z(phy)
# 
# debug(phylo_lmm)
# debug(modify_phylo_retrms)

# Vphy <- Vphy[levels(dat$sp),levels(dat$sp)]



phyrfit <- communityPGLMM(y_all ~ X
	+ (1 | sp__)
	+ (X | sp__)
	+ (1 | site)
	+ (1 | sp__@site)
	, data = dat
	, family = "gaussian"
	, tree = phy
	, REML = TRUE
	, cpp = T
	, verbose = T
)

lme4fit <- phylo_lmm(y_all ~ X
	+ (1 | sp)
	+ (0 + X | sp)
	+ (1 | obs)
	+ (0 + X | obs)
	+ (1 | site)
	+ (1 | sp:site)
	, data = dat
	, phylonm = c("sp", "sp:site")
	, phylo = phy
	, phyloZ = phyZ
	, control = lmerControl(check.nobs.vs.nlev="ignore", check.nobs.vs.nRE = "ignore")
	, REML = TRUE
)

print(phyrfit)
print(summary(lme4fit))
res_list <- list(phyrfit, lme4fit)

saveRDS(res_list, file=paste("datadir/compare",numsite,size,tree_seed,"rds",sep="."))

#rdnosave()

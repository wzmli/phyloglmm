## simulate and fit compound symmetric case
library(ape)
library(Matrix)
library(lme4)
library(dplyr)
#library(pez)
library(phyr)

dat <- (dat
	%>% mutate(obs = sp
		, site = site_name
		)
	%>% arrange(sp,site)
)

print(dat)
phyZ <- phylo.to.Z(phy)
# 
# debug(phylo_lmm)
# debug(modify_phylo_retrms)

Vphy <- Vphy[levels(dat$sp),levels(dat$sp)]


tempmod <- phylo_lmm(Y ~ X
	+ (1 | sp:site)
	+ (1 + X | sp)
	+ (1 + X | obs)
	+ (1 | site)
	, data=dat
	, phylonm = c("sp", "sp:site")
	, phylo = phy
	, phyloZ = phyZ
	, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
	, REML = TRUE
)


rho.B01 <- 0
rho.slopetip <- 0

t1 <- sd.B0
t2 <- rho.B01*sd.B1
t3 <- sqrt(sd.B1^2 - t2^2)

t4 <- sd.tip
t5 <- rho.slopetip*sd.slope
t6 <- sqrt(sd.slope^2 - t5^2)

new_y <- simulate(tempmod
	, newparams=list(theta=c(ss
	, t1, t2, t3
	, t4, t5, t6
	, sd.site)/sd.resid
	, beta = c(beta0,beta1)
	, sigma = sd.resid
	)
)

dat$new_y <- new_y[[1]]

#phyrfit <- communityPGLMM(new_y ~ X
#	+ (1 | sp__)
#	+ (X | sp__)
	+ (1 | site)
	+ (1 | sp__@site)
	, data = dat
	, family = "gaussian"
	, tree = phy
	, REML = TRUE
	, cpp = T
	, verbose = T
)

lme4fit <- phylo_lmm(new_y ~ X
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

#saveRDS(res_list, file=paste("datadir/compare",numsite,size,seed,"rds",sep="."))

#rdnosave()

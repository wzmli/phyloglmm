## simulate and fit compound symmetric case
library(ape)
library(Matrix)
library(lme4)
library(dplyr)
library(pez)


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

# random intercept with species independent
sp.int <- list(1, sp = dat$sp, covar = diag(nspp))

# random intercept with species showing phylogenetic covariances
phy.int <- list(1, sp = dat$sp, covar = Vphy)

# random slope with species independent
sp.X <- list(dat$X, sp = dat$sp, covar = diag(nspp))

# random slope with species showing phylogenetic covariances
phy.X <- list(dat$X, sp = dat$sp, covar = Vphy)

# sp:site
phy.interaction <- list(1, sp = dat$sp, covar = Vphy, site = dat$site)

site.int <- list(1, site=dat$site, covar = diag(nsite))

tempmod <- phylo_lmm(Y ~ X
	+ (1 | sp:site)
	+ (1 + X | sp)
	+ (1 + X | obs)
	+ (1 | site)
	, data=dat
	, phylonm = c("sp","sp:site")
	, phylo = phy
	, phyloZ=phyZ
	, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
	, REML = TRUE)


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
, sigma = sd.resid))

dat$new_y <- new_y[[1]]

fittime <- system.time(
	fitpez <- communityPGLMM(new_y ~ X
		, data = dat
		, family = "gaussian"
		, sp = dat$sp
		, site = dat$site
		, random.effects = list(phy.interaction
			 , phy.int, phy.X
			 , sp.int, sp.X
			, site.int
			)
		, REML = TRUE
		, verbose = FALSE
	)
)

print(fittime)
print(summary(fitpez))

saveRDS(list(fittime,fitpez),file=paste("datadir/pez",numsite,size,seed,"rds",sep="."))

#rdnosave()

## simulate and fit compound symmetric case
library(ape)
library(Matrix)
library(lme4)
library(dplyr)
library(pez)


dat <- (dat
	%>% mutate(obs = sp
		, site = factor(site)
		)
	%>% arrange(site,sp)
)

print(dat$site)

print(dat$sp)

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

print(phy.interaction)

site.int <- list(1, site=dat$site, covar = diag(nsite))

print(site.int)


fittime <- system.time(
	fitpez <- communityPGLMM(y_all ~ X
		, data = dat
		, family = "gaussian"
		, sp = dat$sp
		, site = dat$site
		, random.effects = list(phy.interaction
			 , phy.int, phy.X
			 , sp.int, sp.X
			, site.int
			)
		, REML = FALSE
		, verbose = FALSE
	)
)

print(fittime)
print(summary(fitpez))

saveRDS(list(fittime,fitpez),file=paste("datadir/pez/pez",numsite,size,tree_seed,"rds",sep="."))

#rdnosave()

### Fit using pez

library(pez)

# Format input and perform communityPGLMM(
# set up random effects

# random intercept with species independent
re.1 <- list(1, sp = dat$sp, covar = diag(nspp))

# random intercept with species showing phylogenetic covariances
re.2 <- list(1, sp = dat$sp, covar = Vphy)

# random slope with species independent
re.3 <- list(dat$X, sp = dat$sp, covar = diag(nspp))

# random slope with species showing phylogenetic covariances
re.4 <- list(dat$X, sp = dat$sp, covar = Vphy)

# random effect for site
re.site <- list(1, site = dat$site, covar = diag(nsite))


peztime <- system.time(
	pezfit <- communityPGLMM(Y ~ X
		, data = dat
		, family = "gaussian"
		, sp = dat$sp
		, site = dat$site
		, random.effects = list(re.2
			, re.4
			)
		, REML = FALSE
		, verbose = FALSE
		)
)

print(peztime)

print(summary(pezfit))


pez_list <- list(pezfit, peztime)

saveRDS(pez_list, file=paste("datadir/pez",size,seed,"rds",sep="."))


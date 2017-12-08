#### Fitting phyloglmm via lme4 and pez

library(ape)
library(Matrix)
library(lme4)
library(dplyr)
library(pez)


phyZ <- phylo.to.Z(phy)

dat <- (dat
%>% mutate(obs = sp)
)

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


if(platform == "pez"){
	fittime <- system.time(
		fitmod <- communityPGLMM(Y ~ X
			, data = dat
			, family = "gaussian"
			, sp = dat$sp
			, site = dat$site
			, random.effects = list(re.2
				, re.4)
			, REML = FALSE
			, verbose = FALSE
		)
	)
}
	
if(platform == "lme4"){
	fittime <- system.time(
		fitmod <- phylo_lmm(Y ~ X
		  + (1 |sp)
			+ (1 + X|sp)
			, data=dat
			, phylonm = c("sp","sp:site")
			, phylo = phy
			, phyloZ=phyZ
			, nsp = nspp
			, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
		)
	)
}


print(fittime)
print(summary(fitmod))

saveRDS(list(fittime,fitmod),file=paste(datadir,paste(platform,seed,"Rds",sep="."),sep=""))


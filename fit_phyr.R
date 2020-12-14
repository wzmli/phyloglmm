library(ape)
library(phyr)
library(Matrix)
library(dplyr)

fittime <- system.time(
	fitphyr <- communityPGLMM(y_all ~ X
		+ (1 | sp__)
		+ (X | sp__)
		+ (1 | site)
		+ (1 | sp__@site)
		, data = dat
		, family = "gaussian"
		, tree = phy
		, REML = FALSE
		, cpp = TRUE
		)
)

print(fittime)

print(fitphyr)

fit_list <- list(fitphyr, fittime)

saveRDS(fit_list, file=paste("datadir/phyr/phyr",numsite, size, tree_seed,"rds",sep="."))

# rdnosave()


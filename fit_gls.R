### Fitting phyloglmm with gls
library(ape)
library(dplyr)
library(nlme)

dat <- data.frame(dat)
rownames(dat) <- dat$sp

print(phy)

phy_var <- diag(vcv(phy))

varweights <- varFixed(~phy_var)

tt <- system.time(fit_gls <- gls(y_main~X
	, data=dat
	, correlation=corBrownian(phy=phy)
	, weights = varweights
#	, correlation=corPagel(0.5,phy=phy)
	, verbose=FALSE
)
)

print(fit_gls)
print(summary(fit_gls))

gls_list <- list(fit_gls,tt)
saveRDS(gls_list,file=paste("datadir/gls/gls",numsite,size,tree_seed,"rds",sep="."))

#rdnosave()

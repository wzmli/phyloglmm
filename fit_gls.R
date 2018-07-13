### Fitting phyloglmm with gls
library(ape)
library(dplyr)
library(nlme)

dat <- data.frame(dat)
rownames(dat) <- dat$sp

tt <- system.time(fit_gls <- gls(Y~X
	, data=dat
#	, correlation=corBrownian(phy=phy) 
#	, correlation=corPagel(0.5,phy=phy)
	, verbose=FALSE
)
)

print(fit_gls)
print(summary(fit_gls))

gls_list <- list(fit_gls,tt)
saveRDS(gls_list,file=paste("test/TIPonly",size,seed,"rds",sep="."))

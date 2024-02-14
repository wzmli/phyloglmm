### Fitting phyloglmm with gls
library(ape)
library(dplyr)
library(nlme, warn.conflicts = FALSE)
library(shellpipes)

loadEnvironments()


print(dat)

print(phy)

phy_var <- diag(vcv(phy))

varweights <- varFixed(~phy_var)

tt <- system.time(fit_gls <- gls(y_main~X
	, data=dat
	, correlation=corBrownian(phy=phy,
                                  form = ~sp)
	, weights = varweights
#	, correlation=corPagel(0.5,phy=phy)
	, verbose=FALSE
)
)

print(fit_gls)
print(summary(fit_gls))

gls_list <- list(fit_gls,tt)

rdsSave(gls_list,target=paste0("datadir/gls/",targetname()))
# saveRDS(gls_list,file=paste("datadir/gls/"rds",sep="."))

#rdnosave()

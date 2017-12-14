### Fitting phyloglmm with gls
library(ape)
library(dplyr)
library(nlme)

dat <- (dat
	%>% rowwise()
	%>% mutate(phylo=paste("t",sp,sep="")
	, obs=phylo
	)
)

dat <- data.frame(dat)

print(phy)

fit_gls <- gls(Y~1
	, data=dat
	, correlation=corBrownian(phy=phy) 
#	, correlation=corPagel(0.5,phy=phy)
	, verbose=FALSE
)

print(fit_gls)
print(summary(fit_gls))

saveRDS(fit_gls,file=paste("datadir/gls","small",seed,"rds",sep="."))

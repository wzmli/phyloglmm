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

fit_gls <- gls(Y~X
	, data=dat
#	, correlation=corBrownian(phy=phy)
	, verbose=FALSE
)

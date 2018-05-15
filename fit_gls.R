### Fitting phyloglmm with gls
library(ape)
library(dplyr)
library(nlme)

dat <- (dat
	%>% rowwise()
	%>% mutate(sp=paste("t",sp,sep=""))
)

dat <- data.frame(dat)
sp_order <- data.frame(sp=phy$tip.label)
ordered_dat <- left_join(sp_order,dat)
rownames(ordered_dat) <- phy$tip.label

tt <- system.time(fit_gls <- gls(Y~X
	, data=ordered_dat
	, correlation=corBrownian(phy=phy) 
#	, correlation=corPagel(0.5,phy=phy)
	, verbose=FALSE
)
)

print(fit_gls)
print(summary(fit_gls))


gls_list <- list(fit_gls,tt)
saveRDS(gls_list,file=paste("datadir/gls_phyfull",size,seed,"rds",sep="."))

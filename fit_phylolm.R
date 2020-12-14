### Fitting phyloglmm with phylolm
library(ape)
library(dplyr)
library(phylolm)

dat <- data.frame(dat)
#rownames(dat) <- dat$sp

print(head(dat))

tt <- system.time(fit_phylolm <- phylolm(y_main~X
    , data=dat
    , phy=phy 
    , model = "BM"
    , measurement_error = TRUE
)
)

print(fit_phylolm)
print(summary(fit_phylolm))


phylolm_list <- list(fit_phylolm,tt)
saveRDS(phylolm_list,file=paste("datadir/phylolm/phylolm",numsite,size,tree_seed,"rds",sep="."))

#rdnosave()

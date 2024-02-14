### Fitting phyloglmm with phylolm
library(ape)
library(dplyr)
library(phylolm)
library(shellpipes)

loadEnvironments()
dat <- data.frame(dat)
#rownames(dat) <- dat$sp

print(head(dat))

rownames(dat) <- dat$tipname
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

rdsSave(phylolm_list,target=paste0("datadir/phylolm/",targetname()))

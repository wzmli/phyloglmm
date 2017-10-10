## Collect simulated results

library(pez)
library(dplyr)

small_path <- "./datadir/small/"
medium_path <- "./datadir/med/"
large_path <- "./datadir/large/"

small_lme4 <- list.files(path = small_path , pattern = "lme4")
medium_lme4 <- list.files(path = medium_path , pattern = "lme4")
large_lme4 <- list.files(path = large_path , pattern = "lme4")

small_pez <- list.files(path = small_path , pattern = "pez")
medium_pez <- list.files(path = medium_path , pattern = "pez")

get_lme4pars <- function(n,path){
	lme4obj <- readRDS(paste(path, n , sep = ""))
	lme4time <- lme4obj[[1]][[1]]
	covmat <- as.data.frame(lme4::VarCorr(lme4obj[[2]]))
	tempdat <- data.frame(time = lme4time
		, sp = covmat$sdcor[1]
		, X = covmat$sdcor[2]
		, cor = covmat$sdcor[3]
		, Residual = covmat$sdcor[4]
		, dataSize = path
		)
	rownames(tempdat) <- NULL
	return(tempdat)
}


get_pezpars <- function(n,path){
	pezobj <- readRDS(paste(path, n , sep = ""))
	peztime <- pezobj[[1]][[1]]
	pezfit <- pezobj[[2]]
	tempdat <- data.frame(time = peztime
		, sp = sqrt(pezfit$s2r[1])
		, Xsp = sqrt(pezfit$s2r[2])
		, Residual = sqrt(pezfit$s2resid)
		, dataSize = path
)
rownames(tempdat) <- NULL
return(tempdat)
}

pez_res_small <- lapply(small_pez,get_pezpars,path=small_path)
pez_res_medium <- lapply(medium_pez,get_pezpars,path=medium_path)

lme4_res_small <- lapply(small_lme4,get_lme4pars,path=small_path)
lme4_res_medium <- lapply(medium_lme4,get_lme4pars,path=medium_path)
lme4_res_large <- lapply(large_lme4,get_lme4pars,path=large_path)

ps <- rbind_list(pez_res_small)
pm <- rbind_list(pez_res_medium)

lsmall <- rbind_list(lme4_res_small)
lmed <- rbind_list(lme4_res_medium)
llarge <- rbind_list(lme4_res_large)

print(summary(ps))
print(summary(pm))
print(summary(lsmall))
print(summary(lmed))
print(summary(llarge))


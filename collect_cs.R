## Collect compound symmetric cov fits

library(pez)
library(dplyr)

small_path <- "./datadir/cs_small/"
medium_path <- "./datadir/cs_med/"
large_path <- "./datadir/cs_large/"

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
		, sp = covmat$sdcor[2]
		, sp_site = covmat$sdcor[1]
		, Residual = covmat$sdcor[3]
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
		, sp_site = sqrt(pezfit$s2n[1])
		, Residual = sqrt(pezfit$s2resid)
		, dataSize = path
		)
	rownames(tempdat) <- NULL
	return(tempdat)
}

pez_res_small <- lapply(small_pez,get_pezpars,path=small_path)
pez_res_med <- lapply(medium_pez,get_pezpars,path=medium_path)

lme4_res_small <- lapply(small_lme4,get_lme4pars,path=small_path)
lme4_res_med <- lapply(medium_lme4,get_lme4pars,path=medium_path)
lme4_res_large <- lapply(large_lme4,get_lme4pars,path=large_path)

ps <- bind_rows(pez_res_small)
pm <- bind_rows(pez_res_med)

lsmall <- bind_rows(lme4_res_small)
lmed <- bind_rows(lme4_res_med)
llarge <- bind_rows(lme4_res_large)


print(summary(ps))
print(summary(pm))
print(summary(lsmall))
print(summary(lmed))
print(summary(llarge))


df_list <- list(ps,pm,lsmall,lmed,llarge)
saveRDS(df_list,"./results/mscs_simtest.RDS")


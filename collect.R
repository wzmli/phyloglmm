#### Collect simulated results ----
## load packages
library(pez)
library(lme4)
library(dplyr)
library(ggplot2)
library(brms)
library(tidyr)

#### Collect gls results ----

gls_path <- "./datadir/gls/"
gls_res <- list.files(path = gls_path)
gls_results <- function(tt){
  gls_df <- data.frame(resid = numeric(900)
    , phylo_int = numeric(900)
    , phylo_X = NA
    , phylo_cor = NA
    , model = numeric(900)
    , time = numeric(900)
  )
  for(i in 1:length(tt)){
    gls_obj <- readRDS(paste(gls_path,tt[i],sep=""))
    
    gls_df[i,"phylo_int"] <- attr(resid(gls_obj[[1]]),"std")[[1]]
    gls_df[i,"model"] <- tt[i]
    gls_df[i,"time"] <- gls_obj[[2]][[1]]
    }
  return(gls_df)
}

gls_data <- bind_rows(gls_results(gls_res))

#### Collect lme4 single site results ----

lme4_path <- "./datadir/lme4/"
lme4ss_res <- list.files(path = lme4_path, pattern = "ss")
lme4ss_results <- function(tt){
  lme4ss_df <- data.frame(resid = numeric(900)
    , phylo_X = numeric(900)
    , phylo_int = numeric(900)
    , phylo_cor = numeric(900)
    , model = numeric(900)
    , time = numeric(900)
  )
  for(i in 1:length(tt)){
    lme4_obj <- readRDS(paste(lme4_path,tt[i],sep=""))
    covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
    lme4ss_df[i,"resid"] <- (covmat 
      %>% filter(grp=="Residual") 
      %>% select(sdcor) 
      %>% as.numeric()
    )
    lme4ss_df[i, "phylo_X"] <- (covmat 
      %>% filter((grp=="sp") & (var1 =="X")) 
      %>% select(sdcor) 
      %>% as.numeric()
    )
    lme4ss_df[i, "phylo_int"] <- (covmat 
      %>% filter((grp=="sp") 
        & (var1 == "(Intercept)")
        & (is.na(var2))
        ) 
      %>% select(sdcor) 
      %>% as.numeric()
    )
    lme4ss_df[i,"phylo_cor"] <- (covmat 
      %>% filter((grp=="sp") 
        & (var1 == "(Intercept)")
        & (var2 == "X")
        ) 
      %>% select(sdcor) 
      %>% as.numeric()
    )
    lme4ss_df[i,"model"] <- tt[i]
    lme4ss_df[i,"time"] <- lme4_obj[[2]][[1]]
    }
  return(lme4ss_df)
}

lme4ss_data <- lme4ss_results(lme4ss_res)

### collect brms ----

brms_path <- "./datadir/brms/"
brmsss_res <- list.files(path = brms_path, pattern = "ss")
brmsss_results <- function(tt){
  brms_df <- data.frame(resid = numeric(200)
    , phylo_X = numeric(200)
    , phylo_int = numeric(200)
    , phylo_cor = numeric(200)
    , model = numeric(200)
    , time = numeric(200)
  )
  for(i in 1:length(tt)){
    brms_obj <- readRDS(paste(brms_path,tt[i],sep=""))
    sd_dat <- as.data.frame(posterior_samples(brms_obj[[1]],c("^sigma","^sd_","^cor_")))
    brms_df[i,"resid"] <- median(sd_dat[,"sigma"])
    brms_df[i, "phylo_X"] <- median(sd_dat[,"sd_sp__X"])
    brms_df[i, "phylo_int"] <- median(sd_dat[,"sd_sp__Intercept"])
    brms_df[i,"phylo_cor"] <- median(sd_dat[,"cor_sp__Intercept__X"])
    brms_df[i,"model"] <- tt[i]
    brms_df[i,"time"] <- brms_obj[[2]][[1]]
  }
  return(brms_df)
}

brmsss_data <- brmsss_results(brmsss_res)

#### Collect phylolm single site results ----
## need to check out NA cases

phylolm_path <- "./datadir/phylolm/"
phylolm_res <- list.files(path = phylolm_path)
phylolm_results <- function(tt){
  phylolm_df <- data.frame(resid = numeric(900)
    , phylo_int = numeric(900)
    , phylo_X = NA
    , phylo_cor = NA
    , model = numeric(900)
    , time = numeric(900)
  )
  for(i in 1:length(tt)){
    phylolm_obj <- readRDS(paste(phylolm_path,tt[i],sep=""))
    phylolm_df[i,"resid"] <- sqrt(as.numeric(phylolm_obj[[1]]["sigma2_error"]))
    phylolm_df[i,"phylo_int"] <- sqrt(as.numeric(phylolm_obj[[1]]["sigma2"]))
    phylolm_df[i,"model"] <- tt[i]
    phylolm_df[i,"time"] <- phylolm_obj[[2]][[1]]
    }
  return(phylolm_df)
}

phylolm_data <- phylolm_results(phylolm_res)

ssdat <- rbind(gls_data, phylolm_data, lme4ss_data, brmsss_data)


### Collect multiple sites ----

lme4ms_res <- list.files(path = lme4_path, pattern = "ms")
lme4ms_results <- function(tt){
	lme4ms_df <- data.frame(resid = numeric(300)
		, phylo_X = numeric(300)
		, phylo_int = numeric(300)
		, phylo_cor = numeric(300)
		, phylo_interaction = numeric(300)
		, species_X = numeric(300)
		, species_int = numeric(300)
		, species_cor = numeric(300)
		, model = numeric(300)
		, time = numeric(300)
		)
	for(i in 1:length(tt)){
	  lme4_obj <- readRDS(paste(lme4_path,tt[i],sep=""))
	  covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
	  lme4ms_df[i,"resid"] <- (covmat 
      %>% filter(grp=="Residual") 
	    %>% select(sdcor) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i, "phylo_X"] <- (covmat 
	    %>% filter((grp=="sp") & (var1 =="X")) 
	    %>% select(sdcor) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i, "phylo_int"] <- (covmat 
	    %>% filter((grp=="sp") 
	      & (var1 == "(Intercept)")
	      & (is.na(var2))
	      ) 
	    %>% select(sdcor) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i,"phylo_cor"] <- (covmat 
	    %>% filter((grp=="sp") 
	      & (var1 == "(Intercept)")
	      & (var2 == "X")
	      ) 
	    %>% select(sdcor) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i,"phylo_interaction"] <- (covmat
	    %>% filter((grp=="sp:site"))
	    %>% select(sdcor)
	    %>% as.numeric()
	  )
	  lme4ms_df[i, "species_X"] <- (covmat 
	    %>% filter((grp=="obs") & (var1 =="X"))
	    %>% select(sdcor)
	    %>% as.numeric()
	  )
	  lme4ms_df[i, "species_int"] <- (covmat 
	    %>% filter((grp=="obs") 
	      & (var1 == "(Intercept)")
	      & (is.na(var2))
	      ) 
	    %>% select(sdcor) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i,"species_cor"] <- (covmat 
	    %>% filter((grp=="sp") 
	      & (var1 == "(Intercept)")
	      & (var2 == "X")
	      ) 
	    %>% select(sdcor) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i,"model"] <- tt[i]
	  lme4ms_df[i,"time"] <- lme4_obj[[2]][[1]]
	}
	return(lme4ms_df)
}

lme4ms_data <- lme4ms_results(lme4ms_res)

msdat <- lme4ms_data
#### Save results ----

# data_list <- list(gls_data, lme4ss_data, lme4ss_slope_data, lme4ms_data, pez_data, lme4ms_slope_data , pez_slope_data, lme4cs_data, pez_cs_data, lme4cs_slope_data, pez_cs_slope_data)
# saveRDS(data_list,file="./datadir/result_list.RDS")
# rdsave(ssdat, msdat)

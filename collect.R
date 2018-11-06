#### Collect simulated results ----
## load packages
library(pez)
library(lme4)
library(dplyr)
library(ggplot2)
library(brms)
library(tidyr)
library(phylolm)
library(phyr)

#### Collect gls results ----

gls_path <- "./datadir/gls/"
gls_res <- list.files(path = gls_path)
gls_results <- function(tt){
  gls_df <- data.frame(resid = numeric(300)
    , phylo_int = numeric(300)
    , phylo_X = NA
    , phylo_cor = NA
    , B0 = numeric(300)
    , B1 = numeric(300)
    , model = numeric(300)
    , time = numeric(300)
    , convergence = NA
  )
  for(i in 1:length(tt)){
    gls_obj <- readRDS(paste(gls_path,tt[i],sep=""))
    gls_df[i,"phylo_int"] <- as.numeric(gls_obj[[1]]["sigma"])
    B0 <- coef(summary(gls_obj[[1]]))["(Intercept)","Value"]
    B0se <- coef(summary(gls_obj[[1]]))["(Intercept)","Std.Error"]
    B1 <- coef(summary(gls_obj[[1]]))["X","Value"]
    B1se <- coef(summary(gls_obj[[1]]))["X","Std.Error"]
    gls_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
    gls_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
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
  lme4ss_df <- data.frame(resid = numeric(300)
    , phylo_X = numeric(300)
    , phylo_int = numeric(300)
    , phylo_cor = numeric(300)
    , B0 = numeric(300)
    , B1 = numeric(300)
    , model = numeric(300)
    , time = numeric(300)
    , convergence = NA
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
    B0 <- coef(summary(lme4_obj[[1]]))["(Intercept)","Estimate"]
    B0se <- coef(summary(lme4_obj[[1]]))["(Intercept)","Std. Error"]
    B1 <- coef(summary(lme4_obj[[1]]))["X","Estimate"]
    B1se <- coef(summary(lme4_obj[[1]]))["X","Std. Error"]
    lme4ss_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
    lme4ss_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
    lme4ss_df[i,"model"] <- tt[i]
    lme4ss_df[i,"time"] <- lme4_obj[[2]][[1]]
    lme4ss_df[i,"convergence"] <- lme4_obj[[1]]@optinfo$conv$opt
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
    , B0 = numeric(200)
    , B1 = numeric(200)
    , model = numeric(200)
    , time = numeric(200)
    , convergence = numeric(200)
  )
  for(i in 1:length(tt)){
    brms_obj <- readRDS(paste(brms_path,tt[i],sep=""))
    sd_dat <- as.data.frame(posterior_samples(brms_obj[[1]],c("^sigma","^sd_","^cor_")))
    brms_df[i,"resid"] <- median(sd_dat[,"sigma"])
    brms_df[i, "phylo_X"] <- median(sd_dat[,"sd_sp__X"])
    brms_df[i, "phylo_int"] <- median(sd_dat[,"sd_sp__Intercept"])
    brms_df[i,"phylo_cor"] <- median(sd_dat[,"cor_sp__Intercept__X"])
    b_dat <- as.data.frame(posterior_samples(brms_obj[[1]], c("^b")))
    brms_df[i, "B0"] <- as.numeric(between(0
            , quantile(b_dat[,"b_Intercept"], 0.025)
            , quantile(b_dat[,"b_Intercept"], 0.975)
          ))
    brms_df[i,"B1"] <- as.numeric(between(0
            , quantile(b_dat[,"b_X"], 0.025)
            , quantile(b_dat[,"b_X"], 0.975)
          ))
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
  phylolm_df <- data.frame(resid = numeric(300)
    , phylo_int = numeric(300)
    , phylo_X = NA
    , phylo_cor = NA
    , B0 = numeric(300)
    , B1 = numeric(300)
    , model = numeric(300)
    , time = numeric(300)
    , convergence = NA
  )
  for(i in 1:length(tt)){
    phylolm_obj <- readRDS(paste(phylolm_path,tt[i],sep=""))
    phylolm_df[i,"resid"] <- sqrt(as.numeric(phylolm_obj[[1]]["sigma2_error"]))
    phylolm_df[i,"phylo_int"] <- sqrt(as.numeric(phylolm_obj[[1]]["sigma2"]))
    B0 <- coef(summary(phylolm_obj[[1]]))["(Intercept)","Estimate"]
    B0se <- coef(summary(phylolm_obj[[1]]))["(Intercept)","StdErr"]
    B1 <- coef(summary(phylolm_obj[[1]]))["X","Estimate"]
    B1se <- coef(summary(phylolm_obj[[1]]))["X","StdErr"]
    phylolm_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
    phylolm_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
    phylolm_df[i,"model"] <- tt[i]
    phylolm_df[i,"time"] <- phylolm_obj[[2]][[1]]
    }
  return(phylolm_df)
}

phylolm_data <- phylolm_results(phylolm_res)

ssdat <- rbind(gls_data, phylolm_data, lme4ss_data)#, brmsss_data)


### Collect multiple sites ----

lme4_path <- "./datadir/lme4/"
lme4ms_res <- list.files(path = lme4_path, pattern = "ms")
lme4ms_results <- function(tt){
	lme4ms_df <- data.frame(resid = numeric(468)
		, phylo_X = numeric(468)
		, phylo_int = numeric(468)
		, phylo_cor = numeric(468)
		, phylo_interaction = numeric(468)
		, species_X = numeric(468)
		, species_int = numeric(468)
		, species_cor = numeric(468)
		, site_int = numeric(468)
		, B0 = numeric(468)
		, B1 = numeric(468)
		, model = numeric(468)
		, time = numeric(468)
		, convergence = NA
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
	    %>% filter((grp=="obs") 
	      & (var1 == "(Intercept)")
	      & (var2 == "X")
	      ) 
	    %>% select(sdcor) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i,"site_int"] <- (covmat 
	    %>% filter(grp=="site")
	    %>% select(sdcor) 
	    %>% as.numeric()
	  )
	  B0 <- coef(summary(lme4_obj[[1]]))["(Intercept)","Estimate"]
	  B0se <- coef(summary(lme4_obj[[1]]))["(Intercept)","Std. Error"]
	  B1 <- coef(summary(lme4_obj[[1]]))["X","Estimate"]
	  B1se <- coef(summary(lme4_obj[[1]]))["X","Std. Error"]
	  lme4ms_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
	  lme4ms_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
	  lme4ms_df[i,"model"] <- tt[i]
	  lme4ms_df[i,"time"] <- lme4_obj[[2]][[1]]
	  lme4ms_df[i,"convergence"] <- lme4_obj[[1]]@optinfo$conv$opt
	}
	return(lme4ms_df)
}

lme4ms_data <- lme4ms_results(lme4ms_res)

pez_path <- "./datadir/pez/"
pez_res <- list.files(path = pez_path)
pez_results <- function(tt){
  pez_df <- data.frame(resid = numeric(207)
    , phylo_X = numeric(207)
    , phylo_int = numeric(207)
    , phylo_cor = NA
    , phylo_interaction = numeric(207)
    , species_X = numeric(207)
    , species_int = numeric(207)
    , species_cor = NA
    , site_int = numeric(207)
    , B0 = numeric(207)
    , B1 = numeric(207)
    , model = numeric(207)
    , time = numeric(207)
    , convergence = NA
  )
  for(i in 1:length(tt)){
    pez_obj <- readRDS(paste(pez_path,tt[i],sep=""))
    pez_df[i,"resid"] <- sqrt(unlist(pez_obj[[2]]["s2resid"]))
	  pez_df[i,"phylo_interaction"] <- sqrt(unlist(pez_obj[[2]]["s2n"]))
    pez_df[i,"phylo_int"] <- sqrt(unlist(pez_obj[[2]]["s2r"]))[1]
    pez_df[i,"phylo_X"] <- sqrt(unlist(pez_obj[[2]]["s2r"]))[2]
    pez_df[i,"species_int"] <- sqrt(unlist(pez_obj[[2]]["s2r"]))[3]
    pez_df[i,"species_X"] <- sqrt(unlist(pez_obj[[2]]["s2r"]))[4]
    pez_df[i,"site_int"] <- sqrt(unlist(pez_obj[[2]]["s2r"]))[5]
    B0 <- unlist(pez_obj[[2]]["B"])[1]
    B0se <- unlist(pez_obj[[2]]["B.se"])[1]
    B1 <- unlist(pez_obj[[2]]["B"])[2]
    B1se <- unlist(pez_obj[[2]]["B.se"])[2]
    pez_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
    pez_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
    pez_df[i,"model"] <- tt[i]
    pez_df[i,"time"] <- pez_obj[[1]][[1]]
    pez_df[i,"convergence"] <- pez_obj[[2]]["convcode"]
  }
  return(pez_df)
}

pez_data <- pez_results(pez_res)



phyr_path <- "./datadir/phyr/"
phyr_res <- list.files(path = phyr_path)
phyr_results <- function(tt){
  phyr_df <- data.frame(resid = numeric(210)
                       , phylo_X = numeric(210)
                       , phylo_int = numeric(210)
                       , phylo_cor = NA
                       , phylo_interaction = numeric(210)
                       , species_X = numeric(210)
                       , species_int = numeric(210)
                       , species_cor = NA
                       , site_int = numeric(210)
                       , B0 = numeric(210)
                       , B1 = numeric(210)
                       , model = numeric(210)
                       , time = numeric(210)
                       , convergence = NA
  )
  for(i in 1:length(tt)){
    phyr_obj <- readRDS(paste(phyr_path,tt[i],sep=""))
    phyr_df[i,"resid"] <- sqrt(unlist(phyr_obj[[1]]["s2resid"]))
    phyr_df[i,"phylo_interaction"] <- sqrt(unlist(phyr_obj[[1]]["s2n"]))
    phyr_df[i,"phylo_int"] <- sqrt(unlist(phyr_obj[[1]]["s2r"]))[2]
    phyr_df[i,"phylo_X"] <- sqrt(unlist(phyr_obj[[1]]["s2r"]))[4]
    phyr_df[i,"species_int"] <- sqrt(unlist(phyr_obj[[1]]["s2r"]))[1]
    phyr_df[i,"species_X"] <- sqrt(unlist(phyr_obj[[1]]["s2r"]))[3]
    phyr_df[i,"site_int"] <- sqrt(unlist(phyr_obj[[1]]["s2r"]))[5]
    B0 <- unlist(phyr_obj[[1]]["B"])[1]
    B0se <- unlist(phyr_obj[[1]]["B.se"])[1]
    B1 <- unlist(phyr_obj[[1]]["B"])[2]
    B1se <- unlist(phyr_obj[[1]]["B.se"])[2]
    phyr_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
    phyr_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
    phyr_df[i,"model"] <- tt[i]
    phyr_df[i,"time"] <- phyr_obj[[2]][[1]]
    phyr_df[i,"convergence"] <- phyr_obj[[1]]["convcode"]
  }
  return(phyr_df)
}

phyr_data <- phyr_results(phyr_res)




# 
# lme4pez_path <- "./datadir/lme4pez/"
# lme4pez_res <- list.files(path = lme4pez_path, pattern = "ms")
# lme4pez_results <- function(tt){
#   lme4ms_df <- data.frame(resid = numeric(40)
#                           , phylo_X = numeric(40)
#                           , phylo_int = numeric(40)
#                           , phylo_cor = NA
#                           , phylo_interaction = numeric(40)
#                           , species_X = numeric(40)
#                           , species_int = numeric(40)
#                           , species_cor = NA
#                           , site_int = numeric(40)
#                           , B0 = numeric(40)
#                           , B1 = numeric(40)
#                           , model = numeric(40)
#                           , time = numeric(40)
#   )
#   for(i in 1:length(tt)){
#     lme4_obj <- readRDS(paste(lme4pez_path,tt[i],sep=""))
#     covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
#     lme4ms_df[i,"resid"] <- (covmat 
#                              %>% filter(grp=="Residual") 
#                              %>% select(sdcor) 
#                              %>% as.numeric()
#     )
#     lme4ms_df[i, "phylo_X"] <- (covmat 
#                                 %>% filter((grp=="sp") & (var1 =="X")) 
#                                 %>% select(sdcor) 
#                                 %>% as.numeric()
#     )
#     lme4ms_df[i, "phylo_int"] <- (covmat 
#                                   %>% filter((grp=="sp.1") 
#                                              & (var1 == "(Intercept)")
#                                              & (is.na(var2))
#                                   ) 
#                                   %>% select(sdcor) 
#                                   %>% as.numeric()
#     )
#     lme4ms_df[i,"phylo_interaction"] <- (covmat
#                                          %>% filter((grp=="sp:site"))
#                                          %>% select(sdcor)
#                                          %>% as.numeric()
#     )
#     lme4ms_df[i, "species_X"] <- (covmat 
#                                   %>% filter((grp=="obs") & (var1 =="X"))
#                                   %>% select(sdcor)
#                                   %>% as.numeric()
#     )
#     lme4ms_df[i, "species_int"] <- (covmat 
#                                     %>% filter((grp=="obs.1") 
#                                                & (var1 == "(Intercept)")
#                                                & (is.na(var2))
#                                     ) 
#                                     %>% select(sdcor) 
#                                     %>% as.numeric()
#     )
#     lme4ms_df[i,"site_int"] <- (covmat 
#                                 %>% filter(grp=="site")
#                                 %>% select(sdcor) 
#                                 %>% as.numeric()
#     )
#     B0 <- coef(summary(lme4_obj[[1]]))["(Intercept)","Estimate"]
#     B0se <- coef(summary(lme4_obj[[1]]))["(Intercept)","Std. Error"]
#     B1 <- coef(summary(lme4_obj[[1]]))["X","Estimate"]
#     B1se <- coef(summary(lme4_obj[[1]]))["X","Std. Error"]
#     lme4ms_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
#     lme4ms_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
#     lme4ms_df[i,"model"] <- tt[i]
#     lme4ms_df[i,"time"] <- lme4_obj[[2]][[1]]
#   }
#   return(lme4ms_df)
# }
# 
# lme4pez_data <- lme4pez_results(lme4pez_res)



msdat <- rbind(lme4ms_data,pez_data, phyr_data)
#### Save results ----
data_list <- list(ssdat,msdat)

# data_list <- list(gls_data, lme4ss_data, lme4ss_slope_data, lme4ms_data, pez_data, lme4ms_slope_data , pez_slope_data, lme4cs_data, pez_cs_data, lme4cs_slope_data, pez_cs_slope_data)
saveRDS(data_list,file="./datadir/collect.RDS")
# rdsave(ssdat, msdat)

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
library(glmmTMB)
library(MCMCglmm)

#### Collect gls results ----

gls_path <- "./datadir/gls/"
gls_res <- list.files(path = gls_path)
gls_results <- function(tt){
  gls_df <- data.frame(resid = numeric(400)
    , phylo_int = numeric(400)
    , phylo_X = NA
    , phylo_cor = NA
    , B0 = numeric(400)
    , B1 = numeric(400)
    , model = numeric(400)
    , time = numeric(400)
    , convergence = NA
  )
  for(i in 1:length(tt)){
    gls_obj <- readRDS(paste(gls_path,tt[i],sep=""))
    gls_df[i,"phylo_int"] <- as.numeric(gls_obj[[1]]["sigma"])^2
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
  lme4ss_df <- data.frame(resid = numeric(400)
    , phylo_X = numeric(400)
    , phylo_int = numeric(400)
    , phylo_cor = numeric(400)
    , B0 = numeric(400)
    , B1 = numeric(400)
    , model = numeric(400)
    , time = numeric(400)
    , convergence = NA
  )
  for(i in 1:length(tt)){
    lme4_obj <- readRDS(paste(lme4_path,tt[i],sep=""))
    covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
    lme4ss_df[i,"resid"] <- (covmat 
      %>% filter(grp=="Residual") 
      %>% dplyr:::select(vcov) 
      %>% as.numeric()
    )
    lme4ss_df[i, "phylo_X"] <- (covmat 
      %>% filter((grp=="sp") & (var1 =="X")) 
      %>% dplyr:::select(vcov) 
      %>% as.numeric()
    )
    lme4ss_df[i, "phylo_int"] <- (covmat 
      %>% filter((grp=="sp") 
        & (var1 == "(Intercept)")
        & (is.na(var2))
        ) 
      %>% dplyr:::select(vcov) 
      %>% as.numeric()
    )
    lme4ss_df[i,"phylo_cor"] <- (covmat 
      %>% filter((grp=="sp") 
        & (var1 == "(Intercept)")
        & (var2 == "X")
        ) 
      %>% dplyr:::select(vcov) 
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


## single site glmmTMB

glmmTMB_path <- "./datadir/glmmTMB/"
glmmTMBss_res <- list.files(path = glmmTMB_path, pattern = "ss")
glmmTMBss_results <- function(tt){
  glmmTMBss_df <- data.frame(resid = numeric(length(tt))
    , phylo_X = NA
    , phylo_int = NA
    , phylo_cor = NA
    , B0 = NA
    , B1 = NA
    , model = NA
    , time = NA
    , convergence = NA
  )
  for(i in 1:length(tt)){
    glmmTMB_obj <- readRDS(paste(glmmTMB_path,tt[i],sep=""))
    covobj <- VarCorr(glmmTMB_obj[[1]])[["cond"]]
    glmmTMBss_df[i,"resid"] <- attr(covobj,"sc")^2
    glmmTMBss_df[i, "phylo_X"] <- covobj[["sp"]]["X","X"]
    glmmTMBss_df[i, "phylo_int"] <- covobj[["sp"]]["(Intercept)","(Intercept)"]
    glmmTMBss_df[i, "phylo_cor"] <- covobj[["sp"]]["(Intercept)","X"]
    
    
    B0 <- coef(summary(glmmTMB_obj[[1]]))[["cond"]]["(Intercept)","Estimate"]
    B0se <- coef(summary(glmmTMB_obj[[1]]))[["cond"]]["(Intercept)","Std. Error"]
    B1 <- coef(summary(glmmTMB_obj[[1]]))[["cond"]]["X","Estimate"]
    B1se <- coef(summary(glmmTMB_obj[[1]]))[["cond"]]["X","Std. Error"]
    glmmTMBss_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
    glmmTMBss_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
    glmmTMBss_df[i,"model"] <- tt[i]
    glmmTMBss_df[i,"time"] <- glmmTMB_obj[[2]][[1]]
    glmmTMBss_df[i,"convergence"] <- glmmTMB_obj[[1]]$fit$convergence
  }
  return(glmmTMBss_df)
}

glmmTMBss_data <- glmmTMBss_results(glmmTMBss_res)
# glmmTMBss_data <- glmmTMBss_data %>% filter(convergence == 0)
### collect brms ----

brms_path <- "./datadir/brms/"
brmsss_res <- list.files(path = brms_path, pattern = "ss")
brmsss_results <- function(tt){
  brms_df <- data.frame(resid = numeric(300)
    , phylo_X = numeric(300)
    , phylo_int = numeric(300)
    , phylo_cor = numeric(300)
    , B0 = numeric(300)
    , B1 = numeric(300)
    , model = numeric(300)
    , time = numeric(300)
    , convergence = numeric(300)
  )
  for(i in 1:length(tt)){
    brms_obj <- readRDS(paste(brms_path,tt[i],sep=""))
    sd_dat <- as.data.frame(posterior_samples(brms_obj[[1]],c("^sigma","^sd_","^cor_")))
    brms_df[i,"resid"] <- median(sd_dat[,"sigma"])^2
    brms_df[i, "phylo_X"] <- median(sd_dat[,"sd_sp__X"])^2
    brms_df[i, "phylo_int"] <- median(sd_dat[,"sd_sp__Intercept"])^2
    brms_df[i,"phylo_cor"] <- median(sd_dat[,"cor_sp__Intercept__X"] * sd_dat[,"sd_sp__Intercept"] * sd_dat[,"sd_sp__X"])
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

## Collect MCMCglmm

MCMCglmm_path <- "./datadir/MCMCglmm/"
MCMCglmmss_res <- list.files(path = MCMCglmm_path, pattern = "ss")
MCMCglmmss_results <- function(tt){
  MCMCglmm_df <- data.frame(resid = numeric(300)
    , phylo_X = numeric(300)
    , phylo_int = numeric(300)
    , phylo_cor = numeric(300)
    , B0 = numeric(300)
    , B1 = numeric(300)
    , model = numeric(300)
    , time = numeric(300)
    , convergence = numeric(300)
  )
  for(i in 1:length(tt)){
    MCMCglmm_obj <- readRDS(paste(MCMCglmm_path,tt[i],sep=""))
    var_dat <- as.data.frame(MCMCglmm_obj[[1]]$VCV)
    MCMCglmm_df[i,"resid"] <- median(var_dat[,"units"])
    MCMCglmm_df[i, "phylo_X"] <- median(var_dat[,"X:X.sp"])
    MCMCglmm_df[i, "phylo_int"] <- median(var_dat[,"(Intercept):(Intercept).sp"])
    MCMCglmm_df[i,"phylo_cor"] <- median(var_dat[,"(Intercept):X.sp"])
    b_dat <- as.data.frame(MCMCglmm_obj[[1]]$Sol)
    MCMCglmm_df[i, "B0"] <- as.numeric(between(0
            , quantile(b_dat[,"(Intercept)"], 0.025)
            , quantile(b_dat[,"(Intercept)"], 0.975)
          ))
    MCMCglmm_df[i,"B1"] <- as.numeric(between(0
            , quantile(b_dat[,"X"], 0.025)
            , quantile(b_dat[,"X"], 0.975)
          ))
    MCMCglmm_df[i,"model"] <- tt[i]
    MCMCglmm_df[i,"time"] <- MCMCglmm_obj[[2]][[1]]
  }
  return(MCMCglmm_df)
}

MCMCglmmss_data <- MCMCglmmss_results(MCMCglmmss_res)

#### Collect phylolm single site results ----
## need to check out NA cases

phylolm_path <- "./datadir/phylolm/"
phylolm_res <- list.files(path = phylolm_path)
phylolm_results <- function(tt){
  phylolm_df <- data.frame(resid = numeric(400)
    , phylo_int = numeric(400)
    , phylo_X = NA
    , phylo_cor = NA
    , B0 = numeric(400)
    , B1 = numeric(400)
    , model = numeric(400)
    , time = numeric(400)
    , convergence = NA
  )
  for(i in 1:length(tt)){
    phylolm_obj <- readRDS(paste(phylolm_path,tt[i],sep=""))
    phylolm_df[i,"resid"] <- as.numeric(phylolm_obj[[1]]["sigma2_error"])
    phylolm_df[i,"phylo_int"] <- as.numeric(phylolm_obj[[1]]["sigma2"])
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

ssdat <- rbind(gls_data, phylolm_data, lme4ss_data, glmmTMBss_data, brmsss_data, MCMCglmmss_data)


### Collect multiple sites ----

lme4_path <- "./datadir/lme4/"
lme4ms_res <- list.files(path = lme4_path, pattern = "ms")
lme4ms_results <- function(tt){
	lme4ms_df <- data.frame(resid = numeric(length(tt))
		, phylo_X = NA
		, phylo_int = NA
		, phylo_cor = NA
		, phylo_interaction = NA
		, species_X = NA
		, species_int = NA
		, species_cor = NA
		, site_int = NA
		, B0 = NA
		, B1 = NA
		, model = NA
		, time = NA
		, convergence = NA
		)
	for(i in 1:length(tt)){
	  lme4_obj <- readRDS(paste(lme4_path,tt[i],sep=""))
	  covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
	  lme4ms_df[i,"resid"] <- (covmat 
      %>% filter(grp=="Residual") 
	    %>% dplyr:::select(vcov) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i, "phylo_X"] <- (covmat 
	    %>% filter((grp=="sp") & (var1 =="X")) 
	    %>% dplyr:::select(vcov) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i, "phylo_int"] <- (covmat 
	    %>% filter((grp=="sp") 
	      & (var1 == "(Intercept)")
	      & (is.na(var2))
	      ) 
	    %>% dplyr:::select(vcov) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i,"phylo_cor"] <- (covmat 
	    %>% filter((grp=="sp") 
	      & (var1 == "(Intercept)")
	      & (var2 == "X")
	      ) 
	    %>% dplyr:::select(vcov) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i,"phylo_interaction"] <- (covmat
	    %>% filter((grp=="sp:site"))
	    %>% dplyr:::select(vcov)
	    %>% as.numeric()
	  )
	  lme4ms_df[i, "species_X"] <- (covmat 
	    %>% filter((grp=="obs") & (var1 =="X"))
	    %>% dplyr:::select(vcov)
	    %>% as.numeric()
	  )
	  lme4ms_df[i, "species_int"] <- (covmat 
	    %>% filter((grp=="obs") 
	      & (var1 == "(Intercept)")
	      & (is.na(var2))
	      ) 
	    %>% dplyr:::select(vcov) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i,"species_cor"] <- (covmat 
	    %>% filter((grp=="obs") 
	      & (var1 == "(Intercept)")
	      & (var2 == "X")
	      ) 
	    %>% dplyr:::select(vcov) 
	    %>% as.numeric()
	  )
	  lme4ms_df[i,"site_int"] <- (covmat 
	    %>% filter(grp=="site")
	    %>% dplyr:::select(vcov) 
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
  pez_df <- data.frame(resid = numeric(length(tt))
    , phylo_X = NA
    , phylo_int = NA
    , phylo_cor = NA
    , phylo_interaction = NA
    , species_X = NA
    , species_int = NA
    , species_cor = NA
    , site_int = NA
    , B0 = NA
    , B1 = NA
    , model = NA
    , time = NA
    , convergence = NA
  )
  for(i in 1:length(tt)){
    pez_obj <- readRDS(paste(pez_path,tt[i],sep=""))
    pez_df[i,"resid"] <- unlist(pez_obj[[2]]["s2resid"])
	  pez_df[i,"phylo_interaction"] <- unlist(pez_obj[[2]]["s2n"])
    pez_df[i,"phylo_int"] <- unlist(pez_obj[[2]]["s2r"])[1]
    pez_df[i,"phylo_X"] <- unlist(pez_obj[[2]]["s2r"])[2]
    pez_df[i,"species_int"] <- unlist(pez_obj[[2]]["s2r"])[3]
    pez_df[i,"species_X"] <- unlist(pez_obj[[2]]["s2r"])[4]
    pez_df[i,"site_int"] <- unlist(pez_obj[[2]]["s2r"])[5]
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
  phyr_df <- data.frame(resid = numeric(length(tt))
                       , phylo_X = NA
                       , phylo_int = NA
                       , phylo_cor = NA
                       , phylo_interaction = NA
                       , species_X = NA
                       , species_int = NA
                       , species_cor = NA
                       , site_int = NA
                       , B0 = NA
                       , B1 = NA
                       , model = NA
                       , time = NA
                       , convergence = NA
  )
  for(i in 1:length(tt)){
    phyr_obj <- readRDS(paste(phyr_path,tt[i],sep=""))
    phyr_df[i,"resid"] <- unlist(phyr_obj[[1]]["s2resid"])
    phyr_df[i,"phylo_interaction"] <- unlist(phyr_obj[[1]]["s2n"])
    phyr_df[i,"phylo_int"] <- unlist(phyr_obj[[1]]["s2r"])[2]
    phyr_df[i,"phylo_X"] <- unlist(phyr_obj[[1]]["s2r"])[4]
    phyr_df[i,"species_int"] <- unlist(phyr_obj[[1]]["s2r"])[1]
    phyr_df[i,"species_X"] <- unlist(phyr_obj[[1]]["s2r"])[3]
    phyr_df[i,"site_int"] <- unlist(phyr_obj[[1]]["s2r"])[5]
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





glmmTMB_path <- "./datadir/glmmTMB/"
glmmTMB_res <- list.files(path = glmmTMB_path, pattern = "ms")
glmmTMBms_results <- function(tt){
  glmmTMBms_df <- data.frame(resid = length(tt)
    , phylo_X = NA
    , phylo_int = NA
    , phylo_cor = NA
    , phylo_interaction = NA
    , species_X = NA
    , species_int = NA
    , species_cor = NA
    , site_int = NA
    , B0 = NA
    , B1 = NA
    , model = NA
    , time = NA
    , convergence = NA
  )
  for(i in 1:length(tt)){
    print(i)
    glmmTMB_obj <- readRDS(paste(glmmTMB_path,tt[i],sep=""))
    covobj <- VarCorr(glmmTMB_obj[[1]])[["cond"]]
    glmmTMBms_df[i,"resid"] <- attr(covobj,"sc")
    glmmTMBms_df[i, "phylo_X"] <- covobj[["sp"]]["X","X"]
    glmmTMBms_df[i, "phylo_int"] <- covobj[["sp"]]["(Intercept)","(Intercept)"]
    glmmTMBms_df[i, "phylo_cor"] <- covobj[["sp"]]["(Intercept)","X"]
    glmmTMBms_df[i,"phylo_interaction"] <- covobj[["sp:site"]]["(Intercept)","(Intercept)"]
    glmmTMBms_df[i, "species_X"] <- covobj[["obs"]]["X","X"]
    glmmTMBms_df[i, "species_int"] <- covobj[["obs"]]["(Intercept)","(Intercept)"]
    glmmTMBms_df[i, "species_cor"] <- covobj[["obs"]]["(Intercept)","X"]
    glmmTMBms_df[i, "site_int"] <- covobj[["site"]]["(Intercept)","(Intercept)"]
    B0 <- coef(summary(glmmTMB_obj[[1]]))[["cond"]]["(Intercept)","Estimate"]
    B0se <- coef(summary(glmmTMB_obj[[1]]))[["cond"]]["(Intercept)","Std. Error"]
    B1 <- coef(summary(glmmTMB_obj[[1]]))[["cond"]]["X","Estimate"]
    B1se <- coef(summary(glmmTMB_obj[[1]]))[["cond"]]["X","Std. Error"]
    glmmTMBms_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
    glmmTMBms_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
    glmmTMBms_df[i,"model"] <- tt[i]
    glmmTMBms_df[i,"time"] <- glmmTMB_obj[[2]][[1]]
    glmmTMBms_df[i,"convergence"] <- glmmTMB_obj[[1]]$fit$convergence
    
  }
  return(glmmTMBms_df)
}

glmmTMBms_data <- glmmTMBms_results(glmmTMB_res)


msdat <- rbind(lme4ms_data,pez_data, phyr_data, glmmTMBms_data)
#### Save results ----
data_list <- list(ssdat,msdat)

# data_list <- list(gls_data, lme4ss_data, lme4ss_slope_data, lme4ms_data, pez_data, lme4ms_slope_data , pez_slope_data, lme4cs_data, pez_cs_data, lme4cs_slope_data, pez_cs_slope_data)
saveRDS(data_list,file="./datadir/collect_rerun.RDS")
# rdsave(ssdat, msdat)


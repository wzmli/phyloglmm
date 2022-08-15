#### Collect simulated results ----
## load packages
library(lme4)
library(dplyr)
library(ggplot2)
library(tidyr)
library(glmmTMB)
library(broom)
library(broom.mixed)

verbose <- TRUE
target_effsize <- 4000

ss_names = c("resid", "phylo_int", "phylo_X", "phylo_cor",
             "B0", "B1", "model", "rawtime", "time", "convergence")
ms_names <- c(ss_names, "phylo_interaction",
              "species_X", "species_int", "species_cor",
              "site_int")

blank_df <- function(nrow, type="ss") {
  colnames <- if (type=="ss") ss_names else ms_names
  df <- as.data.frame(matrix(NA, nrow, length(colnames)))
  names(df) <- colnames
  return(df)
}

#### Collect gls results ----
if (verbose) cat("gls\n")


if (verbose) cat("lme4ms\n")
lme4_path <- "./datadir/lme4/"
lme4ms_res <- list.files(path = lme4_path, pattern = "[.]ms[.]")
lme4mms_res <- list.files(path = lme4_path, pattern = "mms")
lme4ms_results <- function(tt){
	lme4ms_df <- blank_df(length(tt), type="ms")
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
	  lme4ms_df[i,"time"] <- lme4_obj[[2]][["elapsed"]]
	  lme4ms_df[i,"convergence"] <- lme4_obj[[1]]@optinfo$conv$opt
	}
	return(lme4ms_df)
}

lme4ms_data <- lme4ms_results(lme4ms_res)
lme4mms_data <- lme4ms_results(lme4mms_res)
# lme4mms_data <- lme4mms_data %>% 

if (verbose) cat("glmmTMBms\n")
glmmTMB_path <- "./datadir/glmmTMB/"
glmmTMB_res <- list.files(path = glmmTMB_path, pattern = "[.]ms[.]")
glmmTMBmms_res <- list.files(path = glmmTMB_path, pattern = "mms")

glmmTMBms_results <- function(tt){
  glmmTMBms_df <- blank_df(length(tt), type = "ms")
  for(i in 1:length(tt)){
    print(i)
    glmmTMB_obj <- readRDS(paste(glmmTMB_path,tt[i],sep=""))
    covobj <- VarCorr(glmmTMB_obj[[1]])[["cond"]]
    glmmTMBms_df[i,"resid"] <- attr(covobj,"sc")^2
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
    glmmTMBms_df[i,"time"] <- glmmTMB_obj[[2]][["elapsed"]]
    glmmTMBms_df[i,"convergence"] <- glmmTMB_obj[[1]]$fit$convergence

  }
  return(glmmTMBms_df)
}

glmmTMBms_data <- glmmTMBms_results(glmmTMB_res)
glmmTMBmms_data <- glmmTMBms_results(glmmTMBmms_res)


msdat <- rbind(lme4ms_data,glmmTMBms_data
, lme4mms_data,glmmTMBmms_data)

saveRDS(msdat, file="./datadir/ms_lme4verse.RDS")

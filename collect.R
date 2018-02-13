#### Collect simulated results ----
## load packages
library(pez)
library(lme4)
library(dplyr)
library(ggplot2)
library(tidyr)

#### Collect gls results ----

gls_path <- "./datadir/gls/"
gls_res <- list.files(path = gls_path)
gls_results <- function(tt){
  gls_list <- list()
  gls_df <- data.frame(resid = numeric(200)
    , model = numeric(200)
    , time = numeric(200)
  )
  for(i in tt){
    gls_list[[i]] <- gls_df
    ff <- list.files(path=paste(gls_path,i,"/",sep=""),patter="rds")
    for(j in 1:length(ff)){
      gls_obj <- readRDS(paste(gls_path,i,"/",ff[j],sep=""))
      gls_list[[i]][j,"resid"] <- attr(resid(gls_obj[[1]]),"std")[[1]]
      gls_list[[i]][j,"model"] <- i
      gls_list[[i]][j,"time"] <- gls_obj[[2]][[1]]
    }
  }
  return(gls_list)
}

gls_data <- bind_rows(gls_results(gls_res))

#### Collect lme4 single site results ----

lme4ss_path <- "./datadir/lme4/ss/"
lme4ss_res <- list.files(path = lme4ss_path, pattern = "resid")
lme4ss_results <- function(pp,tt){
  lme4ss_list <- list()
  lme4ss_df <- data.frame(resid = numeric(200)
    , phylo = numeric(200)
    , model = numeric(200)
    , time = numeric(200)
  )
  for(i in tt){
    lme4ss_list[[i]] <- lme4ss_df
    ff <- list.files(path=paste(pp,i,"/",sep=""),patter="rds")
    for(j in 1:length(ff)){
      lme4_obj <- readRDS(paste(pp,i,"/",ff[j],sep=""))
      covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
      lme4ss_list[[i]][j,"resid"] <- covmat %>% filter(grp=="Residual") %>% select(sdcor) %>% as.numeric()
      lme4ss_list[[i]][j,"phylo"] <- covmat %>% filter(grp=="sp") %>% select(sdcor) %>% as.numeric()
      lme4ss_list[[i]][j,"model"] <- i
      lme4ss_list[[i]][j,"time"] <- lme4_obj[[2]][[1]]
    }
  }
  return(lme4ss_list)
}

lme4ss_data <- bind_rows(lme4ss_results(lme4ss_path,lme4ss_res))

#### Collect lme4 single site slope results ----
## need to check out NA cases

lme4ss_slope_path <- "./datadir/lme4/ss/slope/"
lme4ss_slope_res <- list.files(path = lme4ss_slope_path, pattern = "cor")
lme4ss_slope_results <- function(pp,tt){
  lme4ss_list <- list()
  lme4ss_df <- data.frame(resid = numeric(200)
    , phylo = numeric(200)
    , phyloX = numeric(200)
    , cor = numeric(200)
    , model = numeric(200)
    , time = numeric(200)
  )
  for(i in tt){
    lme4ss_list[[i]] <- lme4ss_df
    ff <- list.files(path=paste(pp,i,"/",sep=""),patter="rds")
    for(j in 1:length(ff)){
      lme4_obj <- readRDS(paste(pp,i,"/",ff[j],sep=""))
      covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
      lme4ss_list[[i]][j,"resid"] <- covmat %>% filter(grp == "Residual") %>% select(sdcor) %>% as.numeric()
      lme4ss_list[[i]][j,"phylo"] <- covmat %>% filter(var1 == "(Intercept)") %>% filter(is.na(var2)) %>% select(sdcor) %>% as.numeric
      lme4ss_list[[i]][j,"phyloX"] <- covmat %>% filter(var1 == "X") %>% select(sdcor) %>% as.numeric
      lme4ss_list[[i]][j,"cor"] <- covmat %>% filter(var1 == "(Intercept)") %>% filter(var2 == "X") %>% select(sdcor) %>% as.numeric
      lme4ss_list[[i]][j,"model"] <- i
      lme4ss_list[[i]][j,"time"] <- lme4_obj[[2]][[1]]
    }
  }
  return(lme4ss_list)
}

lme4ss_slope_data <- bind_rows(lme4ss_slope_results(lme4ss_slope_path, lme4ss_slope_res))

#### Collect lme4 multiple site results ----

lme4ms_path <- "./datadir/lme4/ms/mm_lme4/"
lme4ms_res <- list.files(path = lme4ms_path, pattern = "lme4")

lme4ms_data <- bind_rows(lme4ss_results(lme4ms_path,lme4ms_res))

#### Collect pez ms ----


pez_path <- "./datadir/pez/"
pez_res <- c("pez_small","pez_med")
pez_results <- function(pp,tt){
  pez_list <- list()
  pez_df <- data.frame(resid = numeric(200)
    , phylo = numeric(200)
    , phyloX = numeric(200)
    , model = numeric(200)
    , time = numeric(200)
  )
  for(i in tt){
    pez_list[[i]] <- pez_df
    ff <- list.files(path=paste(pp,i,"/",sep=""),patter="rds")
    for(j in 1:length(ff)){
      pez_obj <- readRDS(paste(pp,i,"/",ff[j],sep=""))
      pez_list[[i]][j,"resid"] <- sqrt(pez_obj[[1]]$s2resid)
      pez_list[[i]][j,"phylo"] <- sqrt(pez_obj[[1]]$s2r[1])
      pez_list[[i]][j,"phyloX"] <- sqrt(pez_obj[[1]]$s2r[2])
      pez_list[[i]][j,"model"] <- i
      pez_list[[i]][j,"time"] <- pez_obj[[2]][[1]]
    }
  }
  return(pez_list)
}

pez_data <- bind_rows(pez_results(pez_path,pez_res))

#### Collect lme4 multiple site slope results ----

lme4ms_slope_path <- "./datadir/lme4/ms/mm_lme4/mm_slope/"
lme4ms_slope_res <- list.files(path = lme4ms_slope_path, pattern = "cor")

lme4ms_slope_data <- bind_rows(lme4ss_slope_results(lme4ms_slope_path, lme4ms_slope_res))

#### Collect pez slope results ----

pez_slope_path <- "./datadir/pez/slope/"
pez_slope_res <- c("small_uncor","med_uncor")
pez_slope_data <- bind_rows(pez_results(pez_slope_path,pez_slope_res))

#### Collect lme4 cs results ----

lme4cs_path <- "./datadir/lme4/ms/lme4_cs/"
lme4cs_res <- list.files(path = lme4cs_path, pattern = "cs")
lme4cs_results <- function(pp,tt){
  lme4cs_list <- list()
  lme4cs_df <- data.frame(resid = numeric(200)
    , phylo = numeric(200)
    , sp_by_site = numeric(200)
    , model = numeric(200)
    , time = numeric(200)
  )
  for(i in tt){
    lme4cs_list[[i]] <- lme4cs_df
    ff <- list.files(path=paste(pp,i,"/",sep=""),patter="rds")
    for(j in 1:length(ff)){
      lme4_obj <- readRDS(paste(pp,i,"/",ff[j],sep=""))
      covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
      lme4cs_list[[i]][j,"resid"] <- covmat %>% filter(grp=="Residual") %>% select(sdcor) %>% as.numeric()
      lme4cs_list[[i]][j,"phylo"] <- covmat %>% filter(grp=="sp") %>% select(sdcor) %>% as.numeric()
      lme4cs_list[[i]][j,"sp_by_site"] <- covmat %>% filter(grp=="sp:site") %>% select(sdcor) %>% as.numeric()
      lme4cs_list[[i]][j,"model"] <- i
      lme4cs_list[[i]][j,"time"] <- lme4_obj[[2]][[1]]
    }
  }
  return(lme4cs_list)
}

#lme4cs_data <- bind_rows(lme4cs_results(lme4cs_path,lme4cs_res))

#### Collect pez cs results ----

pez_cs_path <- "./datadir/pez/pez_cs/"
pez_res <- c("pez_small","pez_med")
pez_cs_results <- function(pp,tt){
  pez_list <- list()
  pez_df <- data.frame(resid = numeric(200)
                       , phylo = numeric(200)
                       , phyloX = numeric(200)
                       , sp_by_site = numeric(200)
                       , model = numeric(200)
                       , time = numeric(200)
  )
  for(i in tt){
    pez_list[[i]] <- pez_df
    ff <- list.files(path=paste(pp,i,"/",sep=""),patter="rds")
    for(j in 1:length(ff)){
      pez_obj <- readRDS(paste(pp,i,"/",ff[j],sep=""))
      pez_list[[i]][j,"resid"] <- sqrt(pez_obj[[2]]$s2resid)
      pez_list[[i]][j,"phylo"] <- sqrt(pez_obj[[2]]$s2r[1])
      pez_list[[i]][j,"phyloX"] <- sqrt(pez_obj[[2]]$s2r[2])
      pez_list[[i]][j,"sp_by_site"] <- sqrt(pez_obj[[2]]$s2n[1])
      pez_list[[i]][j,"model"] <- i
      pez_list[[i]][j,"time"] <- pez_obj[[1]][[1]]
    }
  }
  return(pez_list)
}

#pez_cs_data <- bind_rows(pez_cs_results(pez_cs_path,pez_res))

#### lme4 cs slopes NOT DONE, MISSING uncor ----

lme4cs_slope_path <- "./datadir/lme4/ms/lme4_cs/slope/"
lme4cs_slopecor_res <- c("large_cor","med_cor","small_cor")
lme4cs_slopeuncor_res <- c("large_uncor","med_uncor","small_uncor")

lme4cs_slopecor_results <- function(pp,tt){
  lme4cs_list <- list()
  lme4cs_df <- data.frame(resid = numeric(200)
    , phylo = numeric(200)
    , phyloX = numeric(200)
    , cor = numeric(200)
    , sp_by_site = numeric(200)
    , model = numeric(200)
    , time = numeric(200)
  )
  for(i in tt){
    lme4cs_list[[i]] <- lme4cs_df
    ff <- list.files(path=paste(pp,i,"/",sep=""),patter="rds")
    for(j in 1:length(ff)){
      lme4_obj <- readRDS(paste(pp,i,"/",ff[j],sep=""))
      covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
      lme4cs_list[[i]][j,"resid"] <- covmat %>% filter(grp == "Residual") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"phylo"] <- covmat %>% filter(var1 == "(Intercept)") %>% filter(is.na(var2)) %>% filter(grp=="sp") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"phyloX"] <- covmat %>% filter(var1 == "X") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"cor"] <- covmat %>% filter(var1 == "(Intercept)") %>% filter(var2 == "X") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"sp:site"] <- covmat %>% filter(grp == "sp:site") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"model"] <- i
      lme4cs_list[[i]][j,"time"] <- lme4_obj[[2]][[1]]
    }
  }
  return(lme4cs_list)
}

lme4cs_slopeuncor_results <- function(pp,tt){
  lme4cs_list <- list()
  lme4cs_df <- data.frame(resid = numeric(200)
                          , phylo = numeric(200)
                          , phyloX = numeric(200)
                          , cor = numeric(200)
                          , sp_by_site = numeric(200)
                          , model = numeric(200)
                          , time = numeric(200)
  )
  for(i in tt){
    lme4cs_list[[i]] <- lme4cs_df
    ff <- list.files(path=paste(pp,i,"/",sep=""),patter="rds")
    for(j in 1:length(ff)){
      lme4_obj <- readRDS(paste(pp,i,"/",ff[j],sep=""))
      covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[1]]))
      lme4cs_list[[i]][j,"resid"] <- covmat %>% filter(grp == "Residual") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"phylo"] <- covmat %>% filter(var1 == "(Intercept)") %>% filter(grp == "sp.1") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"phyloX"] <- covmat %>% filter(var1 == "X") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"cor"] <- covmat %>% filter(var1 == "(Intercept)") %>% filter(var2 == "X") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"sp:site"] <- covmat %>% filter(grp == "sp.site") %>% select(sdcor) %>% as.numeric
      lme4cs_list[[i]][j,"model"] <- i
      lme4cs_list[[i]][j,"time"] <- lme4_obj[[2]][[1]]
    }
  }
  return(lme4cs_list)
}

#lme4cs_slopecor_data <- bind_rows(lme4cs_slopecor_results(lme4cs_slope_path,lme4cs_slopecor_res))
#lme4cs_slopeuncor_data <- bind_rows(lme4cs_slopeuncor_results(lme4cs_slope_path,lme4cs_slopeuncor_res))

#lme4cs_slope_data <- rbind(lme4cs_slopecor_data,lme4cs_slopeuncor_data)
#### Collect pez cs slope ----
## phylo and phyloX out of order?
pez_cs_slope_path <- "./datadir/pez/pez_cs/slope/"

pez_cs_slope_res <- list.files(pez_cs_slope_path,pattern="uncor")

#pez_cs_slope_data <- bind_rows(pez_cs_results(pez_cs_slope_path,pez_cs_slope_res))


#### Save results ----

data_list <- list(gls_data, lme4ss_data, lme4ss_slope_data, lme4ms_data, pez_data, lme4ms_slope_data , pez_slope_data) #, lme4cs_data, pez_cs_data, lme4cs_slope_data, pez_cs_slope_data)
saveRDS(data_list,file="./datadir/result_list.RDS")
# rdsave(gls_data, lme4ss_data, lme4ss_slope_data, lme4ms_data, pez_data, lme4ms_slope_data , pez_slope_data)
#,lme4cs_data, pez_cs_data, lme4cs_slope_data, pez_cs_slope_data)

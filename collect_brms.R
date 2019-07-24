library(brms)
library(dplyr)

brms_path <- "./datadir/brms/"
brmsss_res <- list.files(path = brms_path, pattern = "ss")
brmsss_results <- function(tt){
  brms_df <- data.frame(resid = length(tt)
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

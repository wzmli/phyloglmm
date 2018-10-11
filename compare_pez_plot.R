#### Collect simulated results ----
## load packages
library(pez)
library(lme4)
library(dplyr)
library(ggplot2)
library(brms)
library(tidyr)
library(phylolm)
library(cowplot)

lme4_path <- "./datadir/"
lme4ms_res <- list.files(path = lme4_path, pattern = "rds")
lme4ms_results <- function(tt){
  lme4ms_df <- data.frame(resid = numeric(10)
                          , phylo_X = numeric(10)
                          , phylo_int = numeric(10)
                          , phylo_cor = numeric(10)
                          , phylo_interaction = numeric(10)
                          , species_X = numeric(10)
                          , species_int = numeric(10)
                          , species_cor = numeric(10)
                          , site_int = numeric(10)
                          , B0 = numeric(10)
                          , B1 = numeric(10)
                          , model = numeric(10)
                          , platform = "lme4"
  )
  for(i in 1:length(tt)){
    lme4_obj <- readRDS(paste(lme4_path,tt[i],sep=""))
    covmat <- as.data.frame(lme4::VarCorr(lme4_obj[[2]]))
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
                                  %>% filter((grp=="sp.1") 
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
                                         %>% filter((grp=="sp.site"))
                                         %>% select(sdcor)
                                         %>% as.numeric()
    )
    lme4ms_df[i, "species_X"] <- (covmat 
                                  %>% filter((grp=="obs") & (var1 =="X"))
                                  %>% select(sdcor)
                                  %>% as.numeric()
    )
    lme4ms_df[i, "species_int"] <- (covmat 
                                    %>% filter((grp=="obs.1") 
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
    lme4ms_df[i,"site_int"] <- (covmat 
                                %>% filter(grp=="site")
                                %>% select(sdcor) 
                                %>% as.numeric()
    )
    B0 <- coef(summary(lme4_obj[[2]]))["(Intercept)","Estimate"]
    B0se <- coef(summary(lme4_obj[[2]]))["(Intercept)","Std. Error"]
    B1 <- coef(summary(lme4_obj[[2]]))["X","Estimate"]
    B1se <- coef(summary(lme4_obj[[2]]))["X","Std. Error"]
    lme4ms_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
    lme4ms_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
    lme4ms_df[i,"model"] <- tt[i]
  }
  return(lme4ms_df)
}

lme4ms_data <- lme4ms_results(lme4ms_res)



pez_path <- "./datadir/"
pez_res <- list.files(path = pez_path, pattern = "rds")
pez_results <- function(tt){
  pez_df <- data.frame(resid = numeric(10)
                       , phylo_X = numeric(10)
                       , phylo_int = numeric(10)
                       , phylo_cor = NA
                       , phylo_interaction = numeric(10)
                       , species_X = numeric(10)
                       , species_int = numeric(10)
                       , species_cor = NA
                       , site_int = numeric(10)
                       , B0 = numeric(10)
                       , B1 = numeric(10)
                       , model = numeric(10)
                       , convcode = numeric(10)
                       , platform = "pez"
  )
  for(i in 1:length(tt)){
    pez_obj <- readRDS(paste(pez_path,tt[i],sep=""))
    pez_df[i,"resid"] <- sqrt(unlist(pez_obj[[1]]["s2resid"]))
    pez_df[i,"phylo_interaction"] <- sqrt(unlist(pez_obj[[1]]["s2n"]))
    pez_df[i,"phylo_int"] <- sqrt(unlist(pez_obj[[1]]["s2r"]))[1]
    pez_df[i,"phylo_X"] <- sqrt(unlist(pez_obj[[1]]["s2r"]))[2]
    pez_df[i,"species_int"] <- sqrt(unlist(pez_obj[[1]]["s2r"]))[3]
    pez_df[i,"species_X"] <- sqrt(unlist(pez_obj[[1]]["s2r"]))[4]
    pez_df[i,"site_int"] <- sqrt(unlist(pez_obj[[1]]["s2r"]))[5]
    pez_df[i,"convcode"] <- pez_obj[[1]]["convcode"]
    B0 <- unlist(pez_obj[[1]]["B"])[1]
    B0se <- unlist(pez_obj[[1]]["B.se"])[1]
    B1 <- unlist(pez_obj[[1]]["B"])[2]
    B1se <- unlist(pez_obj[[1]]["B.se"])[2]
    pez_df[i,"B0"] <- as.numeric(between(0, B0-1.96*B0se, B0+1.96*B0se))
    pez_df[i,"B1"] <- as.numeric(between(0, B1-1.96*B1se, B1+1.96*B1se))
    pez_df[i,"model"] <- tt[i]
  }
  return(pez_df)
}

pez_data <- pez_results(pez_res)
# 
# msdat_raw <- rbind(lme4ms_data,pez_data)

msdat <- (lme4ms_data 
	%>% left_join(.,pez_data,by = "model")
	%>% separate(model, c("platform", "sites", "size", "seed", "saveformat"), "[.]") 
)

gg_temp <- (ggplot(msdat, aes(colour=factor(convcode)))
  + geom_abline(slope = 1, intercept = 0)
  + theme_bw()
)

gg_phylo_int <- gg_temp + geom_point(aes(x=phylo_int.x, y=phylo_int.y))
gg_phylo_X <- gg_temp + geom_point(aes(x=phylo_X.x, y=phylo_X.y))

gg_species_int <- gg_temp + geom_point(aes(x=species_int.x, y=species_int.y))
gg_species_X <- gg_temp + geom_point(aes(x=species_X.x, y=species_X.y))


gg_site_int <- gg_temp + geom_point(aes(x=site_int.x, y=site_int.y))
gg_interaction <- gg_temp + geom_point(aes(x=phylo_interaction.x, y=phylo_interaction.y))
gg_resid <- gg_temp + geom_point(aes(x=resid.x, y=resid.y))

print(gg_site_int)

print(plot_grid(gg_phylo_int, gg_species_int, gg_phylo_X, gg_species_X, gg_interaction, gg_resid, ncol=2 ))


head(lme4ms_data)
rlme4 <-(lme4ms_data
  %>% select(-c(model, platform, phylo_cor, species_cor))
  %>% gather(key = type, value = sd)
  %>% filter(sd > 0.0)
)

print(ggplot(rlme4, aes(y=sd, x=type))+geom_violin()+coord_flip())

print(rlme4 %>% filter(type %in% c("B1","B0")) %>% summary())

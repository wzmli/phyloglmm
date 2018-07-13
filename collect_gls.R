library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)

obs <- 50

results_dir_glslist <- list.files(path="test",pattern="RIonly")

glsdf <- data.frame(phylo.int = numeric(obs)
  , model = numeric(obs)
  , time = numeric(obs)
)

df_glslist <- list(glsdf,glsdf,glsdf) #list(glsdf,glsdf,glsdf,glsdf,glsdf)

for(i in 1:length(results_dir_glslist)){
  results_dir <- results_dir_glslist[[i]]
  ff <- list.files(path=paste("test/",results_dir,"/",sep=""),pattern="rds")
  for(j in 1:length(ff)){
    gls_list <- readRDS(paste("test/",results_dir,"/",ff[[j]],sep=""))
    df_glslist[[i]][j,1] <- attr(resid(gls_list[[1]]),"std")[[1]]
    df_glslist[[i]][j,2] <- results_dir
    df_glslist[[i]][j,3] <- gls_list[[2]][[1]]
  }
}

gls_df <- bind_rows(df_glslist)

gg <- (ggplot(gls_df,aes(y=phylo.int,x=model,group=model))
  + geom_violin()
  + geom_hline(yintercept = 4)
  + theme_bw()
  + ylab("phylo sd")
  )

print(gg)


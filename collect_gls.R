library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)

results_dir_glslist <- list.files(path="datadir",pattern="gls")
results_dir_lme4list <- list.files(path="datadir",pattern="lme4")

glsdf <- data.frame(phylo.int = numeric(200)
  , model = numeric(200)
  , time = numeric(200)
)

lme4df <- data.frame(phylo.int = numeric(200)
                    , residual = numeric(200)
                    , model = numeric(200)
                    , time = numeric(200)
)
df_glslist <- list(glsdf,glsdf,glsdf,glsdf,glsdf)
df_lme4list <- list(lme4df,lme4df,lme4df,lme4df,lme4df)

for(i in 1:length(results_dir_glslist)){
  results_dir <- results_dir_glslist[[i]]
  ff <- list.files(path=paste("datadir/",results_dir,"/",sep=""),pattern="rds")
  for(j in 1:length(ff)){
    gls_list <- readRDS(paste("datadir/",results_dir,"/",ff[[j]],sep=""))
    df_glslist[[i]][j,1] <- attr(resid(gls_list[[1]]),"std")[[1]]
    df_glslist[[i]][j,2] <- results_dir
    df_glslist[[i]][j,3] <- gls_list[[2]][[1]]
  }
}

for(i in 1:length(results_dir_lme4list)){
  results_dir <- results_dir_lme4list[[i]]
  ff <- list.files(path=paste("datadir/",results_dir,"/",sep=""),pattern="rds")
  for(j in 1:length(ff)){
    lme4_list <- readRDS(paste("datadir/",results_dir,"/",ff[[j]],sep=""))
    vardf <- as.data.frame(VarCorr(lme4_list[[1]]))
    df_lme4list[[i]][j,1] <- vardf$sdcor[1]
    df_lme4list[[i]][j,2] <- vardf$sdcor[2]
    df_lme4list[[i]][j,3] <- results_dir
    df_lme4list[[i]][j,4] <- lme4_list[[2]][[1]]
  }
}

gls_df <- bind_rows(df_glslist)
lme4_df <- bind_rows(df_lme4list)

gg <- (ggplot(gls_df,aes(y=phylo.int,x=model,group=model))
  + geom_violin()
  + geom_hline(yintercept = 4)
  + theme_bw()
  + ylab("phylo sd")
  )

print(gg)

gglme4 <- (ggplot(lme4_df,aes(y=phylo.int,x=model,group=model))
          + geom_violin()
          + geom_hline(yintercept = 4)
          + theme_bw()
          + ylab("phylo sd")
)

print(gglme4)

gglme4_2 <- (ggplot(lme4_df,aes(y=residual,x=model,group=model))
           + geom_violin()
           + geom_hline(yintercept = 10)
           + theme_bw()
           + ylab("Residual")
)

print(gglme4_2)

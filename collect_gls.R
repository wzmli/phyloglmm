library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)

obs <- 50

results_dir_glslist <- list.files(path="test",pattern="0")

glsdf <- data.frame(phylo.int = numeric(obs)
  , intercept = numeric(obs)
  , slope = numeric(obs)
  , model = numeric(obs)
  , time = numeric(obs)
)

df_glslist <- list(glsdf,glsdf,glsdf,glsdf,glsdf,glsdf,glsdf)

for(i in 1:length(results_dir_glslist)){
  results_dir <- results_dir_glslist[[i]]
  ff <- list.files(path=paste("test/",results_dir,"/",sep=""),pattern="rds")
  for(j in 1:length(ff)){
    gls_list <- readRDS(paste("test/",results_dir,"/",ff[[j]],sep=""))
    df_glslist[[i]][j,1] <- attr(resid(gls_list[[1]]),"std")[[1]]
    df_glslist[[i]][j,2] <- coef(gls_list[[1]])[1]
    df_glslist[[i]][j,3] <- coef(gls_list[[1]])[2]
    df_glslist[[i]][j,4] <- results_dir
    df_glslist[[i]][j,5] <- gls_list[[2]][[1]]
  }
}

gls_df <- (bind_rows(df_glslist)
     %>% select(-time)
#     %>% group_by(model)
     %>% gather(key="parameter",value="estimate",-model)
     %>% mutate(model = factor(model,levels=c("TIP50","TIP200","RI50","RI200","RI500","RI1000","RSI200")))
)

gg <- (ggplot(gls_df,aes(y=estimate,x=model,group=model))
  + geom_violin()
  + facet_wrap(~parameter,nrow=3)
#  + geom_hline(yintercept = 4)
  + theme_bw()
#  + ylab("sd")
  )

print(gg)


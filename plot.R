## plots for the ms

library(dplyr)
library(ggplot2)
library(tidyr)
theme_set(theme_bw())

gls_data2 <- (gls_data
  %>% transmute(sd=resid
    , sdtype = "phylo"
    , model = model
    , time = time)
  %>% separate(model,c("platform","size","type"),"_")
  %>% select(platform, type, size, sdtype, sd, time)
)

gg_gls <- (ggplot(data=gls_data2,aes(x=size,y=sd))
  + facet_grid(type~., scale="free_y")
  + geom_violin()
  + ggtitle("GLS")
)

# print(gg_gls)

lme4ss_data2 <- (lme4ss_data
  %>% gather(key=sdtype, value=sd, -c(time,model))
  %>% separate(model,c("platform","size","type"),"_")
  %>% select(platform, type, size, sdtype, sd, time)
)

gg_lme4ss <- (ggplot(data=lme4ss_data2,aes(x=size,y=sd))
  + facet_grid(type~sdtype, scale="free_y")
  + geom_violin()
  + ggtitle("lme4ss")
)

# print(gg_lme4ss)

ss_data <- rbind(gls_data2,lme4ss_data2)

gg_ss <- (ggplot(data=ss_data,aes(x=size,y=sd,fill=platform))
  + facet_grid(type~sdtype, scale="free_y")
  + geom_violin(position=position_dodge(width=0),alpha=0.4)
  + ggtitle("Single Site using GLS and lme4")
)

print(gg_ss)


## single site lme4 slope plots

ss_slope_data <- (lme4ss_slope_data
  %>% gather(key=sdtype, value=sd, -c(time,model))
  %>% separate(model,c("platform","size","type"),"_")
  %>% select(platform, type, size, sdtype, sd, time)
)

gg_ss_slope <- (ggplot(data=ss_slope_data,aes(x=size,y=sd))
  + facet_grid(sdtype~type, scale="free_y")
  + geom_violin()
  + ggtitle("LME4 Single site slope (Need to include/checkout NA removed cases")
)

print(gg_ss_slope)

## multiple sites

ms_data <- (pez_data
  %>% select(-phyloX)
  %>% rbind(.,lme4ms_data)
  %>% gather(key=sdtype, value=sd, -c(time,model))
  %>% separate(model,c("platform","size"),"_")
  %>% select(platform, size, sdtype, sd, time)
)

gg_ms <- (ggplot(data=ms_data,aes(x=size,y=sd,fill=platform))
  + facet_grid(.~sdtype, scale="free_y")
  + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
  + ggtitle("Multiple Site using pez and lme4")
)

print(gg_ms)

gg_mstime <- (ggplot(data=ms_data,aes(x=size,y=time,fill=platform))
  + geom_violin(position=position_dodge(width=0),alpha=0.4)
  + scale_y_log10()
  + ggtitle("Multiple Site timing")
)

print(gg_mstime)

## ms slope

ms_slope_data <- (pez_slope_data
  %>% mutate(cor = NA
        , model=paste("pez_",model,sep="")
      )
  %>% rbind(.,lme4ms_slope_data)
  %>% gather(key=sdtype, value=sd, -c(time,model))
  %>% separate(model,c("platform","size","type"),"_")
  %>% select(platform, size, sdtype, type, sd, time)
)

gg_ms_slope <- (ggplot(data=ms_slope_data,aes(x=size,y=sd,fill=platform))
  + facet_grid(type~sdtype, scale="free_y")
  + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
  + scale_y_log10()
  + ggtitle("Multiple Site using pez and lme4")
)

print(gg_ms_slope)

gg_ms_slope_time <- (ggplot(data=ms_slope_data,aes(x=size,y=time,fill=platform))
              + facet_grid(type~., scale="free_y")
              + geom_violin(position=position_dodge(width=0),alpha=0.4)
              + scale_y_log10()
              + ggtitle("Multiple Site slope timing")
)

print(gg_ms_slope_time)


## multiple sites cs


cs_data <- (pez_cs_data
  %>% select(-phyloX)
  %>% mutate(model=ifelse(model=="pez_small","pez_cs_small","pez_cs_med"))
  %>% rbind(.,lme4cs_data)
  %>% gather(key=sdtype, value=sd, -c(time,model))
  %>% separate(model,c("platform","type","size"),"_")
  %>% select(platform, size, sdtype, sd, time)
)

gg_cs <- (ggplot(data=cs_data,aes(x=size,y=sd,fill=platform))
  + facet_grid(.~sdtype, scale="free_y")
  + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
  + ggtitle("Compound symmetric using pez and lme4")
)

print(gg_cs)

gg_cstime <- (ggplot(data=cs_data,aes(x=size,y=time,fill=platform))
  + geom_violin(position=position_dodge(width=0),alpha=0.4)
  + scale_y_log10()
  + ggtitle("CS timing")
)

print(gg_cstime)

## cs slope

lme4cs_slope_data[,"sp_by_site"] <- lme4cs_slope_data[,"sp:site"]
lme4cs_slope_data <- (lme4cs_slope_data 
  %>% mutate(model = paste("lme4_",model,sep=""))
  %>% select(resid,phylo,phyloX,cor,sp_by_site,model,time)
)
cs_slope_data <- (pez_cs_slope_data
  %>% mutate(cor = NA
        , model=paste("pez_",model,sep="")
        , temp = phylo        ## hacking, I feel like pez gets the ordering confuse
        , phylo = phyloX
        , phyloX = temp
      )
  %>% select(-temp)
  %>% rbind(.,lme4cs_slope_data)
  %>% gather(key=sdtype, value=sd, -c(time,model))
  %>% separate(model,c("platform","size","type"),"_")
  %>% select(platform, size, sdtype, type, sd, time)
)

gg_cs_slope <- (ggplot(data=cs_slope_data,aes(x=size,y=sd,fill=platform))
                + facet_grid(type~sdtype, scale="free_y")
                + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
                # + scale_y_log10()
                + ggtitle("CS slope using pez and lme4")
)

print(gg_cs_slope)

gg_cs_slope_time <- (ggplot(data=cs_slope_data,aes(x=size,y=time,fill=platform))
                     + facet_grid(type~., scale="free_y")
                     + geom_violin(position=position_dodge(width=0),alpha=0.4)
                     + scale_y_log10()
                     + ggtitle("CS slope timing")
)

print(gg_cs_slope_time)


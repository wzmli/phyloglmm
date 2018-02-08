## plots for the ms

library(dplyr)
library(ggplot2)
library(tidyr)
theme_set(theme_bw())

data_list <- readRDS("./datadir/result_list.RDS")

gls_data <- data_list[[1]]
lme4ss_data <- data_list[[2]]
lme4ss_slope_data <- data_list[[3]]
lme4ms_data <- data_list[[4]]
pez_data <- data_list[[5]]
lme4ms_slope_data <- data_list[[6]]
pez_slope_data <- data_list[[7]]
#lme4cs_data <- data_list[[8]]
#pez_cs_data <- data_list[[9]]
#lme4cs_slope_data <- data_list[[10]]
#pez_cs_slope_data <- data_list[[11]]


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

ss_data <- (rbind(gls_data2,lme4ss_data2)
  %>% mutate(sdtype = ifelse(sdtype=="phylo","Phylogenetic","Tip")
    , type = ifelse(type=="noresid","Without Tip Var","With Tip Var")
    , Platform = platform
    , size = factor(size,levels=c("small","med","large"),labels=c("Small", "Medium", "Large"))
  )
)

ss_dummy <- data.frame(sdtype = factor(c("Phylogenetic","Tip","Phylogenetic","Tip"))
                       , type = factor(c("Without Tip Var","Without Tip Var","With Tip Var","With Tip Var"))
                       , Y=c(4,NA,4,10))
gg_ss <- (ggplot(data=ss_data,aes(x=size,y=sd,fill=Platform))
  + facet_grid(type~sdtype, scale="free_y")
  + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
  # + geom_hline(aes(yintercept = c(4)), linetype=c("solid"))
  # + geom_hline(aes(yintercept = c(10)), linetype=c("dashed"))
  + geom_hline(data=ss_dummy,aes(yintercept = Y))
  + scale_y_log10(limits=c(0.5,50),breaks=c(1,4,10,50))
  + ggtitle("Single Site using GLS and lme4")
  + ylab("Standard Deviation")
)

print(gg_ss)


## single site lme4 slope plots

ss_slope_data <- (lme4ss_slope_data
  %>% gather(key=sdtype, value=sd, -c(time,model))
  %>% separate(model,c("platform","size","type"),"_")
  %>% select(platform, type, size, sdtype, sd, time)
  %>% mutate(sdtype = factor(sdtype,levels=c("phylo", "phyloX", "cor", "resid")
                             , labels=c("Phylogenetic", "Phylogenetic Slope", "Correlation", "Tip"))
             , type = ifelse(type=="cor","Correlated Slope","Uncorrelated Slope")
             , Platform = platform
             , size = factor(size,levels=c("small","med","large"),labels=c("Small", "Medium", "Large"))
  )
)

ss_slope_dummy <- data.frame(
  sdtype = factor(c("Phylogenetic", "Phylogenetic Slope", "Correlation", "Tip"
                    , "Phylogenetic", "Phylogenetic Slope", "Correlation", "Tip"))
  , type = factor(c("Correlated Slope", "Correlated Slope", "Correlated Slope", "Correlated Slope"
                    ,"Uncorrelated Slope" , "Uncorrelated Slope","Uncorrelated Slope","Uncorrelated Slope"))
  , Y=c(4,2,0.7,10,4,2,NA,10))

gg_ss_slope <- (ggplot(data=ss_slope_data,aes(x=size,y=sd))
  + facet_grid(type~sdtype, scale="free_y")
  + geom_violin()
  + geom_hline(data=ss_slope_dummy,aes(yintercept = Y))
  + scale_y_log10(limits=c(0.5,20),breaks=c(0.1,1,2,4,10,15))
  + ylab("Standard Deviation")
  + ggtitle("LME4 Single site slope (Need to draw a separate correlation plot")
)

print(gg_ss_slope)

## multiple sites

ms_data <- (pez_data
  %>% select(-phyloX)
  %>% rbind(.,lme4ms_data)
  %>% gather(key=sdtype, value=sd, -c(time,model))
  %>% separate(model,c("platform","size"),"_")
  %>% select(platform, size, sdtype, sd, time)
  %>% mutate(sdtype = ifelse(sdtype=="phylo","Phylogenetic","Tip")
    , Platform = factor(platform)
    , Platform = factor(Platform,levels=c("pez","lme4"))
    , size = factor(size,levels=c("small","med","large"),labels=c("Small", "Medium", "Large"))
  )
)  

ms_dummy <- data.frame(sdtype = factor(c("Phylogenetic","Tip"))
                       , Y=c(4,10))

gg_ms <- (ggplot(data=ms_data,aes(x=size,y=sd,fill=Platform))
  + facet_grid(.~sdtype, scale="free_y")
  + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
  + geom_hline(data=ms_dummy,aes(yintercept=Y))
  + scale_y_continuous(limits=c(0,12),breaks=c(0,2,4,6,8,10,12))
  + ggtitle("Multiple Site using pez and lme4")
  + ylab("Standard Deviation")
)

print(gg_ms)

gg_mstime <- (ggplot(data=ms_data,aes(x=size,y=time,fill=Platform))
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
  %>% mutate(sdtype = factor(sdtype,levels=c("phylo", "phyloX", "cor", "resid")
                             , labels=c("Phylogenetic", "Phylogenetic Slope", "Correlation", "Tip"))
             , type = ifelse(type=="cor","Correlated Slope","Uncorrelated Slope")
             , Platform = factor(platform)
             , Platform = factor(Platform, levels=c("pez","lme4"))
             , size = factor(size,levels=c("small","med","large"),labels=c("Small", "Medium", "Large"))
  )
)

gg_ms_slope <- (ggplot(data=ms_slope_data,aes(x=size,y=sd,fill=Platform))
  + facet_grid(type~sdtype, scale="free_y")
  + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
  # + scale_y_log10()
  # + geom_hline(data=ss_slope_dummy,aes(yintercept = Y))
  + scale_y_log10(limits=c(0.5,20),breaks=c(0.1,1,2,4,10,15))
  + ylab("Standard Deviation")
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


#cs_data <- (pez_cs_data
#  %>% select(-phyloX)
#  %>% mutate(model=ifelse(model=="pez_small","pez_cs_small","pez_cs_med"))
#  %>% rbind(.,lme4cs_data)
#  %>% gather(key=sdtype, value=sd, -c(time,model))
#  %>% separate(model,c("platform","type","size"),"_")
#  %>% select(platform, size, sdtype, sd, time)
#)

#gg_cs <- (ggplot(data=cs_data,aes(x=size,y=sd,fill=platform))
#  + facet_grid(.~sdtype, scale="free_y")
#  + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
#  + ggtitle("Compound symmetric using pez and lme4")
#)

#print(gg_cs)

#gg_cstime <- (ggplot(data=cs_data,aes(x=size,y=time,fill=platform))
#  + geom_violin(position=position_dodge(width=0),alpha=0.4)
#  + scale_y_log10()
#  + ggtitle("CS timing")
#)

#print(gg_cstime)

## cs slope

#lme4cs_slope_data[,"sp_by_site"] <- lme4cs_slope_data[,"sp:site"]
#lme4cs_slope_data <- (lme4cs_slope_data 
#  %>% mutate(model = paste("lme4_",model,sep=""))
#  %>% select(resid,phylo,phyloX,cor,sp_by_site,model,time)
#)
#cs_slope_data <- (pez_cs_slope_data
#  %>% mutate(cor = NA
#        , model=paste("pez_",model,sep="")
#        , temp = phylo        ## hacking, I feel like pez gets the ordering confuse
#        , phylo = phyloX
#        , phyloX = temp
#      )
#  %>% select(-temp)
#  %>% rbind(.,lme4cs_slope_data)
#  %>% gather(key=sdtype, value=sd, -c(time,model))
#  %>% separate(model,c("platform","size","type"),"_")
#  %>% select(platform, size, sdtype, type, sd, time)
#)

#gg_cs_slope <- (ggplot(data=cs_slope_data,aes(x=size,y=sd,fill=platform))
#                + facet_grid(type~sdtype, scale="free_y")
#                + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
#                # + scale_y_log10()
#                + ggtitle("CS slope using pez and lme4")
#)

#print(gg_cs_slope)

#gg_cs_slope_time <- (ggplot(data=cs_slope_data,aes(x=size,y=time,fill=platform))
#                     + facet_grid(type~., scale="free_y")
#                     + geom_violin(position=position_dodge(width=0),alpha=0.4)
#                     + scale_y_log10()
#                     + ggtitle("CS slope timing")
#)

#print(gg_cs_slope_time)


## plots for the ms

library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)
theme_set(theme_bw())
zmargin <- theme(panel.margin=grid::unit(0,"lines"))

data_list <- readRDS("./datadir/collect.RDS")
ssdat_raw <- data_list[[1]]
msdat_raw <- data_list[[2]]

ssdat <- (ssdat_raw
	%>% separate(model,c("platform", "sites", "size", "seed","saveformat"), "[.]")
	%>% select(time, platform, size, seed, resid, phylo_int, phylo_X, phylo_cor, B0, B1)
	%>% gather(key = sdtype, value = sd, -c(platform, size, seed, time, B0, B1))
	%>% mutate(size = factor(size
	      , levels=c("small","med","large","xlarge"), labels=c("25","50","100","1000")
			)
			, platform = factor(platform, levels=c("gls", "phylolm","lme4", "brms"))
		)
	%>% filter(platform != "brms")
)


gg_ss <- (ggplot(data=ssdat, aes(x=size, y=sd, fill=platform))
  + facet_wrap(~sdtype, scale="free_y")
  + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
  + ggtitle("Single site")
)

print(gg_ss)

gg_sstime <- (ggplot(data=ssdat, aes(x=size, y=time, fill=platform))
	+ geom_violin(position=position_dodge(width=0.2),alpha=0.4)
	+ scale_y_log10()
	+ ggtitle("Single site time")
)

print(gg_sstime)

ss_coverage <- (ssdat
	%>% group_by(platform, size)
	%>% summarise(B0_coverage = mean(B0, na.rm=TRUE)
		, B1_coverage = mean(B1, na.rm=TRUE)
		)
	%>% gather(key=fixed_parameter, value=coverage, -c(platform, size))
)

gg_sscoverage <- (ggplot(data=ss_coverage
	, aes(x=size, y=coverage, shape=fixed_parameter, colour=platform)
	)
	+ facet_wrap(~platform, nrow = 1)
	+ geom_point(size=4)
	+ geom_hline(aes(yintercept=0.95))
	+ ggtitle("Single Site coverage")
)

print(gg_sscoverage)



msdat <- (msdat_raw
	%>% separate(model,c("platform", "sites", "size", "seed", "saveformat"),"[.]")
	%>% select(-c(sites,seed,saveformat))
	%>% gather(key=sdtype, value=sd, -c(platform,size,time,B0,B1))
	%>% mutate(size = factor(size,
			levels=c("small","med","large","xlarge"), labels=c("25","50","100","1000")
			)
			, platform = factor(platform, levels=c("pez","phyr","lme4"))
		)
)

gg_ms <- (gg_ss
	%+% msdat
	+ ggtitle("Multiple Sites (pez is problematic)")
)

print(gg_ms)

print(gg_ms
	%+% (msdat %>% filter(platform %in% c("lme4")))
	+ ggtitle("TODO: cowplot y-int simulation parameter")
)

gg_mstime <- (gg_sstime
	%+% msdat
	+ ggtitle("multiple site time")
)

print(gg_mstime)



ms_coverage <- (msdat
  %>% group_by(platform, size)
                %>% summarise(B0_coverage = mean(B0, na.rm=TRUE)
                              , B1_coverage = mean(B1, na.rm=TRUE)
                )
                %>% gather(key=fixed_parameter, value=coverage, -c(platform, size))
)

print(gg_mscoverage <- gg_sscoverage %+% ms_coverage)

quit()
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
  + zmargin
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
  + zmargin
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
  + zmargin
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
  + geom_hline(data=ss_slope_dummy,aes(yintercept = Y))
  + ylab("Standard Deviation")
  + ggtitle("Multiple Site using pez and lme4")
  + zmargin
)

print(gg_ms_slope)

gg_ms_slope_time <- (ggplot(data=ms_slope_data,aes(x=size,y=time,fill=platform))
              + facet_grid(type~., scale="free_y")
              + geom_violin(position=position_dodge(width=0),alpha=0.4)
              + scale_y_log10()
              + ggtitle("Multiple Site slope timing")
              + zmargin
)

print(gg_ms_slope_time)


# multiple sites cs


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
 + zmargin
)

print(gg_cs)

gg_cstime <- (ggplot(data=cs_data,aes(x=size,y=time,fill=platform))
 + geom_violin(position=position_dodge(width=0),alpha=0.4)
 + scale_y_log10()
 + ggtitle("CS timing")
 + zmargin
)

print(gg_cstime)

# cs slope

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
               + zmargin
)

print(gg_cs_slope)

gg_cs_slope_time <- (ggplot(data=cs_slope_data,aes(x=size,y=time,fill=platform))
                    + facet_grid(type~., scale="free_y")
                    + geom_violin(position=position_dodge(width=0),alpha=0.4)
                    + scale_y_log10()
                    + ggtitle("CS slope timing")
                    + zmargin
)

print(gg_cs_slope_time)


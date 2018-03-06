## plots for the ms

library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)
theme_set(theme_bw())
zmargin <- theme(panel.margin=grid::unit(0,"lines"))

data_list <- readRDS("./datadir/result_list.RDS")

lme4ms_slope_data <- data_list[[6]]
pez_slope_data <- data_list[[7]]

ss_slope_dummy <- data.frame(sdtype = factor(c("Phylogenetic Intercept", "Phylogenetic Slope"
    , "Correlation", "Tip", "Phylogenetic Intercept", "Phylogenetic Slope"
    , "Correlation", "Tip")
    )
  , type = factor(c("Correlated Slope", "Correlated Slope"
    , "Correlated Slope", "Correlated Slope" ,"Uncorrelated Slope" 
    , "Uncorrelated Slope","Uncorrelated Slope","Uncorrelated Slope")
    )
  , Y=c(4,2,0.7,10,4,2,NA,10)
)

ms_slope_data <- (pez_slope_data
	%>% mutate(cor = NA
		, model=paste("pez_",model,sep="")
		)
	%>% rbind(.,lme4ms_slope_data)
	%>% gather(key=sdtype, value=sd, -c(time,model))
	%>% separate(model,c("platform","size","type"),"_")
	%>% select(platform, size, sdtype, type, sd, time)
	%>% mutate(sdtype = factor(sdtype,levels=c("phylo", "phyloX", "cor", "resid")
		, labels=c("Phylogenetic Intercept", "Phylogenetic Slope", "Correlation", "Tip"))
		, type = ifelse(type=="cor","Correlated Slope","Uncorrelated Slope")
		, type = ifelse(platform == "pez", "Correlated Slope", type)
		, Platform = factor(platform)
		, Platform = factor(Platform, levels=c("pez","lme4"))
		, size = factor(size,levels=c("small","med","large"),labels=c("Small", "Medium", "Large"))
		)
	%>% filter(type=="Correlated Slope")
)

gg_ms_slope <- (ggplot(data=ms_slope_data,aes(x=size,y=sd,fill=Platform))
	+ facet_wrap(~sdtype, scale="free_y")
	+ geom_violin(position=position_dodge(width=0.2),alpha=0.4)
	# + scale_y_continuous(limits=c(-1,10))
	# + scale_y_log10()
	# + geom_hline(data=ss_slope_dummy,aes(yintercept = Y))
	+ geom_hline(data=ss_slope_dummy,aes(yintercept = Y))
	+ ylab("Standard Deviation")
	# + ggtitle("Multiple Site using pez and lme4")
	+ zmargin
)

print(gg_ms_slope)

gg_ms_slope_time <- (ggplot(data=ms_slope_data,aes(x=size,y=time,fill=platform))
	# + facet_grid(type~., scale="free_y")
	+ geom_violin(position=position_dodge(width=0),alpha=0.4)
	+ scale_y_log10()
	+ ylab("Time (Seconds)")
	# + ggtitle("Multiple Site slope timing")
	+ zmargin
)

print(gg_ms_slope_time)

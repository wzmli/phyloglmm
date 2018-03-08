## plots for the ms

library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)
theme_set(theme_bw())
zmargin <- theme(panel.margin=grid::unit(0,"lines"))

data_list <- readRDS("./datadir/result_list.RDS")

lme4cs_slope_data <- data_list[[10]]
pez_cs_slope_data <- data_list[[11]]

cs_slope_dummy <- data.frame(sdtype = factor(c("Phylogenetic Intercept"
		, "Phylogenetic Slope", "Correlation", "Species by Site Interaction", "Tip", "Phylogenetic Intercept"
		, "Phylogenetic Slope", "Correlation", "Species by Site Interaction", "Tip")
		)
	, type = factor(c("Correlated Slope", "Correlated Slope"
      , "Correlated Slope", "Correlated Slope" , "Correlated Slope" 
      , "Uncorrelated Slope", "Uncorrelated Slope","Uncorrelated Slope"
      , "Uncorrelated Slope", "Uncorrelated Slope")
		)
	, Y=c(4,2,0.7,6,10,4,2,NA,6,10)
)


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
  %>% mutate(sdtype = factor(sdtype,levels=c("phylo", "phyloX", "cor", "sp_by_site", "resid")
       , labels=c("Phylogenetic Intercept", "Phylogenetic Slope", "Correlation", "Species by Site Interaction", "Tip"))
       , type = ifelse(type=="cor","Correlated Slope","Uncorrelated Slope")
       , type = ifelse(platform == "pez", "Correlated Slope", type)
       , Platform = factor(platform)
       , Platform = factor(Platform, levels=c("pez","lme4"))
       , size = factor(size,levels=c("small","med","large"),labels=c("Small", "Medium", "Large"))
)
  %>% filter(type=="Correlated Slope")
)

gg_cs_slope <- (ggplot(data=cs_slope_data,aes(x=size,y=sd,fill=Platform))
+ facet_wrap(~sdtype, scale="free_y")
	+ geom_violin(position=position_dodge(width=0.2),alpha=0.4)
  + geom_hline(data=cs_slope_dummy,aes(yintercept = Y))
  + ylab("Standard Deviation")
	# + scale_y_log10()
	# + ggtitle("CS slope using pez and lme4")
	+ zmargin
)

print(gg_cs_slope)

gg_cs_slope_time <- (ggplot(data=cs_slope_data,aes(x=size,y=time,fill=platform))
	# + facet_grid(type~., scale="free_y")
	+ geom_violin(position=position_dodge(width=0),alpha=0.4)
	+ scale_y_log10()
	+ ylab("Time (Seconds)")
	# + ggtitle("CS slope timing")
	+ zmargin
)

print(gg_cs_slope_time)



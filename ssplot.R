library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)
theme_set(theme_bw())
zmargin <- theme(panel.margin=grid::unit(0,"lines"))

data_list <- readRDS("./datadir/result_list.RDS")


gls_data <- data_list[[1]]
lme4ss_data <- data_list[[2]]
lme4ss_slope_data <- data_list[[3]]

gls_data2 <- (gls_data
	%>% transmute(sd=resid
		, sdtype = "phylo"
		, model = model
		, time = time)
	%>% separate(model,c("platform","size","type"),"_")
	%>% select(platform, type, size, sdtype, sd, time)
	%>% mutate(sdtype = ifelse(type %in% c("nophy","noresidnophy"), "resid", sdtype)
		, type = ifelse(type == "nophy", "noresid",type)
		, type = ifelse(type == "noresidnophy", "resid", type)
	)
)

gg_gls <- (ggplot(data=gls_data2,aes(x=size,y=sd))
	+ facet_grid(type~sdtype, scale="free_y")
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
	
gg_ss	<- (ggplot(data=ss_data,aes(x=size,y=sd,fill=Platform))
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

#print(gg_ss)


## single site lme4 slope plots

ss_slope_data <- (lme4ss_slope_data
	%>% gather(key=sdtype, value=sd, -c(time,model))
	%>% separate(model,c("platform","size","type"),"_")
	%>% select(platform, type, size, sdtype, sd, time)
	%>% mutate(sdtype = factor(sdtype,levels=c("phylo", "phyloX", "cor", "resid")
		, labels=c("Phylogenetic Intercept", "Phylogenetic Slope", "Correlation", "Tip"))
		, type = ifelse(type=="cor","Correlated Slope","Uncorrelated Slope")
		, Platform = platform
		, size = factor(size,levels=c("small","med","large"),labels=c("Small", "Medium", "Large"))
	)
)

ss_slope_dummy <- data.frame(
	sdtype = factor(c("Phylogenetic Intercept", "Phylogenetic Slope"
	, "Correlation", "Tip", "Phylogenetic Intercept", "Phylogenetic Slope"
	, "Correlation", "Tip")
	)
	, type = factor(c("Correlated Slope", "Correlated Slope"
	, "Correlated Slope", "Correlated Slope" ,"Uncorrelated Slope" 
	, "Uncorrelated Slope","Uncorrelated Slope","Uncorrelated Slope")
	)
	, Y=c(4,2,0.7,10,4,2,NA,10)
)

gg_ss_slope <- (ggplot(data=ss_slope_data,aes(x=size,y=sd))
	+ facet_grid(type~sdtype, scale="free_y")
	+ geom_violin()
	+ geom_hline(data=ss_slope_dummy,aes(yintercept = Y))
	+ scale_y_log10(limits=c(0.5,20),breaks=c(0.1,1,2,4,10,15))
	+ ylab("Standard Deviation")
	+ ggtitle("LME4 Single site slope (Need to draw a separate correlation plot")
	+ zmargin
)

#print(gg_ss_slope)

print(head(ss_slope_data))

dat <- (ss_slope_data
	%>% filter(type == "Correlated Slope")                 
)

gg_ss_slope2 <- (ggplot(data=dat, aes(x=size, y=sd))
	+ facet_wrap(~sdtype, scale="free_y")
	+ geom_boxplot()
	+ geom_hline(data=ss_slope_dummy,aes(yintercept = Y),linetype=2)
	+ ylab("Standard Deviation")
	# + ggtitle("LME4 Single site slope (Need to draw a separate correlation plot")
	+ zmargin
)


print(gg_ss_slope2)

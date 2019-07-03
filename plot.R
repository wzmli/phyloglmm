## plots for the ms

library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)
theme_set(theme_bw())
zmargin <- theme(panel.margin=grid::unit(0,"lines"))

data_list <- readRDS("./datadir/collect.RDS")
ssdat_raw <- rbind(data_list[[1]],brmsss_data)
msdat_raw <- data_list[[2]]

tree_seed = 1
source("parameters.R")

sspar_df <- data.frame(
  sdtype = c("resid", "phylo_int", "phylo_X", "phylo_cor")
  , y_int = c(sd.resid, physd.B0, physd.B1, phyrho.B01)
)

ssdat <- (ssdat_raw
	%>% separate(model,c("platform", "sites", "size", "seed","saveformat"), "[.]")
	%>% dplyr:::select(time, platform, size, seed, resid, phylo_int, phylo_X, phylo_cor, B0, B1)
	%>% gather(key = sdtype, value = sd, -c(platform, size, seed, time, B0, B1))
	%>% left_join(.,sspar_df)
	%>% mutate(size = factor(size
	      , levels=c("small","med","large","xlarge"), labels=c("25","50","100","500")
			)
			, Platform = factor(platform, levels=c("gls", "phylolm","lme4", "glmmTMB", "brms"))
			, sdtype = factor(sdtype, levels=c("phylo_int","phylo_cor","phylo_X","resid")
			                  , labels=c(expression(paste("Phylogenetic random intercept ", Sigma[phy[int]]))
			                             , expression(paste("Phylogenetic random intercept-slope correlation ", Sigma[phy[int-slope]]))
			                             , expression(paste("Phylogenetic random slope ", Sigma[phy[slope]]))
			                             , expression(paste("Residual ", sigma[epsilon])))
		)
	)
	# %>% filter(platform != "brms")
)


gg_ss <- (ggplot(data=ssdat, aes(x=size, y=sd, col=Platform))
  + facet_wrap(~sdtype, scale="free_y", labeller = label_parsed)
  # + geom_violin(position=position_dodge(width=0.5),alpha=0.4)
  + geom_boxplot(outlier.colour = NULL, varwidth=TRUE)
  + geom_hline(aes(yintercept = y_int))
  # + scale_fill_manual(values=c("Dark Red","Dark Blue","Dark Green","Purple"))
  # + scale_fill_brewer(palette = "Dark2")
  # + scale_color_brewer(palette = "Dark2")
  + scale_color_manual(values=c("Black","Red","Dark Blue","Dark Green","Orange"))
  # + ggtitle("Single site")
  + theme_bw()
  + ylab("Parameter Estimates")
  + xlab("Number of Species")
)

print(gg_ss)

gg_sstime <- (ggplot(data=ssdat, aes(x=size, y=time, col=Platform))
	# + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
	+ geom_boxplot(outlier.colour = NULL, varwidth=TRUE)
	# + scale_fill_brewer(palette = "Dark2")
	+ scale_y_log10()
	+ scale_color_manual(values=c("Black","Red","Dark Blue","Dark Green","Orange"))
	
	+ ylab("Time (seconds)")
	# + ggtitle("Single site time")
  + theme_bw()
)

print(gg_sstime)

ss_coverage <- (ssdat
	%>% group_by(Platform, size)
	%>% summarise(B0_coverage = mean(B0, na.rm=TRUE)
		, B1_coverage = mean(B1, na.rm=TRUE)
		)
	%>% gather(key=fixed_parameter, value=coverage, -c(Platform, size))
	%>% mutate(Parameter = factor(fixed_parameter, labels=c(expression(beta[0])
	                                                        , expression(beta[1])
	                                                        ))
	)
)

gg_sscoverage <- (ggplot(data=ss_coverage
	, aes(x=size, y=coverage, shape=Parameter, colour=Platform)
	)
	+ facet_wrap(~Platform, nrow = 1)
	+ geom_point(size=4)
	+ scale_shape_discrete("Parameters",labels=c(expression(beta[0])
	                              , expression(beta[1])))
	+ geom_hline(aes(yintercept=0.95))
	# + ggtitle("Single Site coverage")
	+ scale_color_manual(values=c("Black","Red","Dark Blue","Dark Green","Orange"))
	+ annotate("rect", xmin=0, xmax=4
	           , ymin=0.95 - 2*sqrt(0.95*0.05/100)
	           , ymax=0.95 + 2*sqrt(0.95*0.05/100), alpha=0.2)
  + theme(panel.spacing = unit(0,'lines'))
	+ xlab("Number of Species")
	+ ylab("Coverage")
)

print(gg_sscoverage)

mspar_df <- data.frame(
  sdtype = c("resid", "phylo_X", "phylo_int", "phylo_cor", "phylo_interaction"
             , "species_X", "species_int", "species_cor", "site_int")
  , y_int = c(sd.resid, physd.B1, physd.B0, phyrho.B01, sd.interaction
              , sd.B1, sd.B0, rho.B01, sd.site)
)

msdat <- (msdat_raw
	%>% separate(model,c("platform", "sites", "size", "seed", "saveformat"),"[.]")
	%>% dplyr:::select(-c(sites,seed,saveformat))
	%>% filter((convergence != 1) | is.na(convergence))
	%>% dplyr:::select(-convergence)
	%>% gather(key=sdtype, value=sd, -c(platform,size,time,B0,B1))
	%>% filter(sd <20)
	%>% left_join(.,mspar_df)
	%>% mutate(size = factor(size,
			levels=c("small","med","large","xlarge"), labels=c("25","50","100","500")
			)
			, Platform = factor(platform, levels=c("pez","phyr","lme4", "glmmTMB"))
			, sdtype = factor(sdtype, levels=c("phylo_int","phylo_cor","phylo_X"
			                                   , "species_int", "species_cor", "species_X"
			                                   , "phylo_interaction", "site_int", "resid")
			                  , labels=c(expression(paste(Sigma[phy[int]]))
			                             , expression(paste(Sigma[phy[int-slope]]))
			                             , expression(paste(Sigma[phy[slope]]))
			                             , expression(paste(sigma[sp[int]]))
			                             , expression(paste(sigma[sp[int-slope]]))
			                             , expression(paste(sigma[sp[slope]]))
			                             , expression(paste(Sigma[phy[sp:site]]))
			                             , expression(paste(sigma[site]))
			                             , expression(sigma[epsilon]))
			)
		)
)

gg_ms <- (gg_ss
	%+% msdat
	+ scale_color_manual(values=c("Gray","Purple","Dark Blue","Orange"))

)

print(gg_ms)

gg_mstime <- (gg_sstime
	%+% msdat
	# + ggtitle("multiple site time")
	+ scale_color_manual(values=c("Gray","Purple","Dark Blue","Orange"))
	
)

print(gg_mstime)

ms_coverage <- (msdat
                %>% filter(!(sd %in% c(-1,1)))
                %>% group_by(Platform, size)
                %>% summarise(B0_coverage = mean(B0, na.rm=TRUE)
                              , B1_coverage = mean(B1, na.rm=TRUE)
                )
                %>% gather(key=fixed_parameter, value=coverage, -c(Platform, size))
                %>% mutate(Parameter = factor(fixed_parameter, labels=c(expression(beta[0])
                                                                        , expression(beta[1])
                ))
                )
)

gg_mscoverage <- (ggplot(data=ms_coverage
                         , aes(x=size, y=coverage, shape=Parameter, colour=Platform)
)
+ facet_wrap(~Platform, nrow = 1)
+ geom_point(size=4)
+ geom_hline(aes(yintercept=0.95))
+ scale_shape_discrete("Parameters",labels=c(expression(beta[0])
                                             , expression(beta[1])))
# + ggtitle("Single Site coverage")
+ annotate("rect", xmin=0, xmax=5
           , ymin=0.95 - 2*sqrt(0.95*0.05/100)
           , ymax=0.95 + 2*sqrt(0.95*0.05/100), alpha=0.2)
+ theme(panel.spacing = unit(0,'lines'))
)

print(mscoverage<- gg_mscoverage 	+ scale_color_manual(values=c("Gray","Purple","Dark Blue","Orange"))
)

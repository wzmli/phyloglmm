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
  + geom_hline(aes(yintercept = yint))
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


ms_coverage <- (msdat
                %>% filter(!(sd %in% c(-1,1)))
                %>% group_by(platform, size)
                %>% summarise(B0_coverage = mean(B0, na.rm=TRUE)
                              , B1_coverage = mean(B1, na.rm=TRUE)
                )
                %>% gather(key=fixed_parameter, value=coverage, -c(platform, size))
)

print(gg_mscoverage <- gg_sscoverage %+% ms_coverage)


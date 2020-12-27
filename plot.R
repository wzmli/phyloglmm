## plots for the ms

library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)
library(scales)
theme_set(theme_bw())
zmargin <- theme(panel.spacing=grid::unit(0,"lines"))

data_list <- readRDS("datadir/collect_rerun.RDS")
ssdat_raw <- data_list[[1]]
msdat_raw <- data_list[[2]]

tree_seed = 1
source("parameters.R")

sspar_df <- data.frame(
  sdtype = c("resid", "phylo_int", "phylo_X", "phylo_cor")
  , y_int = c(sd.resid^2, physd.B0^2, physd.B1^2, phyrho.B01*physd.B0*physd.B1)
)

ssdat <- (ssdat_raw
	%>% separate(model,c("platform", "sites", "size", "seed","saveformat"), "[.]")
	%>% dplyr:::select(time, platform, size, seed, resid, phylo_int, phylo_X, phylo_cor, B0, B1)
	%>% gather(key = sdtype, value = sd, -c(platform, size, seed, time, B0, B1))
	%>% left_join(.,sspar_df)
	%>% mutate(size = factor(size
	      , levels=c("small","med","large","xlarge"), labels=c("25","50","100","500")
			)
			, Platform = factor(platform, levels=c("gls", "phylolm","lme4", "glmmTMB", "brms","MCMCglmm"))
			, sdtype = factor(sdtype, levels=c("phylo_int","phylo_cor","phylo_X","resid")
			                  , labels=c(expression(paste("Phylogenetic random intercept ", Sigma[phy[int]]))
			                             , expression(paste("Phylogenetic random intercept-slope correlation ", rho[phy[int-slope]]))
			                             , expression(paste("Phylogenetic random slope ", Sigma[phy[slope]]))
			                             , expression(paste("Residual ", sigma[epsilon])))
		)
	)
	# %>% filter(platform != "brms")
)

## FIXME: could match these more programatically (use 'values' in
## scale_colour_manual() ?)
colvec_all <-  c(gls="Black",phylolm="Red",lme4="Dark Blue",
                 glmmTMB="Dark Green",brms="Orange",pez="Gray",phyr="Purple",MCMCglmm="Yellow")
colvec  <- colvec_all[c("gls","phylolm","lme4","glmmTMB","brms","MCMCglmm")]
colvec2 <- colvec_all[c("pez","phyr","lme4","glmmTMB")]

gg_ss <- (ggplot(data=ssdat, aes(x=size, y=sd, col=Platform, fill=Platform))
  + scale_y_continuous(limits = c(0,200), oob=scales::squish)
  + facet_wrap(~sdtype, scale="free_y", labeller = label_parsed)
  # + geom_violin(position=position_dodge(width=0.5),alpha=0.4)
  + geom_boxplot(outlier.colour = NULL, varwidth=TRUE, alpha=0.2)
  + geom_hline(aes(yintercept = y_int), lty=2)
  + scale_fill_manual(values=colvec)
  # + scale_fill_brewer(palette = "Dark2")
  # + scale_color_brewer(palette = "Dark2")
  + scale_color_manual(values=colvec)
  # + ggtitle("Single site")
  + theme_bw()
  + ylab("Parameter Estimates")
  + xlab("Number of Species")
)

print(gg_ss)
ggsave(plot = gg_ss,filename = "figure/ssplot.pdf",width = 10, height=7)

scaleFUN <- function(x) sprintf("%.2f", x)

gg_sstime <- (ggplot(data=ssdat, aes(x=size, y=time, col=Platform,
                                     fill=Platform))
	# + geom_violin(position=position_dodge(width=0.2),alpha=0.4)
	+ geom_boxplot(outlier.colour = NULL, varwidth=TRUE, alpha=0.2)
	# + scale_fill_brewer(palette = "Dark2")
	# + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)
	#     , labels = trans_format("log10", math_format(10^.x))
	#   )
	+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = scaleFUN)
	+ scale_color_manual(values=colvec)
	+ scale_fill_manual(values=colvec)
    + ylab("Time (seconds)")
    + xlab("Number of Species")
	# + ggtitle("Single site time")
  + theme_bw()
)

print(gg_sstime)
ggsave(plot = gg_sstime,filename = "figure/sstime.pdf",width = 10, height=7)


b_ci <- function(x,w) {
    x <- na.omit(x)
    binom.test(sum(x),length(x))$conf.int[w]
}

ss_coverage <- (ssdat
    %>% as_tibble()
    %>% tidyr::gather(key=fixed_parameter, value=cov, B0,B1)
    %>% group_by(Platform, size, fixed_parameter)
    %>% summarise(coverage=mean(cov, na.rm=TRUE),
                  lwr=b_ci(cov,1),
                  upr=b_ci(cov,2))
    %>% mutate(Parameter = factor(fixed_parameter,
                                  labels=c(expression(beta[0])
                                         , expression(beta[1])
                                           )))
)
                  
## ss_coverage <- (ssdat
## 	%>% group_by(Platform, size)
## 	%>% summarise(B0_coverage = mean(B0, na.rm=TRUE)
## 		, B1_coverage = mean(B1, na.rm=TRUE)
## 		)
## 	%>% gather(key=fixed_parameter, value=coverage, -c(Platform, size))
## 	%>% mutate(Parameter = factor(fixed_parameter, labels=c(expression(beta[0])
## 	                                                        , expression(beta[1])
## 	                                                        ))
## 	)
## )

gg_sscoverage <- (ggplot(data=ss_coverage
	, aes(x=size, y=coverage, shape=Parameter, colour=Platform)
	)
	+ facet_wrap(~Platform, nrow = 1)
        + geom_point(size=4, alpha=0.5)
        # + geom_linerange(aes(ymin=lwr,ymax=upr))
	+ scale_shape_discrete("Parameters",labels=c(expression(beta[0])
	                              , expression(beta[1])))
	+ geom_hline(aes(yintercept=0.95))
	# + ggtitle("Single Site coverage")
	+ scale_color_manual(values=colvec)
	+ annotate("rect", xmin=0, xmax=4
	           , ymin=0.95 - 2*sqrt(0.95*0.05/100)
	           , ymax=0.95 + 2*sqrt(0.95*0.05/100), alpha=0.1
	           , fill="blue")
  + theme(panel.spacing = unit(0,'lines'))
	+ xlab("Number of Species")
	+ ylab("Coverage")
)

print(gg_sscoverage)
ggsave(plot = gg_sscoverage,filename = "figure/sscoverage.pdf",width = 10, height=5)

mspar_df <- data.frame(
  sdtype = c("resid", "phylo_X", "phylo_int", "phylo_cor", "phylo_interaction"
             , "species_X", "species_int", "species_cor", "site_int")
  , y_int = c(sd.resid^2, physd.B1^2, physd.B0^2, phyrho.B01*physd.B0*physd.B1, sd.interaction^2
              , sd.B1^2, sd.B0^2, rho.B01*sd.B0*sd.B1, sd.site^2)
)

msdat <- (msdat_raw
	%>% separate(model,c("platform", "sites", "size", "seed", "saveformat"),"[.]")
	%>% dplyr:::select(-c(sites,seed,saveformat))
	%>% filter((convergence != 1) | is.na(convergence))
	%>% dplyr:::select(-convergence)
	%>% gather(key=sdtype, value=sd, -c(platform,size,time,B0,B1))
	# %>% filter(sd <20)
	%>% left_join(.,mspar_df)
	%>% mutate(size = factor(size,
			levels=c("small","med","large","xlarge"), labels=c("25","50","100","500")
			)
			, Platform = factor(platform, levels=c("pez","phyr","lme4", "glmmTMB"))
			, sdtype = factor(sdtype, levels=c("phylo_int","phylo_cor","phylo_X"
			                                   , "species_int", "species_cor", "species_X"
			                                   , "phylo_interaction", "site_int", "resid")
			                  , labels=c(expression(paste("Phylogenetic random intercept ",Sigma[phy[int]]))
			                             , expression(paste("Phylo random intercept-slope correlation ",rho[phy[int-slope]]))
			                             , expression(paste("Phylogenetic random slope ",Sigma[phy[slope]]))
			                             , expression(paste("Species random intercept ",sigma[sp[int]]))
			                             , expression(paste("Species random intercept-slope correlation ",rho[sp[int-slope]]))
			                             , expression(paste("Species random slope ",sigma[sp[slope]]))
			                             , expression(paste("Phylo random species-group interaction ",Sigma[phy[sp:group]]))
			                             , expression(paste("Random group ",sigma[group]))
			                             , expression(paste("Residual ",sigma[epsilon])))
			)
		)
)


gg_ms <- (gg_ss
	%+% msdat
	+ scale_color_manual(values=colvec2)
	+ scale_fill_manual(values=colvec2)
	+ theme(legend.position = "bottom")

)

print(gg_ms)
ggsave(plot = gg_ms,filename = "figure/msplot.pdf",width = 10, height=7)



gg_mstime <- (gg_sstime
	%+% msdat
	# + ggtitle("multiple site time")
	+ scale_color_manual(values=colvec2)
	
)

print(gg_mstime)
ggsave(plot = gg_mstime,filename = "figure/mstime.pdf",width = 10, height=7)


## BMB: see revised code above
ms_coverage <- (msdat
                # %>% filter(!(sd %in% c(-1,1)))
                %>% as_tibble()
                %>% tidyr::gather(key=fixed_parameter, value=cov, B0,B1)
                %>% group_by(Platform, size, fixed_parameter)
                %>% summarise(coverage=mean(cov, na.rm=TRUE),
                              lwr=b_ci(cov,1),
                              upr=b_ci(cov,2))
                %>% mutate(Parameter = factor(fixed_parameter, labels=c(expression(beta[0])
                                                                        , expression(beta[1])
                ))
                )
)

gg_mscoverage <- (ggplot(data=ms_coverage
                         , aes(x=size, y=coverage, shape=Parameter, colour=Platform)
)
+ facet_wrap(~Platform, nrow = 1)
+ geom_point(size=4, alpha=0.5)
# + geom_linerange(aes(ymin=lwr,ymax=upr))
+ geom_hline(aes(yintercept=0.95))
+ scale_shape_discrete("Parameters",labels=c(expression(beta[0])
                                             , expression(beta[1])))
# + ggtitle("Single Site coverage")
+ annotate("rect", xmin=0, xmax=5
           , ymin=0.95 - 2*sqrt(0.95*0.05/100)
           , ymax=0.95 + 2*sqrt(0.95*0.05/100), alpha=0.1, fill="blue")
    + theme(panel.spacing = unit(0,'lines'))
    + xlab("Number of Species")
)

print(mscoverage<- gg_mscoverage 	+ scale_color_manual(values=colvec2))
ggsave(plot = mscoverage,filename = "figure/mscoverage.pdf",width = 10, height=5)

pp <- readRDS("datadir/lme4_ms_small_profile.RDS")
pp2 <- pp %>% mutate(b0 = between(0,B0_lower,B0_upper), b1=between(0,B1_lower,B1_upper))
profile_dat <- data.frame(Platform = "lme4", size=factor(25)
                          ,fixed_parameter=c("B0","B1") ,Parameter=c("beta[0]","beta[1]")
                          , coverage = c(sum(pp2$b0,na.rm=TRUE)/nrow(pp2), sum(pp2$b1,na.rm=TRUE)/nrow(pp2)))
print(gg_mscoverage 	+ scale_color_manual(values=colvec2)
      + geom_point(data=profile_dat,color="red",size=3))


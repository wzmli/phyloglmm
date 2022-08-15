## plots for the ms
library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)
library(scales)
library(colorspace)
library(directlabels)
library(ggrepel)

plot_width <- 7

theme_set(theme_bw())
zmargin <- theme(panel.spacing=grid::unit(0,"lines"))

all_pkgs <- c("gls","phylolm","lme4","glmmTMB","brms","pez","phyr","MCMCglmm")
sub_pkgs <- setdiff(all_pkgs, c("pez","phyr"))
sub_pkgs2 <- setdiff(all_pkgs, c("gls","phylolm","brms","MCMCglmm"))
colvec_all <- setNames(colorspace::qualitative_hcl(n=length(all_pkgs)),
                       all_pkgs)
colvec  <- colvec_all[sub_pkgs]
colvec2  <- colvec_all[sub_pkgs2]
##    c(gls="Black",phylolm="Red",lme4="Dark Blue",
##                 glmmTMB="Dark Green",brms="Orange",pez="Gray",phyr="Purple",MCMCglmm="Yellow")

msdat_raw <- readRDS("datadir/ms_lme4verse.RDS")

tree_seed = 1
source("parameters.R")

sspar_df <- data.frame(
  sdtype = c("resid", "phylo_int", "phylo_X", "phylo_cor")
  , y_int = c(sd.resid^2, physd.B0^2, physd.B1^2, phyrho.B01*physd.B0*physd.B1)
)

trans_platform <- function(platform) {
  ifelse(platform %in% c("lme4","glmmTMB"),
         sprintf("phyloglmm/\n%s", platform),
         platform)
}



colvec_trans <- setNames(colvec, trans_platform(names(colvec)))

mspar_df <- data.frame(
  sdtype = c("resid", "phylo_X", "phylo_int", "phylo_cor", "phylo_interaction"
             , "species_X", "species_int", "species_cor", "site_int")
  , y_int = c(sd.resid^2, physd.B1^2, physd.B0^2, phyrho.B01*physd.B0*physd.B1, sd.interaction^2
              , sd.B1^2, sd.B0^2, rho.B01*sd.B0*sd.B1, sd.site^2)
)


lev_vec <- c("phylo_int","phylo_cor","phylo_X"
           , "species_int", "species_cor", "species_X"
           , "phylo_interaction", "site_int", "resid")

lab_vec <- c(expression(atop("Phylogenetic random intercept ",sigma[phy[int]]^2))
           , expression(atop("Phylo random intercept-slope cov",sigma[phy[int-slope]]))
           , expression(atop("Phylogenetic random slope ",sigma[phy[slope]]^2))
           , expression(atop("Species random intercept ",sigma[sp[int]]^2))
           , expression(atop("Species random intercept-slope cov",sigma[sp[int-slope]]))
           , expression(atop("Species random slope ",sigma[sp[slope]]^2))
           , expression(atop("Phylo random species-group interact",sigma[phy[sp:group]]^2))
           , expression(atop("Random group ",sigma[group]^2))
           , expression(atop("Residual ",sigma[epsilon]^2)))

msdat <- (msdat_raw
    %>% separate(model,c("platform", "sites", "size", "seed", "saveformat"),"[.]")
    %>% dplyr:::select(-c(seed,saveformat))
    %>% filter((convergence != 1) | is.na(convergence))
    %>% dplyr:::select(-c(convergence, rawtime))
    %>% gather(key=sdtype, value=sd, -c(platform,size,time,B0,B1,sites))
    ## %>% filter(sd <20)
    %>% left_join(.,mspar_df)
    %>% mutate(size = factor(size,
                             levels=c("small","med","large","xlarge"),
                             labels=c("25","50","100","500")
                             )
             , Platform = factor(platform,
                                 levels=c("lme4", "glmmTMB"),
                                 labels=trans_platform(c("lme4", "glmmTMB")),
                                 )
             , sdtype = factor(sdtype, levels = lev_vec,
                             , labels = lab_vec
                               )
               )
    %>% mutate(Platform2 = paste0(Platform,"_",sites))
    )



trans_colvec2 <- setNames(colvec2, trans_platform(names(colvec2)))[levels(msdat$Platform)]
gg_ms <- (ggplot(data=msdat, aes(x=size, y=sd, col=Platform2, fill=Platform2))
		  + scale_y_continuous(limits = c(0,200), oob=scales::squish)
		  + facet_wrap(~sdtype, scale="free_y", labeller = label_parsed)
		  # + geom_violin(position=position_dodge(width=0.5),alpha=0.4)
		  + geom_boxplot(outlier.colour = NULL, varwidth=TRUE, alpha=0.5)
		  + geom_hline(aes(yintercept = y_int), lty=2)
		  # + ggtitle("Single site")
		  + theme_bw()
		  + ylab("Parameter Estimates")
		  + xlab("Number of Species")
	+ theme(legend.position = "bottom")
    + theme(strip.text.x = element_text(size = 12))
)

print(gg_ms)
ggsave("lme4verse_msplot.pdf",plot=gg_ms)

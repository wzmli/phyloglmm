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

gls_path <- "./datadir/gls/"
gls_res <- list.files(path = gls_path, pattern = "full")
gls_results <- function(tt){
  gls_list <- list()
  gls_df <- data.frame(Phylogenetic= numeric(200)
                       , model = numeric(200)
                       , time = numeric(200)
  )
  for(i in tt){
    gls_list[[i]] <- gls_df
    ff <- list.files(path=paste(gls_path,i,"/",sep=""),patter="rds")
    for(j in 1:length(ff)){
      gls_obj <- readRDS(paste(gls_path,i,"/",ff[j],sep=""))
      gls_list[[i]][j,"sd"] <- attr(resid(gls_obj[[1]]),"std")[[1]]
      gls_list[[i]][j,"model"] <- i
      gls_list[[i]][j,"time"] <- gls_obj[[2]][[1]]
    }
  }
  return(gls_list)
}

gls_full_data <- (bind_rows(gls_results(gls_res))
  %>% select(-c(Phylogenetic))
  %>% separate(model,c("platform","size","type"),"_")
  %>% mutate(type = "Correlated Slope"
  , sdtype = "Phylogenetic Signal"
  , size = factor(size,levels=c("small","med","large"),labels=c("Small", "Medium", "Large"))
  )
)



phylolm_path <- "./datadir/phylolm/"
phylolm_res <- list.files(path = phylolm_path)
phylolm_results <- function(pp,tt){
  phylolm_list <- list()
  phylolm_df <- data.frame(resid = numeric(200)
                          , phylo = numeric(200)
                          , model = numeric(200)
                          , time = numeric(200)
  )
  for(i in tt){
    phylolm_list[[i]] <- phylolm_df
    ff <- list.files(path=paste(pp,i,"/",sep=""),pattern="rds")
    for(j in 1:length(ff)){
      phylolm_obj <- readRDS(paste(pp,i,"/",ff[j],sep=""))
      phylolm_list[[i]][j,"resid"] <- sqrt(phylolm_obj[[1]]$sigma2_error)
      phylolm_list[[i]][j,"phylo"] <- sqrt(phylolm_obj[[1]]$sigma2)
      phylolm_list[[i]][j,"model"] <- i
      phylolm_list[[i]][j,"time"] <- phylolm_obj[[2]][[1]]
    }
  }
  return(phylolm_list)
}

phylolm_full_data <- (bind_rows(phylolm_results(phylolm_path,phylolm_res))
  %>% gather(key=sdtype, value=sd, -c(time,model))
  %>% separate(model,c("platform","size"),"_")
  %>% mutate(type = "Correlated Slope"
      , sdtype = factor(sdtype,levels=c("phylo", "phyloX", "cor", "resid")
                  , labels=c("Phylogenetic Intercept", "Phylogenetic Slope", "Correlation", "Tip"))
      , size = factor(size,levels=c("small","med","large"),labels=c("Small", "Medium", "Large"))
      )
)

## single site lme4 slope plots

ss_slope_data <- (lme4ss_slope_data
	%>% gather(key=sdtype, value=sd, -c(time,model))
	%>% separate(model,c("platform","size","type"),"_")
	%>% select(platform, type, size, sdtype, sd, time)
	%>% mutate(sdtype = factor(sdtype,levels=c("phylo", "phyloX", "cor", "resid")
		, labels=c("Phylogenetic Intercept", "Phylogenetic Slope", "Correlation", "Tip"))
		, type = ifelse(type=="cor","Correlated Slope","Uncorrelated Slope")
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



dat <- (ss_slope_data
	%>% filter(type == "Correlated Slope")      
	%>% bind_rows(phylolm_full_data,gls_full_data)
)

gg_ss_slope2 <- (ggplot(data=dat, aes(x=size, y=sd, fill=platform))
	+ facet_wrap(~sdtype, scale="free_y")
	+ geom_violin(position=position_dodge(width=0.5),alpha=0.4)
	+ geom_hline(data=ss_slope_dummy,aes(yintercept = Y),linetype=2)
	+ ylab("Standard Deviation")
	# + ggtitle("LME4 Single site slope (Need to draw a separate correlation plot")
	+ zmargin
)


print(gg_ss_slope2)

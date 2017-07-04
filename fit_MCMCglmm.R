library(MCMCglmm)
library(dplyr)

dat <- (dat
  %>% rowwise()
  %>% mutate(phylo=paste("t",sp,sep=""))
)

nitt <- 5e3 ## was 5e6

inv.phylo <- inverseA(phy,nodes="TIPS",scale=TRUE)

prior <- list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
MCMC_time <- system.time(
	MCMCglmm_fit <- MCMCglmm(Y~X
		, random=~phylo
		, family="gaussian"
		, ginverse=list(phylo=aaa)
		, prior=prior
		, data=dat
		, nitt=nitt
		, burnin=1000
		, thin=nitt/1000
		, verbose=TRUE)
)

print(summary(MCMCglmm_fit))

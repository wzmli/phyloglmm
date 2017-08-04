library(MCMCglmm)
library(Matrix)
library(dplyr)

dat <- (dat
  %>% rowwise()
  %>% mutate(phylo=paste("t",sp,sep="")
      , obs=phylo
      )
)

dat <- data.frame(dat)

inv.phylo <- inverseA(phy,nodes="TIPS",scale=TRUE)

prior <- list(G=list(G1=list(V=2,nu=1)
                     )
              , R=list(V=1,nu=1)
              )
MCMC_time <- system.time(
	MCMCglmm_fit <- MCMCglmm(Y~1
		, random=~ phylo
		, family="gaussian"
		, ginverse=list(phylo=inv.phylo$Ainv)
		, prior=prior
		, data=dat
		, nitt=nitt
		, burnin=1000
		, thin=nitt/1000
		, verbose=TRUE)
)

print(summary(MCMCglmm_fit))

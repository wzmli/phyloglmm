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

prior <- list(G=list(G1=list(V=1,nu=0.02)
#                      , G2=list(V=1,nu=0.02)
#                      , G3=list(V=1,nu=0.02)
                     , G4=list(V=1,nu=0.02)
                     )
              , R=list(V=1,nu=0.02)
              )
MCMC_time <- system.time(
	MCMCglmm_fit <- MCMCglmm(Y~X
		, random=~ phylo + us(X):phylo 
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

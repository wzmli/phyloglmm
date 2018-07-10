library(MCMCglmm)
library(Matrix)
library(dplyr)

dat <- (dat
    %>% rowwise()
	 %>% mutate(phylo=paste("t",sp,sep="")
      , obs=phylo
      )
#	 %>% filter(site == 1) 
)

dat <- data.frame(dat)

inv.phylo <- inverseA(phy,nodes="TIPS",scale=TRUE)

prior <- list(G=list(G1=list(V=10,nu=0.1)
                     )
              , R=list(V=1,nu=1)
              )


MCMC_time <- system.time(
	MCMCglmm_fit <- MCMCglmmhacked(Y~1+X
		, random=~ sp
		, family="gaussian"
		, ginverse=list(sp=inv.phylo$Ainv)
		, prior=prior
		, data=dat
		, nitt=nitt
		, burnin=500
		, thin=nitt/1000
		, singular.ok = TRUE
		, verbose=TRUE)
)

print(summary(MCMCglmm_fit))
print(plot(MCMCglmm_fit))

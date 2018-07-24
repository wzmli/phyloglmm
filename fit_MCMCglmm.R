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

time1 <- proc.time()
inv.phylo <- inverseA(phy,nodes="TIPS",scale=TRUE)

prior <- list(G=list(G1=list(V=diag(2),nu=2)
						, G2 = list(V=1,nu=2)
#						, G3 = list(V=diag(2), nu=2)
                     )
              , R=list(V=1,nu=1)
              )


	MCMCglmm_fit <- MCMCglmm(Y~1+X
		, random=~ us(1+X):sp  + site #+ idh(site):sp
		, family="gaussian"
		, ginverse=list(sp=inv.phylo$Ainv)
		, prior=prior
		, data=dat
		, nitt=nitt
		, burnin=500
		, thin=nitt/10000
		, singular.ok = TRUE
		, verbose=TRUE)

time2 <- proc.time()

tt <- time2-time1
print(summary(MCMCglmm_fit))
print(plot(MCMCglmm_fit))


MCMCglmm_list <- list(MCMCglmm_fit,tt)
saveRDS(MCMCglmm_list,file=paste("datadir/MCMCglmm",size,seed,"rds",sep="."))

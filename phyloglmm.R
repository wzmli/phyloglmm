#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)


t1 <- proc.time()

phyZ <- phylo.to.Z(phy,stand=FALSE)
phyZ <- phyZ[order(rownames(phyZ)),]

dat <- (dat
	%>% mutate(obs = sp)
	%>% ungroup()
	# %>% arrange(sp)
)	

# 
# lme4fit_phy <- phylo_lmm(y_phy ~ X
#                       # + (1 | sp)
#                      	# + (1 | obs)
#                      + (1 + X | sp)
#                      # + (1 + X | obs)
#                      # + (1 | site)
#                      # + (1 | sp:site)
#                      # + (1 | site:sp)
#                      , data=dat
#                      , phylonm = c("sp","sp:site")
#                      , phylo = phy
#                      , phyloZ=phyZ
#                      , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#                      , REML = FALSE
# )
# 
# print(summary(lme4fit_phy))
# 
# lme4fit_nophy <- phylo_lmm(y_nophy ~ X
#                      #	 	+ (1 | sp)
#                      #		+ (1 | obs)
#                      # + (1 + X | sp)
#                      + (1 + X | obs)
#                      # + (1 | site)
#                      # + (1 | sp:site)
#                      # + (1 | site:sp)
#                      , data=dat
#                      , phylonm = c("sp","sp:site")
#                      , phylo = phy
#                      , phyloZ=phyZ
#                      , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#                      , REML = FALSE
# )
# 
# print(summary(lme4fit_nophy))

lme4fit_noint <- phylo_lmm(y_noint ~ X
                           #	 	+ (1 | sp)
                           #		+ (1 | obs)
                           + (1 + X | sp)
                           + (1 + X | obs)
                           # + (1 | site)
                           # + (1 | sp:site)
                           # + (1 | site:sp)
                           , data=dat
                           , phylonm = c("sp","sp:site")
                           , phylo = phy
                           , phyloZ=phyZ
                           , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                           , REML = FALSE
)

print(summary(lme4fit_noint))


lme4fit_intnophy <- phylo_lmm(y_intnophy ~ X
                           #	 	+ (1 | sp)
                           #		+ (1 | obs)
                           # + (1 + X | sp)
                           # + (1 + X | obs)
                           # + (1 | site)
                           # + (1 | sp:site)
                           # + (1 | obs:site)
                           + (1 | site:obs)
                           # + (1 | site:sp)
                           , data=dat
                           , phylonm = c("sp","sp:site")
                           , phylo = phy
                           , phyloZ=phyZ
                           , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                           , REML = FALSE
)

print(summary(lme4fit_intnophy))


lme4fit_intphy <- phylo_lmm(y_intphy ~ X
                         #	 	+ (1 | sp)
                         #		+ (1 | obs)
                         # + (1 + X | sp)
                         # + (1 + X | obs)
                         # + (1 | site)
                         + (1 | sp:site)
                         # + (1 | obs:site)
                         # + (1 | site:sp)
                         , data=dat
                         , phylonm = c("sp","sp:site")
                         , phylo = phy
                         , phyloZ=phyZ
                         , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                         , REML = FALSE
)

print(summary(lme4fit_intphy))


lme4fit_int <- phylo_lmm(y_int ~ X
                            #	 	+ (1 | sp)
                            #		+ (1 | obs)
                            # + (1 + X | sp)
                            # + (1 + X | obs)
                            # + (1 | site)
                            + (1 | sp:site)
                            + (1 | obs:site)
                            # + (1 | site:sp)
                            , data=dat
                            , phylonm = c("sp","sp:site")
                            , phylo = phy
                            , phyloZ=phyZ
                            , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                            , REML = FALSE
)

print(summary(lme4fit_int))
# 
# 
# lme4fit <- phylo_lmm(y ~ X
#                            + (1 + X | sp)
#                            + (1 + X | obs)
#                            + (1 | site)
#                            + (1 | sp:site)
#                            # + (1 | site:sp)
#                            , data=dat
#                            , phylonm = c("sp","sp:site")
#                            , phylo = phy
#                            , phyloZ=phyZ
#                            , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#                            , REML = FALSE
# )
# 
# print(summary(lme4fit))
# 
# 
# 
# lme4fit2 <- phylo_lmm(y ~ X
#                      + (1 + X | sp)
#                      + (1 + X | obs)
#                      + (1 | site)
#                      # + (1 | sp:site)
#                      + (1 | site:sp)
#                      , data=dat
#                      , phylonm = c("sp","sp:site")
#                      , phylo = phy
#                      , phyloZ=phyZ
#                      , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
#                      , REML = FALSE
# )
# 
# print(summary(lme4fit2))

quit()

#debug(phylo_lmm)
#debug(modify_phylo_retrms)

if(numsite == "ss"){
	lme4fit <- phylo_lmm(Y ~ X + (1+X|sp)
		, data=dat
		, phylonm = c("sp","site:sp")
		, phylo = phy
		, phyloZ=phyZ
		, control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
		, REML = FALSE
	)
}

if(numsite == "ms"){
  
t4 <- proc.time()
    lme4fit <- phylo_lmm(new_y ~ X
#	 	+ (1 | sp)
#		+ (1 | obs)
      + (1 + X | sp)
		+ (1 + X | obs)
      + (1 | site)
		+ (1 | sp:site)
      , data=dat
      , phylonm = c("sp","sp:site")
      , phylo = phy
      , phyloZ=phyZ
      , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
      , REML = FALSE
    )
}
t2 <- proc.time()

lme4time <- t2 - t1 
print(lme4time)

print(summary(lme4fit))


lme4_list <- list(lme4fit, lme4time)

saveRDS(lme4_list, file=paste("datadir/lme4/lme4",numsite,size,tree_seed,"rds",sep="."))

# 
# library(MCMCglmm)
# 
# nitt <- 5e4 ## was 5e6
# inv.phylo <- inverseA(phy,nodes="TIPS",scale=FALSE)
# prior <- list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
# MCMC_time <- system.time(
#   MCMCglmm_fit <- MCMCglmm(Y~X,random=~sp+obs,
#                            family="gaussian",ginverse=list(sp=inv.phylo$Ainv),
#                            prior=prior,data=dat,nitt=nitt,burnin=100,
#                            thin=nitt/100,verbose=TRUE))
# 
# ss <- summary(MCMCglmm_fit)
# print(ss)
#rdnosave()




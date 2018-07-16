## compound symmetric via lme4 
library(ape)
library(Matrix)
library(lme4)
library(dplyr)

print(covmat)
print(cov2cor(covmat))

phyZ <- phylo.to.Z(phy)

dat <- (dat
  %>% mutate(obs = sp)
)	

print(head(dat))

tempmod <- phylo_lmm(Y ~ X
  + (1|sp)
  # + (0 + X|sp)
 + (1 + X|sp)
  + (1|sp:site)
  , data=dat
  , phylonm = c("sp","sp:site")
  , phylo = phy
  , phyloZ=phyZ
  , nsp = nsite
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
)

t1 <- sd.B0/sigma(tempmod)
t2 <- rho.B01*sd.B1/sigma(tempmod)
t3 <- sqrt((sd.B1/sigma(tempmod))^2 - t2^2)

new_y <- simulate(tempmod, newparams=list(theta=c(ss/sigma(tempmod)
	, t1
	, sd.B1/sigma(tempmod)
#   , t2
#   , t3
	)
  , beta = c(beta0,beta1)))
dat$new_y <- new_y[[1]]

lme4time <- system.time(
  lme4fit <- phylo_lmm(new_y ~ X
    + (1|sp)
    # + (0 + X | sp)
	 + (1 + X | sp)
	 + (1|sp:site)
    , data=dat
    , phylonm = c("sp","sp:site")
    , phylo = phy
    , phyloZ=phyZ
    , nsp = nsite
    , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  )
)

print(lme4time)

print(summary(lme4fit))

lme4_list <- list(lme4fit, lme4time)

saveRDS(lme4_list, file=paste("datadir/lme4_cs",size,seed,"rds",sep="."))

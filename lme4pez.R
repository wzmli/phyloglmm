#### Fitting phyloglmm via lme4

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

print(covmat)
print(cov2cor(covmat))

t1 <- proc.time()

phyZ <- phylo.to.Z(phy,stand=FALSE)

dat <- (dat
        %>% mutate(obs = sp
                   , site = site_name
                   , y_na = NA)
)	

t2 <- proc.time()
#debug(phylo_lmm)
#debug(modify_phylo_retrms)

if(numsite == "ss"){
  lme4fit <- phylo_lmm(Y ~ X + (1+X|sp)
                       , data=dat
                       , phylonm = c("sp","site:sp")
                       , phylo = phy
                       , phyloZ=phyZ
                       , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"
                       )
                       , REML = TRUE
  )
  t3 <- proc.time()
  lme4time <- t3-t1
}

if(numsite == "ms"){
  
  tempmod <- phylo_lmm(Y ~ X
                       + (1 | sp:site)
                       + (1 + X | sp)
                       + (1 + X | obs)
                       + (1 | site)
                       , data=dat
                       , phylonm = c("sp","sp:site")
                       , phylo = phy
                       , phyloZ=phyZ
                       # , nsp = nsite
                       , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                       , REML = TRUE
  )
  
  t1 <- sd.B0 
  t2 <- rho.B01*sd.B1 
  t3 <- sqrt(sd.B1^2 - t2^2)
  
  t4 <- sd.tip 
  t5 <- rho.slopetip*sd.slope 
  t6 <- sqrt(sd.slope^2 - t5^2)
  new_y <- simulate(tempmod
                    , newparams=list(theta=c(ss
                                             , t1, t2, t3
                                             , t4, t5, t6
                                             , sd.site
                    )/sd.resid
                    , beta = c(beta0,beta1)
                    , sigma = sd.resid))
  dat$new_y <- new_y[[1]]
  
  t4 <- proc.time()
  lme4fit <- phylo_lmm(new_y ~ X
                       + (1 | sp)
                       + (0 + X| sp)
                       + (1 | obs)
                       + (0 + X | obs)
                       + (1 | site)
                       + (1 | sp:site)
                       , data=dat
                       , phylonm = c("sp","sp:site")
                       , phylo = phy
                       , phyloZ=phyZ
                       , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                       , REML = TRUE
  )
  t5 <- proc.time()
  lme4time <- t2-t1 + (t5-t4)
}

print(lme4time)

print(summary(lme4fit))


lme4_list <- list(lme4fit, lme4time)

saveRDS(lme4_list, file=paste("datadir/lme4pez",numsite,size,seed,"rds",sep="."))

#rdnosave()

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
                     , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
							, REML = FALSE
)

t1 <- sd.B0/sigma(tempmod)
t2 <- rho.B01*sd.B1/sigma(tempmod)
t3 <- sqrt((sd.B1/sigma(tempmod))^2 - t2^2)

new_y <- simulate(tempmod, newparams=list(theta=c(ss/sigma(tempmod)
                                                  , t1
                                                  , sd.B1/sigma(tempmod)
                                                    , t2
                                                    , t3
)
, beta = c(beta0,beta1)),allow.new.levels=TRUE)
dat$new_y <- new_y[[1]]

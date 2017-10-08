## lme4 test single site simulation test

library(ape)
library(Matrix)
library(lme4)
library(dplyr)

targetname <- unlist(strsplit(rtargetname,"[.]"))
size <- targetname[2]

if(size == "small"){nspp <- 20}
if(size == "med"){nspp <- 100}
if(size == "large"){nspp <- 1000}


sd0vec_ss <- sd1vec_ss <- corvec_ss <- residvec_ss <- numeric(simnum)
for(i in 1:simnum){
  seed <- i
  source("simulate_tree.R",echo=FALSE)
  phyZ <- phylo.to.Z(phy)
  lme4fit <- phylo_lmm(Y ~ X + (1 + X |sp)
                                  , data=dat
                                  , phylonm = "sp"
                                  , phylo = phy
                                  , phyloZ = phyZ
                                  , control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  )
  cov_ssdf <- as.data.frame(lme4::VarCorr(lme4fit))
  sdcor_ss <- cov_ssdf$sdcor
  sd0vec_ss[i] <- sdcor_ss[1]
  sd1vec_ss[i] <- sdcor_ss[2]
  corvec_ss[i] <- sdcor_ss[3]
  residvec_ss[i] <- sdcor_ss[4]
}

df <- data.frame(sd0vec_ss,sd1vec_ss,corvec_ss,residvec_ss)
print(summary(df))

# print(summary(sd0vec_ss))
# print(summary(residvec_ss))



library(lme4)

dat <- readRDS("trans_dune.RDS")

mod1<- lmer(Y ~ 1 + (1 |site:sp) , data=dat
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
)


mod2<- lmer(Y ~ 1 + (1 |sp:site) , data=dat
            , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
)

all.equal(logLik(mod1),logLik(mod2))

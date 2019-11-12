library(dplyr)
library(Matrix)
library(lme4)

set.seed <- 1001

ngroup <- 5
nid <- 9
nrep <- 80

id <- rep(1:nid,each=ngroup*nrep)*10
groups <- rep(1:ngroup,nid*nrep)

group_sd <- 5 
interaction_sd <- 10
res_sd <- 1

## Group random effects
b_group <- rnorm(n=ngroup, mean=0, sd = group_sd)
y_group <- rep(b_group, each=nid*nrep)

## Interaction random effects
interaction_sdvec <- rep(interaction_sd, ngroup)
interaction_varmat <- tcrossprod(interaction_sdvec)
interaction_covmat <- interaction_varmat * Diagonal(ngroup)

interactionSigma <- kronecker(interaction_covmat, diag(nid))

image(Matrix(interactionSigma))

interaction_Cholesky <- Cholesky(interactionSigma)

beta0 <- 0

b_interaction <- sparseMVN::rmvn.sparse(n=1
  , mu=rep(beta0, each = ngroup*nid)
  , CH = interaction_Cholesky
  , prec = FALSE # use covariance_Cholesky when F
)

### Note: b_interaction has a specific order

## Residual 
noise <- rnorm(ngroup*nid*nrep,sd=res_sd)

Y <- noise + rep(c(b_interaction),each=nrep)
dd <- data.frame(id,groups,y=Y+y_group,check=id+groups)
dd2 <- dd %>% arrange(groups) %>% mutate(y2 = Y + y_group,y3=Y)

ff <- lmer(y3~(1|groups:id)
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  , data=dd2
)


ff2 <- lmer(y3~(1|id:groups)
           , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
           , data=dd2
)

ff3 <- lmer(y2~(1|groups)  +(1|groups:id)
           , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
           , data=dd2
)


ff4 <- lmer(y2~(1|groups)  +(1|id:groups)
            , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
            , data=dd2
)

summary(ff)
summary(ff2)
summary(ff3)
summary(ff4)


## We have to arrange by groups
interaction_covmat <- interaction_varmat * Matrix(diag(1:5))
interactionSigma <- kronecker(interaction_covmat, diag(nid))


image(interactionSigma)
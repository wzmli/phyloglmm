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

dd <- data.frame(id,groups)
## head(dd)  ## groups varies fastest
## Group random effects
b_group <- rnorm(n=ngroup, mean=0, sd = group_sd)
y_group <- rep(b_group, each=nid*nrep)
## rep(rep(b_group,each=nrep),nid)  ## this would be the matching order

## Interaction random effects
interaction_sdvec <- rep(interaction_sd, ngroup)
interaction_varmat <- tcrossprod(interaction_sdvec)
interaction_covmat <- interaction_varmat * Diagonal(ngroup)
## more natural? form the quadratic product of S*C*S'
## where S is the sd vector and C is the correlation matrix
## admittedly I often use outer(interaction_sdvec,interaction_sdvec)*C
## but the 'mathy' way would be interaction_sd %*% C %*% t(interaction_sd)
## (but we need to get the dimensions/vector orientations right
## maybe even define quadform(S,C) ?
## if we knew we had homogeneous variances
##  Diagonal(ngroup)*interaction_sd^2 ?

interactionSigma <- kronecker(interaction_covmat, diag(nid))

image(Matrix(interactionSigma))

interaction_Cholesky <- Cholesky(interactionSigma)
## this might sometimes be the long way round (i.e. if we can
## construct the Cholesky directly)
## ?Matrix::triangularMatrix (should be possible to construct this
##  directly)
## the diagonal case is trivial

beta0 <- 0

b_interaction <- sparseMVN::rmvn.sparse(n=1
  , mu=rep(beta0, each = ngroup*nid)  ## be careful if beta0 is non-trivial
  , CH = interaction_Cholesky
  , prec = FALSE # use covariance_Cholesky when F
)

### Note: b_interaction has a specific order

## Residual 
noise <- rnorm(ngroup*nid*nrep,sd=res_sd)

## suggest:
if (FALSE) {
    ## don't replicate b_group ... and then ...
    dd <- within(dd, {
        int <- interaction(group,id)  ## ???
        Y <- noise + b_interaction[int] + b_group[group]
    })
}
Y <- noise + rep(c(b_interaction),each=nrep)
dd <- data.frame(id,groups,y=Y+y_group,check=id+groups, y3=Y)
dd2 <- dd %>% arrange(groups) %>% mutate(y2 = Y + y_group,y3=Y)

tmpf <- function(d) {
    d <- d[c("groups","id","y3")]
    d <- with(d,d[order(groups,id),])
    rownames(d) <- NULL
    return(d)
}
all.equal(tmpf(dd),tmpf(dd2))
               
ff <- lmer(y3~(1|groups:id)
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  , data=dd2
)

ff2 <- lmer(y3~(1|id:groups)
           , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
           , data=dd2
)

## identical
ff_bad <- update(ff, data=dd)
ff2_bad <- update(ff2, data=dd)


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

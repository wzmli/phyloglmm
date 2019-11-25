library(dplyr)
library(Matrix)
library(lme4)
library(tidyverse)

set.seed <- 10012

ngroup <- 50
nid <- 100
nrep <- 1

id <- 1:nid
groups <- 1:ngroup

group_sd <- 5
interaction_sd <- 10
resid_sd <- 1

interactions <- interaction(groups,id)

dd <- (data.frame(ints = levels(interactions))
   %>% rowwise()
   %>% mutate(ints = as.character(ints)
      , group = unlist(strsplit(ints,"[.]"))[1]
      , id = unlist(strsplit(ints,"[.]"))[2]
      )  
)

group_df <- data.frame(group = as.character(1:ngroup)
  , y_group = rnorm(ngroup,mean=0,sd=group_sd)
)

quadform <- function(sd_mat,cormat){
  sd_mat %*% cormat %*% sd_mat
}

## Interaction random effects
interaction_sdvec <- diag(rep(interaction_sd,ngroup)) ## matrix(rep(interaction_sd, ngroup),ncol=1)
interaction_cormat <- Diagonal(ngroup)
interaction_covmat <- quadform(interaction_sdvec,interaction_cormat)

interactionSigma <- kronecker(interaction_covmat, diag(nid))
## What is this order? groups varying fastest or id varying fastest

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

interaction_df <- (data.frame(ints = levels(interactions))
   %>% mutate(ints = as.character(ints)
       , y_interaction = c(b_interaction)
   )
)

ddjoin <- (dd[rep(1:nrow(dd),each=nrep),]
   %>% left_join(.,interaction_df)
   %>% left_join(.,group_df)
   %>% ungroup()
   %>% mutate(noise = rnorm(nid*ngroup*nrep,sd=resid_sd)
      , y = noise + y_group + y_interaction
      , y_int = noise + y_interaction
   )
)


ff <- lmer(y_int~(1|group:id)
  , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
  , data=ddjoin
)

ff2 <- glmer(y_int~(1|id:group)
           , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
           , data=ddjoin
)


ff3 <- lmer(y~(1|group)  +(1|group:id)
           , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
           , data=ddjoin
)


ff4 <- lmer(y~(1|group)  +(1|id:group)
            , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
            , data=ddjoin
)

summary(ff)
summary(ff2)
summary(ff3)
summary(ff4)


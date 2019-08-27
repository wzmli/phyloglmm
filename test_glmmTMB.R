## setting parameters 

tree_seed <- 820
set.seed(tree_seed)

nsite <- 1
nrep <- 1 
nspp <- 500

# residual sd (set to zero for binary data)
sd.site <- 0
sd.resid <- 2
Xsd <- 0
sd.interaction <- 0.000001

# # fixed effects
beta0 <- 0
beta1 <- 0

# magnitude of phylogenetic random effects
physd.B0 <- 4
physd.B1 <- 0
phyrho.B01 <- 0

## species level random effects 
sd.B0 <- 0
sd.B1 <- 0
rho.B01 <- 0

source('simulate_tree.R', echo=TRUE)
## lme4 code to set up random effect matrix
source('new_phylo_setup.R', echo=TRUE)
source('glmmTMBhacked.R', echo=TRUE)

library(dplyr)
library(purrr)
library(glmmTMB)
library(lme4)
library(Matrix)

## help("image-methods")
sparsity <-  function(m) {
    tot_filled <- length(m@i)
    tot_dim <- prod(m@Dim)
    return(c(n_elem=tot_filled,size=tot_dim,sparsity=tot_filled/tot_dim))
}

ii <- function(m,ylab="species\n(=observation)",...) {
    s <- sparsity(m)
    image(m,
          sub=sprintf("n: %d; size: %d; sparsity: %1.3f",
                      s[1],s[2],s[3]),
          ...)
}

## construct branch matrix ( ~ 4 seconds)
phyZ_time <- system.time(phyZ <- phylo.to.Z(phy,stand=FALSE))


## create a copy of species index
## regressor method looks for (...|sp) and changes the random effects according to phyloZ
## correlation method constructs the random effect Z and if ``lambhack==TRUE'' it will get multiplied by chol(Sigma)
dat <- (dat
        %>% mutate(obs = sp
        )
)	

regressor_time <- system.time(
  regressor_mod <- glmmTMBhacked(Y ~ X  
    + (1 | sp)
    # + (0 + X |sp)
    , data=dat
    , phyloZ = phyZ
    , phylonm = c("sp", "sp:site")
    , doFit=TRUE
    , control=glmmTMBControl(optCtrl=list(trace=1,iter.max=1e5,eval.max=1e5))
    , REML = FALSE
    , lambhack=FALSE ## Lambda switch 
    ) 
)
ii(phyZ,xlab="branch")

chol_time <- system.time(cholphyt <- t(chol(phyZ %*% t(phyZ))))
str(cholphyt)
ii(cholphyt,xlab="species")

corr_time <- system.time(
  corr_mod <- glmmTMBhacked(Y ~ X  
    + (1 | obs)
    # + (0 + X |obs)
    , data=dat
    , phyloZ = cholphyt
    , phylonm = c("sp", "sp:site")
    , doFit=TRUE
    , control=glmmTMBControl(optCtrl=list(trace=1,iter.max=1e5,eval.max=1e5))
    , REML = FALSE
    , lambhack=TRUE
  ) 
)

support_time <- c(phyZ_time[3],chol_time[3])
fit_time <- c(regressor_time[3],corr_time[3])

sumtable <- data.frame(method=c("regressor","correlation")
                      , support = c("phyloZ", "cholesky")
                      , support_time
                      , fit_time
                      , total_time = support_time + fit_time
)

print(VarCorr(regressor_mod))
print(VarCorr(corr_mod))

print(sumtable)

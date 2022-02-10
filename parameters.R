### parameters for phylogenetic tree simulations
# tree_seed <- 1205
set.seed(tree_seed)


# seed <- 0122

nrep <- 1 ## only lme4 can do number reps atm

# residual variance (set to zero for binary data)
sd.resid <- 1
sd.site <- 5
Xsd <- 1

# sd.resid <- 0.000000001
# # fixed effects
beta0 <- 0
beta1 <- 0

# magnitude of phylogenetic random effects
physd.B0 <- 10
physd.B1 <- 10
phyrho.B01 <- 0.7

# physd.B1 <- 0
# phyrho.B01 <- 0

# magnitude of random effects

## sd.B0 <- 3
## sd.B1 <- 3
## rho.B01 <- 0.4

sd.B0 <- 0
sd.B1 <- 0
rho.B01 <- 0

# compound symmetric parameters
sd.interaction  <- 5 ## sp:site

# MCMC iterations

nitt <- 20000  ##2e4 ## was 5e6
stan_nitt <- 10000

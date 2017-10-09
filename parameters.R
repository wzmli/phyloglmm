### parameters for phylogenetic tree simulations

seed <- 100
simnum <- 100
# Generate simulated data for nspp species and nsite sites
nspp <- 1000

nsite <- 1

# residual variance (set to zero for binary data)
sd.resid <- 10
#
# # fixed effects
beta0 <- 0
beta1 <- 0

# magnitude of random effects
sd.B0 <- 2
sd.B1 <- 6

rho.B01 <- 0.0

# whether or not to include phylogenetic signal in B0 and B1
signal.B0 <- TRUE
signal.B1 <- TRUE

# MCMC iterations

nitt <- 2e6 ## was 5e6


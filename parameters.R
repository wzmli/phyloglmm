### parameters for phylogenetic tree simulations

seed <- 2830

# Generate simulated data for nspp species and nsite sites
nspp <- 500

#nsite <- 50
nsite <- 1


# residual variance (set to zero for binary data)
sd.resid <- 10
#
# # fixed effects
beta0 <- 0
beta1 <- 0

# magnitude of random effects
sd.B0 <- 4
sd.B1 <- 8

# whether or not to include phylogenetic signal in B0 and B1
signal.B0 <- TRUE

#signal.B1 <- TRUE
signal.B1 <- FALSE

# single site
#single.site <- FALSE
single.site <- TRUE

# MCMC iterations

nitt <- 1e4 ## was 5e6



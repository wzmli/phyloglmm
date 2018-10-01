### parameters for phylogenetic tree simulations

set.seed(seed)
# simnum <- 200
# Generate simulated data for nspp species and nsite sites

#nspp <- 50
#nsite <- 10

# residual variance (set to zero for binary data)
sd.resid <- 0.0000001
# sd.resid <- 0.000000001
# # fixed effects
beta0 <- 0
beta1 <- 0

# magnitude of phylogenetic random effects
sd.B0 <- 10
sd.B1 <- 2
rho.B01 <- 0.0

# magnitude of random effects

sd.site <- 8
sd.tip <- 100
sd.slope <- 3
rho.slopetip <- 0.0

# whether or not to include phylogenetic signal in B0 and B1
signal.B0 <- TRUE
signal.B1 <- TRUE

# compound symmetric parameters 
ss <- 6 ## sp:site 


# MCMC iterations

nitt <- 80000  ##2e4 ## was 5e6
stan_nitt <- 10000

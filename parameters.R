### parameters for phylogenetic tree simulations
set.seed(seed)
simnum <- 200
# Generate simulated data for nspp species and nsite sites
#nspp <- 20

nsite <- 20
#nsite <- 1

# residual variance (set to zero for binary data)
 sd.resid <- 10
# sd.resid <- 0.000000001
# # fixed effects
beta0 <- 0
beta1 <- 0

# magnitude of random effects
sd.B0 <- 4
sd.B1 <- 2
#sd.B1 <- 0.00000001

rho.B01 <- 0.7
rho.B01 <- 0

# whether or not to include phylogenetic signal in B0 and B1
signal.B0 <- TRUE
signal.B1 <- TRUE

# compound symmetric parameters 
ss <- 6 ## sp:site 


# MCMC iterations

nitt <- 8e4 ## was 5e6


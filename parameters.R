### parameters for phylogenetic tree simulations
library(shellpipes)

nrep <- 1 ## only lme4 can do number reps atm

if(exists("numsite")){
	if(numsite == "mms"){
		nrep <- 3
	}
}


# residual variance (set to zero for binary data)
sd.resid <- 1
sd.site <- 5
Xsd <- 1

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

## ? lots of messing about here
sd.B0 <- 3
sd.B1 <- 3
rho.B01 <- 0.4

# compound symmetric parameters
sd.interaction  <- 5 ## sp:site

# MCMC iterations

nitt <- 20000  ##2e4 ## was 5e6
stan_nitt <- 10000

saveEnvironment()

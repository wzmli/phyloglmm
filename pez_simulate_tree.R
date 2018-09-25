nspp <- 20
nsite <- 15

# Simulate a phylogenetic tree
phy <- compute.brlen(rtree(n = nspp), method = "Grafen", power = 1)

# Generate two environmental variables that are randomly distributed
# among sites. Only values of env1 are known in the dataset.
sd.env <- 1
env1 <- rnorm(nsite, sd=sd.env)
env2 <- rnorm(nsite, sd=sd.env)
env2 <- rep(env2, each=nspp)

# Generate two traits that differ phylogenetically among species. 
# Both traits govern the response of species to env1, 
# but only the values of trait1 are known in the dataset.
sd.trait <- 1
trait1 <- rTraitCont(phy, model = "BM", sigma = sd.trait)
trait2 <- rTraitCont(phy, model = "BM", sigma = sd.trait)
trait2 <- rep(trait2, times=nsite)

# Simulate a dataset. The terms 'f11*d$trait1*env1' and 'f21*d$trait2*env1'
# give the trait-by-environment interactions.
b0 <- 0
b1 <- 1
b2 <- 1
c1 <- 1
c2 <- 1
f11 <- .5
f21 <- .5
sd.e <- 1

Y.e <- rnorm(nspp*nsite, sd=sd.e)

d <- data.frame(sp=rep(phy$tip.label, times=nsite),
	site=rep(1:nsite, each=nspp), trait1=rep(trait1,
	times=nsite), env1=rep(env1, each=nspp))
d$Y <- b0 + b1*d$trait1 + b2*trait2 + c1*d$env1 
	+ c2*env2 + f11*d$trait1*d$env1 + f21*trait2*d$env1 
	+ Y.e

# Simulate a phylogeny that has a lot of phylogenetic signal (power = 1.3)
phy <- compute.brlen(rtree(n = nspp), method = "Grafen", power = 1.3)

# Simulate species means
sd.sp <- 1
mean.sp <- rTraitCont(phy, model = "BM", sigma=sd.sp^2)
# Replicate values of mean.sp over sites
Y.sp <- rep(mean.sp, times=nsite)

# Simulate site means
sd.site <- 1
mean.site <- rnorm(nsite, sd=sd.site)
# Replicate values of mean.site over sp
Y.site <- rep(mean.site, each=nspp)

# Compute the covariance matrix for phylogenetic attraction.
sd.attract <- 1
Vphy <- vcv(phy)

# Standardize Vphy to have determinant = 1.
Vphy <- Vphy/(det(Vphy)^(1/nspp))

# Construct the overall covariance matrix for phylogenetic attraction.
V <- kronecker(diag(nrow=nsite, ncol=nsite), Vphy)
Y.attract <- array(t(rmvnorm(n=1, sigma=sd.attract^2*V)))

# Simulate residual errors
sd.e <- 1
Y.e <- rnorm(nspp*nsite, sd=sd.e)

# Construct the dataset
d <- data.frame(sp=rep(phy$tip.label, times=nsite), site=rep(1:nsite, each=nspp))

# Simulate abundance data
d$Y <- Y.sp + Y.site + Y.attract + Y.e

# Analyze the model
communityPGLMM(Y ~ 1 
	+ (1|sp__) 
	+ (1|site) 
	+ (1|sp__@site), 
	data=d, tree=phy)

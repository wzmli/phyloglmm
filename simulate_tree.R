### simulate phylogenetic tree

library(ape)
## library(MASS)  ## for mvrnorm() ## don't load so we don't screw up dplyr::select()
set.seed(seed)

phy <- rtree(n = nspp)

## Compute branch lengths using other methods (optional)
phy <- compute.brlen(phy, method = "Grafen", power = 0.5)

# standardize the phylogenetic covariance matrix to have determinant 1 (optional)
Vphy <- vcv(phy)
# Vphy <- Vphy/(det(Vphy)^(1/nspp))

# Generate environmental site variable and standardize it
X <- matrix(1:nsite, nrow = 1, ncol = nsite)
X <- (X - mean(X))/sd(X)

cormat <- matrix(c(1,rho.B01,rho.B01,1),2,2)
sdvec <- c(sd.B0,sd.B1)
varmat <- sdvec %*% t(sdvec)
covmat <- varmat * cormat

print(covmat)

Sigma <- kronecker(covmat,Vphy)

#if((signal.B0==signal.B1) & (signal.B0 == FALSE)){
#	Sigma <- kronecker(covmat,diag(nspp))
#}

b.all <- MASS::mvrnorm(n=1
	, mu=rep(c(beta0,beta1),each=nspp)
	, Sigma=Sigma
)

b0 <- b.all[1:nspp] 
b1 <- b.all[(nspp+1):(2*nspp)]

y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp,
            ncol = nsite) + matrix(outer(b1, X), nrow = nspp, ncol = nsite)
e <- rnorm(nspp * nsite, sd = sd.resid) # add residual variance 
y <- y + matrix(e, nrow = nspp, ncol = nsite)
y <- matrix(y, nrow = nspp * nsite, ncol = 1)

Y <- y
Y <- matrix(Y, nrow = nspp, ncol = nsite)

# name the simulated species 1:nspp and sites 1:nsites
rownames(Y) <- 1:nspp
colnames(Y) <- 1:nsite


# Transform data matrices into "long" form, and generate a data frame
YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)

XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow =
               nspp * nsite, ncol = 1)

site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol =
                                           1)), nrow = nspp * nsite, ncol = 1)
sp <- matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp),
             nrow = nspp * nsite, ncol = 1)

dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))





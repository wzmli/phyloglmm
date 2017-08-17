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


if(nsite == 1){
	X <- rep(0,nspp)
#	signal.B1 <- FALSE
}
# Perform a Cholesky decomposition of Vphy. This is used to
# generate phylogenetic signal: a vector of independent normal random
# variables, when multiplied by the transpose of the Cholesky
# deposition of Vphy will have covariance matrix equal to Vphy.

iD <- t(chol(Vphy))

cormat <- matrix(c(1,rho.B01,rho.B01,1),2,2)

## we could set up the entire random-effect var-cov matrix
## if we segregate intercepts and slopes as separate blocks
## i.e. in the combined vector of random effects, b0 comes
## first, then b1, not intermingled
if (rho.B01==0) {
    b0mat <- if (signal.B0) Vphy else sd.B0^2*diag(nspp)
    b1mat <- if (signal.B1) Vphy else sd.B1^2*diag(nspp)
    Sigma <- Matrix::bdiag(b0mat,b1mat)
} else {
    Sigma <- kronecker(cormat,Vphy)
}

b.all <- MASS::mvrnorm(n=1,
              mu=rep(c(beta0,beta1),each=nspp),
              Sigma=Sigma)
b0 <- b.all[1:nspp]
b1 <- b.all[(nspp+1):(2*nspp)]

# Set up species-specific regression coefficients as random effects 
## MLi: CHANGE, we want bivariate normal if we want intercept slope correlation
## It is better if we can just set correlation = 0 instead of simulating them independently


## if (signal.B0 == TRUE) {
##     b0 <- beta0 + iD %*% rnorm(nspp, sd = sd.B0)
## } else {
##     b0 <- beta0 + rnorm(nspp, sd = sd.B0)
## }
## if (signal.B1 == TRUE) {
##     b1 <- beta1 + iD %*% rnorm(nspp, sd = sd.B1)
## } else {
##     b1 <- beta1 + rnorm(nspp, sd = sd.B1)
## }

    
## assume we have signal.B0 and signal.B1
## rho.B01 is the correlation

sdvec <- c(sd.B0,sd.B1)
Sigma <- outer(sdvec,sdvec)
B <- MASS::mvrnorm(nspp,c(0,0),Sigma)
dfun <- function(x,signal) {
    if (signal) (iD %*% x) else x
}
b0 <- b0 + dfun(B[,1],signal.B0)
b1 <- b1 + dfun(B[,2],signal.B1)

## test: if we simulate many sets of values and compute
## the correlation for c(b0,b1), we should see diagonal blocks
## (of size nspxnsp) that reflect the among-species correlations,
## and off-diagonal blocks that are rho.B01*diag(nsp)

## for example: simulating bivariate normal stuff

# Simulate species abundances among sites to give matrix Y that
# contains species in rows and sites in columns
y <- rep(b0, each=nsite)
y <- y + rep(b1, each=nsite) * rep(X, nspp)
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

noise <- rnorm(nspp*nsite,mean=0,sd=1)

dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp), noise = noise)

print(head(dat))




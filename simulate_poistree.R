### simulate phylogenetic tree

library(ape)
## library(MASS)  ## for mvrnorm() ## don't load so we don't screw up dplyr::select()
seed <- 5191
set.seed(seed)


nspp <- 500
phy <- rtree(n = nspp)

# standardize the phylogenetic covariance matrix to have determinant 1 (optional)
Vphy <- vcv(phy)
# Vphy <- Vphy/(det(Vphy)^(1/nspp))

# Generate environmental site variable and standardize it
X <- rnorm(n=nspp,sd=Xsd)

sd.B0=2
sd.B1=3

cormat <- matrix(c(1,rho.B01,rho.B01,1),2,2)
sdvec <- c(sd.B0,sd.B1)
varmat <- sdvec %*% t(sdvec)
covmat <- varmat * cormat

Sigma <- kronecker(covmat,Vphy)

b.all <- MASS::mvrnorm(n=1
                       , mu=rep(c(beta0,beta1),each=nspp)
                       , Sigma=Sigma
)

b0 <- b.all[1:nspp]
b1 <- b.all[(nspp+1):(2*nspp)]

res <- rnorm(length(b0),mean=0,sd=1)

mu <- exp(b0+b1*X + res)

print(summary(mu))

Y <- rpois(length(mu),mu)

dat <- data.frame(Y, X = X, sp = rownames(Vphy))
print(dat)




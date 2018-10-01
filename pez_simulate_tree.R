library(ape)

# Generate simulated data for nspp species and nsite sites
nspp <- 15
nsite <- 10

# residual variance (set to zero for binary data)
sd.resid <- 10

# fixed effects
beta0 <- 0
beta1 <- 0

# magnitude of random effects
sd.B0 <- 6
sd.B1 <- 4

# whether or not to include phylogenetic signal in B0 and B1
signal.B0 <- TRUE
signal.B1 <- TRUE

# simulate a phylogenetic tree
phy <- rtree(n = nspp)

# standardize the phylogenetic covariance matrix to have determinant 1
Vphy <- vcv(phy)
# Vphy <- Vphy/(det(Vphy)^(1/nspp))

# Generate environmental site variable
X <- matrix(1:nsite, nrow = 1, ncol = nsite)
X <- (X - mean(X))/sd(X)

# Perform a Cholesky decomposition of Vphy. This is used to
# generate phylogenetic signal: a vector of independent normal random
# variables, when multiplied by the transpose of the Cholesky
# deposition of Vphy will have covariance matrix equal to Vphy.

iD <- t(chol(Vphy))

# Set up species-specific regression coefficients as random effects
if (signal.B0 == TRUE) {
  b0 <- beta0 + iD %*% rnorm(nspp, sd = sd.B0) + rnorm(nspp, sd = sd.B0)
} else {
  b0 <- beta0 + rnorm(nspp, sd = sd.B0)
}
if (signal.B1 == TRUE) {
  b1 <- beta1 + iD %*% rnorm(nspp, sd = sd.B1) + rnorm(nspp, sd = sd.B1)
} else {
  b1 <- beta1 + rnorm(nspp, sd = sd.B1)
}

# Simulate species abundances among sites to give matrix Y that
# contains species in rows and sites in columns
y <- rep(b0, each=nsite)
y <- y + rep(b1, each=nsite) * rep(X, nspp)
y <- y + rnorm(nspp*nsite) #add some random 'error'
y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp,
            ncol = nsite) + matrix(outer(b1, X), nrow = nspp, ncol = nsite)
e <- rnorm(nspp * nsite, sd = sd.resid)
y <- y + matrix(e, nrow = nspp, ncol = nsite)
y <- matrix(y, nrow = nspp * nsite, ncol = 1)

# Y <- rbinom(n = length(y), size = 1, prob = exp(y)/(1 + exp(y)))
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
sp <- paste0("t", matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp),
                         nrow = nspp * nsite, ncol = 1))

dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))

# Random effects
# random intercept with species independent
# random intercept with species showing phylogenetic covariances
# random slope with species independent
# random slope with species showing phylogenetic covariances
# random effect for site
phyrfit <- phyr::communityPGLMM(Y ~ X + (1|sp__) + (X|sp__), data = dat, family = "gaussian",
               tree = phy, REML = TRUE)
source("new_phylo_setup.R")
phyZ <- phylo.to.Z(phy,stand=FALSE)
lme4fit <- phylo_lmm(Y ~ X
                     + (1 | sp)
                     + (0 + X | sp)
                     # + (1 | obs)
                     # + (0 + X | obs)
                     # + (1 | site)
                     # + (1 | sp:site)
                     , data = dat
                     , phylonm = c("sp", "sp:site")
                     , phylo = phy
                     , phyloZ = phyZ
                     , control = lmerControl(check.nobs.vs.nlev="ignore", check.nobs.vs.nRE = "ignore")
                     , REML = TRUE
)
print(summary(lme4fit))

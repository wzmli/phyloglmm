library(lme4)

set.seed(101)

sd.B0 <- 2
sd.B1 <- 4
rho.B01 <- 0.5
nid <- 100
nrep <- 10

X <- rnorm(nid*nrep, 0, 1)

cormat <- matrix(c(1, rho.B01, rho.B01, 1), 2, 2)

sdvec <- c(sd.B0, sd.B1)
varmat <- sdvec %*% t(sdvec)
covmat <- varmat * cormat

## id varies fastest; terms as blocks
bSigma <- kronecker(covmat,diag(nid))
## terms varies fastest; ids as blocks
# bSigma <- kronecker(diag(nid),covmat)
image(Matrix(bSigma))
dd <- data.frame(sd1=numeric(2000),sd2=numeric(2000),cc=numeric(2000),type=rep(c("single","double"),each=1000))

b1 <- MASS::mvrnorm(n=1
                   , mu=rep(c(0, 0), each =nid)
                   , Sigma = bSigma)
                   # , empirical = TRUE)

dd[i,"sd1"] <- sd(head(b1,nid))
dd[i,"sd2"] <- sd(tail(b1,nid))
dd[i,"cc"] <- cor(head(b1,nid),tail(b1,nid))

# ## extract values (terms as blocks)
# Y.re1 <- rep(head(b1, nid), nrep)    ## intercept: assume reps varies fastest, ids varies slowest
# X.re1 <- rep(tail(b1, nid), nrep)*X  ## slope*X: ditto

# cov(matrix(c(b1),ncol=2))
# cov(cbind(head(b1,nid),tail(b1,nid)))
# cov(b2)

set.seed(i)
b2 <- MASS::mvrnorm(n=nid
                   , mu=c(0, 0)
                   , Sigma = covmat
                   , empirical = TRUE)

# Y.re2 <- rep(b2[,1],nrep)
# X.re2 <- rep(b2[,2],nrep)*X

# cov(cbind(Y.re1,X.re1))
# cov(cbind(Y.re2,X.re2))


dd[1000+i,"sd1"] <- sd(b2[,1])
dd[1000+i,"sd2"] <- sd(b2[,2])
dd[1000+i,"cc"] <- cor(b2[,1],b2[,2])
}

# plot(b2[,1],b2[,2])
# points(head(b1,nid),tail(b1,nid),col=4)

e <- rnorm(nid*nrep,sd=10)
y1 <- Y.re1 + X.re1 + e
y2 <- Y.re2 + X.re2 + e

dd <- data.frame(y1,y2,id=1:nid,X)

ff1 <- lmer(y1~X+(1+X|id),data=dd,control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"))
ff2 <- refit(ff1,y2)
summary(ff1)
summary(ff2)

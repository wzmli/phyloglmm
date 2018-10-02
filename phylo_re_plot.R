## What are phylogenetic signals?
library(ape)
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtree)
library(cowplot)

set.seed(11)
nspp <- 10
nrep <- 5
phy <- rtree(n = nspp)

x <- rnorm(nspp*nrep, sd=1)

gg_tree <- (ggtree(phy)
  + geom_tiplab()
)

# print(plot(phy))

# print(phy$tip.label)

Vphy <- vcv(phy)

physd.y <- 20
physd.x <- 5

sd.y <- 0
sd.x <- 0 

sd.resid <- 0

phycormat <- matrix(c(1,1,1,1),nrow=2)
physdvec <- c(physd.y, physd.x)
phyvarmat <- physdvec %*% t(physdvec)
phycovmat <- phyvarmat * phycormat

physigma <- kronecker(phycovmat, Vphy) 

betas <- c(0,0)

b_phy <- MASS::mvrnorm(n=1, mu=rep(betas,each=nspp), Sigma=physigma)

Y.phy <- rep(head(b_phy,nspp), each = nrep)
X.phy <- rep(tail(b_phy,nspp), each = nrep) 

print(data.frame(Y.phy, X.phy))


cormat <- matrix(c(1,1,1,1),nrow=2)

print(cormat)
sdvec <- c(sd.y, sd.x)
varmat <- sdvec %*% t(sdvec)
covmat <- varmat * cormat

print(covmat)

sigma_b <- kronecker(covmat, diag(nspp)) 

b <- MASS::mvrnorm(n=1, mu=rep(betas,each=nspp), Sigma=sigma_b)

Y.int <- rep(head(b,nspp), each = nrep)

X.int <- rep(tail(b,nspp), each = nrep)

print(data.frame(Y.int,X.int))

Y.e <- rnorm(nspp*nrep, sd=sd.resid)
X.e <- rnorm(nspp*nrep, sd=1)

Y.phyint <- Y.phy + Y.int
X.phyint <- X.phy + X.int
X.phyint.e <- X.phyint + X.e
Y.phyint.e <- Y.phyint + Y.e
Y <- Y.phyint.e + X.phyint.e

reorder_tip <- phy$tip.label[nspp:1]

dat <- data.frame(sp = rep(rownames(Vphy), each = nrep)
  , Y.phy
  , Y.phyint
  , Y.phyint.e
  , Y
  , X.phy
  , X.phyint.e
)


print(data.frame(Y.phyint, X.phyint))

print(plot(Y.phyint, X.phyint))

dat <- dat %>% mutate(sp = factor(sp,levels=rownames(Vphy)), obs=sp)

dat2 <- (dat
  %>% gather(key = sd, value, -c(sp,obs,X.phyint.e, X.phy))
  %>% mutate(sd = factor(sd,levels=c("Y.phy","Y.phyint","Y.phyint.e","Y")))
)

# print(dat2)


dat2 <- dat2 %>% mutate(ll = paste(sp,1:(nspp*nrep),sep="_"))

print(dat2)

gg <- (ggplot(dat2, aes(x=sd, y=value, group=ll, colour=obs))
  + geom_point()
  + geom_line(alpha=0.5)
  + theme_bw()
  + theme(legend.position = "none")
)

print(gg_tree)

print(gg)

gg2 <- (ggplot(dat, aes(y=Y.phy, x="1", color=obs))
	+ geom_point()
	+ theme_bw()
)

print(gg2)

print(plot_grid(gg_tree, gg, gg2, nrow=1))

print(plot_grid(gg,gg2,nrow=1,rel_widths = c(1,0.5)))

gg3 <- (ggplot(dat, aes(x=X.phyint.e, y=Y, color=sp))
	+ geom_point()
	+ theme_bw()
)

print(gg3)
#print(plot_grid(gg_tree,gg,nrow=2))


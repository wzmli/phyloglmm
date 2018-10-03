library(shiny)
library(ggplot2)
library(ape)
library(cowplot)
library(gridExtra)
library(ggtree)
library(dplyr)
library(tidyr)

treeseed <- 62
set.seed(treeseed)
nspp <- 3
nrep <- 100

phy <- rtree(n = nspp)
print(plot(phy))
Vphy <- ape::vcv(phy)
physd <- 3


phycormat <- diag(nrep)
physdvec <- rep(physd,nrep)
phyvarmat <- physdvec %*% t(physdvec)
phycovmat <- phyvarmat * phycormat

physigma <- kronecker(phycovmat, Vphy) 

b_phy <- MASS::mvrnorm(n=1, mu=rep(0,nspp*nrep), Sigma=physigma)

dat <- data.frame(Y = b_phy, sp = rep(rownames(Vphy),nrep), sim=rep(1:nrep,each=nspp))
dat2 <- dat %>% mutate(sp = factor(sp, levels=rownames(Vphy)))

gg <- (ggplot(dat, aes(y=Y, x=sp, color=sp, group=sim)) 
  + geom_point() 
  + geom_line()
  + scale_x_discrete(limits=rownames(Vphy))
)

print(gg)

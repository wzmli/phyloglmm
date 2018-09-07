## Example for the paper 

library(ape)
seed <- 101
set.seed(seed)

nspp <- 3
phy <- rtree(n = nspp)
print(plot(phy, show.tip.label = FALSE))
print(edgelabels(c("L1",'L2','L3','L4')))
print(tiplabels(c("t1", "t2", "t3")))
print(nodelabels(c("N1", "N2")))
print(phy$edge)
Vphy <- vcv(phy)


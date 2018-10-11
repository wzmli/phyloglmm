## scaling VCV
library(ape)
library(MCMCglmm)

seed <- 202
set.seed(seed)
nspp <- 10

phy <- rtree(nspp)
plot(phy)
vc <- vcv(phy)
vc
print(det(vc))
# print(diag(vc))

uphy <- compute.brlen(phy, method = "Grafen", power = 1)
plot(uphy)
uvc <- vcv(uphy)
print(det(uvc))

inverseA(phy,scale=TRUE)
inverseA(phy,scale=FALSE)
aa <- inverseA(uphy)

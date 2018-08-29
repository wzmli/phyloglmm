## scaling VCV
library(ape)
library(MCMCglmm)

seed <- 202
set.seed(seed)
nspp <- 5

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

inverseA(phy,scale=FALSE)
inverseA(uphy)

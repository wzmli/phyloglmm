### Fitting phyloglmm with gls
library(ape)
library(nlme)

fit_gls <- gls(Y~X
  , data=dat
  # , correlation=corBrownian(phy=phy)
  , verbose=FALSE
)

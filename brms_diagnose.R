
bb <- readRDS("datadir/brms/brms.ss.small.1.rds")
library(shinystan)
library(rstan)
library(bayesplot)
launch_shinystan(bb[[1]])
## slow!!!!
## mcmc_pairs(bb[[1]],
##            off_diag_fun = c("hex", "hex"))

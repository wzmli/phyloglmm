# phyloglmm

## package instructions

```r
## current 'working' branch is 'refactor'
remotes::install_github("wzmli/phyloglmm/pkg", force = TRUE)
remotes::install_github("wzmli/phyloglmm/pkg@refactor")
```

Reproducible example from [here](https://github.com/wzmli/phyloglmm/issues/12) (need to download and unpack the `MRE.zip` file linked there).

```r
{
    data <- read.table("MRE/data2.txt") |>
        transform(species_PM2 = factor(species_PM2),
                  X1 = factor(X1,
                              levels = c("1-10", "11-100", "101-1000", "> 1000")))
  tree_MRE <- ape::read.tree("MRE/tree_MRE.txt")
  library(phyloglmm)
  tol <- sqrt(.Machine$double.eps)
  if (any(tree_MRE$edge.length<(-1*tol))) stop("negative edge lengths")
  tree_MRE$edge.length <- pmax(0, tree_MRE$edge.length)
  
  ## enforced match between
  system.time(phyloZ <- phylo.to.Z(tree_MRE))
  phyloZ <- phyloZ[levels(factor(data$species_PM2)), ]
}
MRE_m1_PHY <- phylo_glmmTMB(Y ~
                                X1 + X2 + X3 + X4 +
                                (1|species_PM2) + (1|rand1) + (1|rand3:rand2),
                          , data = data
                          , phyloZ = phyloZ
                          , phylonm = "species_PM2"
                          , REML = FALSE
                          , family = "gaussian"
                          , control = glmmTMBControl(rank_check = "skip")
                            )

```


## pipeline

It looks like

```
make fit.$platform.$numsite.$size.$seed.Rout
```

is supposed to generate a fit using:

- `$platform`: pez phyr lme4 brms gls glmmTMB phylolm ?
- `$numsite`: `ss` (single), `ms` (multiple)
- `$size`: `small`, `medium`, `large`, `xlarge`
- `$seed`: integer value

Fit uncorrelated intercept slope glmms using pez and lme4

## Todo

- single site simulation (problematic, not sure why it doesn't work)
- multivariate simulation to get correlation between slopes and intercepts
- extend example with MCMCglmm, glmmTMB, gls, PGLMM, brms
- what is caper? there is a function call pgls

# TODO

## Create simulates test cases for lme4 

### Single site
nspp = small = 20
med = 100
large = 500

sd.resid = very small or 10

sd.B0 = 4
sd.B1 = very small or 2

- single site intercept zero residual (gls, lme4)
-- single site intercept with residual (gls, lme4)

DONE

##### 

lme4 resid = 10
rho.B01 = 0.7
- phylogenetic signals with uncorrelated slopes
- phylogenetic signals with correlated slopes

## DONE

### Multiple sites
nsite = 20

- with phylogenetic signals
- phylogenetic signals by site with uncorrelated slopes
- phylogenetic signals by site with correlated slopes 

MISSING: resume pez med 134 

Now working on cs case where we ditto above 

finish phylogenetic signal with compound symmetric _sp:site_
done cs uncorrelated slope

There is a small problem with the algebra with cs simulate for thetas

Working on intercept and uncorrelated slope pez

Three sample sizes:
Small: 20 sp by 20 sites
Medium: 100 sp by 20 sites
Large: 1000 sp by 20 sites

Crazy case: 5000 by 20 sites
(Cannot alocate memory 10k +)

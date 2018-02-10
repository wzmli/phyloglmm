# TODO

## 9 Feb notes (BMB)

- comments about tip variance for paper: if we want no variance, we have a few ways of achieving that in principle - (1) blme strategy (or hack blme); (2) glmmTMB with dispersion = ~0; etc..  On the other hand, we don't lose very much by fitting with tip dispersion if there really isn't any, and most of the time there would be some ... (be careful, cf. Boettiger's comments about Pagel's lambda - he says "it doesn't make sense to treat tips differently", it *might*, but we should address this in any case; perhaps ask L. Wolkowich for opinions here?)

- `cowplot` for arranging (see `scales.R`; is there a way to enforce equal plot areas for horizontally arranged plots with different label widths?)

- Brownian motion/Gaussian process: very broadly we could quote Felsenstein's book or Felsenstein's PIC paper or pruning algorithm paper, but it's probably best to just describe it ourselves

- where is the *current* version of the manuscript/write-up? `ms.tex`


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


pez/lme4 small/med phylo correlated slope 



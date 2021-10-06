## plots

* tweak plots, go over all

## 

* rerun all from scratch?
* run pez/phyr longer?
* adjust convergence rules to account for effective sample size

* make rules for making dirs appropriately?
for i in pez phyr brms gls lme4 glmmTMB phylolm; do mkdir -p datadir/$i; done
make fit.pez.ss.small.1.Rout ## OK
make fit.phyr.ss.small.1.Rout ## OK
make fit.MCMCglmm.ss.small.1.Rout
make fit.brms.ss.small.1.Rout
make fit.gls.ss.small.1.Rout
make fit.lme4.ss.small.1.Rout

## problems

* pez: Warning message: arguments `tree` and `tree_site` are deprecated; please use `cov_ranef` instead. 

* phylolm: fails
* lme4: fails
  Subscript `phylonm` has size 2 but must be size 1.

Figure out how new machinery works with multiple phylo grouping vars!!

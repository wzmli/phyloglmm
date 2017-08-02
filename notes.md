Aug 2nd 2017

Multiple site simulation test seems to work really well on lme4, now we have to make sure single site simulations work from 2 weeks ago.
We better make sure we know how to simulate a single site tree properly before hacking away at fits.

MCMCglmm can estimate the residual variance but cannot do the phylo variance. 
glmmPQL and gls cannot fit intercept models, therefore, we need to bump up our simulation and simulate fixed slopes.



Nov 8th 2019

simulations appears to be stable enough, maybe pez is just not good?

Nov 6th 2019

- Need to add "..." to all phyloglmm functions so it can pass other lmer/glmmtmb arguements
- Kronecker product 

July 25th 2018

Standardize phyZ will run into problems with glmmTMB nlminb (NA/NaN function evaluation) when we have a lot of species (nspp > 300).
Otherwise, it fixes the "scaling" problem (yay!). 

July 16th 2018

brms can do everything lme4 does but really slow
- mcmcglmm cannot do nested/compound symmetric/site:sp interaction case without hack (not doing it).
- the problem with MCMCglmm is we cannot use ginverse with interaction

for now, let's run slow brms, pez and lme4

July 12th 2018

Testing what happens when we don't include random slopes, intercept.

Single observation per species
Simulate random intercept only
- random intercept and no residual (tip variation)
- residual only

Simulate random slope no intercept 
- tip 
- intercept
- slope
- tip + intercept 
- tip + slope
- intercept + slope
- tip + intercept + slope


Oct 1st 2017

Hacked lme4 object does not play nicely with "simulate" so we are reverting back to using pez-ish simulation method.

Currently simulating multiple sites but filtering to single site to test if phyloglmm works for single site. It appears to fail if we have larger residual sd than phylogenetic sd. 

- Single site does not appear to look good but multiple site estimation looks okay. Running on yushan now.

-- Results for 1000sp vs 100sites are done and they are great. 

- Create a smaller example to compare with pez. 



----------------------------------------------------------------------


Sept 6th 2017

Li, Ives and Waller 2017

Step 1: fit without phylo.signal and use forward selection on functional traits to get the _best model_

Step 2: forward selection of random effects by site with the functional traits selected from step 1

Step 3: add phylo.signal


----------------------------------------------------------------------



Aug 17th 2017

We figured out and implemented correlated random slope and intercept.
We (MLi) are still having trouble simulating the right data to test fitting machinary. 

TODO:

- The obvious/simple todo is to fix/implement the GLMM version so we can use other distributions in the GLMM family.
- Does the correlation even work? 
----------------------------------------------------------------------

Aug 4th 2017

We need to figure out changes inside glmmTMB object.
The new version of glmmTMB object does not have the pieces to hack the random effect matrix (to make our old code work).

Still fighting with intercept slope correlation, need to figure out the dimensions of Lambdat sparse matrix and hack accordingly. 
We should go back and see if we _actaully_ did multiple sites correctly.
The independent way (current working way) is simply stacking the diagonal sparse matrices together. 
What we actaully want for correlated sparse matrix is do them at the same time.



----------------------------------------------------------------------

Aug 3rd 2017

Adding the "extra" noise term seems to stablize and make things run on all platforms.
TMB seems to be working, but still need code to extract estimates.

TODO:
Try to implement brms.
MCMCglmm , change nu = 1

----------------------------------------------------------------------


Aug 2nd 2017

Multiple site simulation test seems to work really well on lme4, now we have to make sure single site simulations work from 2 weeks ago.
We better make sure we know how to simulate a single site tree properly before hacking away at fits.

MCMCglmm can estimate the residual variance but fail to get the phylo variance. 
glmmPQL and gls cannot fit intercept models, therefore, we need to bump up our simulation and simulate fixed slopes.
A simple solution is to add a _noise_ fixed effect to make glmmPQL and gls to fit, but correlation(brownian) does not appear to do anything. 


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


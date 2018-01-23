## Journal

## Jan 22nd 2018

MLi's pez nested hack worked. 
If we construct the _nested_ RE manually (i.e. applying the kronecker product ourselves instead of letting pez do it internally) and trick pez to think it is a non-nested term, it will do the cholesky decomposition step. 
As expected, the result matches lme4's perfectly. 

For the pez hack, we need to "debug(communityPGLMM.gaussian)" and replace the variables coded in pez_hack.R

MLi is going to try to _hack_ pez and trick the nested case by applying the kronecker product and feed it in as non-nested term.
MLi think this will trigger the cholesky decomposition step and hope it will give the same results as lme4.

MLi tried to do a _step-wise_ additative RE comarison between lme4 and pez. 
All _non-nested_ RE have a cholesky decomposition step in pez and the results are almost idenical as expected between the two platforms.
The nested case however does not have a cholesky decomposition step.

We want to figure out what is wrong with our _nested_ or _sp:site_ model and why lme4 is producing different results as pez.
MLi (and maybe BB) thinks the phylogenetic repulsion case is hard and low in priority. 
MLi have mixed feelings about phylogenetic repulsion because the examples in Li and Ives (2017) have _zero_ variance (BB thinks this is fine, but we need to dig deeper if we want to get into the details). 


## Jan 17th 2018

All _key_ simulations are complete. We are trying to fit phylglmm to real data from a few studies (mainly the studies that used pez and gls). For the _dunes_ dataset used in Li and Ives (2017), we noticed our estimates were off and lower likelihood. 
We realized phyZ is _unordered_ therefore it does not match the data. After we ordered phyZ and cov alphabetically, phylo_lmm matches pez. 

We realized the results don't match in the nested model. 
We suspect the problem is in _sp:site_. We ran one of the models without the nested term and it matches with pez. 
We should test this more carefully in a stepwise additive term method to see if we can match the simple models.

From now on, write down exactly what we need to do before we do it and keep a honest journal.

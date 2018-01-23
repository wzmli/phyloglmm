## Journal

## Jan 22nd 2018

We want to figure out what is wrong with our _nested_ or _sp:site_ model and why lme4 is producing different results as pez.
MLi (and maybe BB) thinks the phylogenetic repulsion case is hard and low in priority. 
MLi have mixed feelings about phylogenetic repulsion because the examples in Li and Ives (2017) have _zero_ variance (BB thinks this is fine, but we need to dig deeper if we want to get into the details). 

MLi tried to do a _step-wise_ additative comparison between lme4 and pez. 
All _non-nested_ RE have a cholesky decomposition step in pez and the results are almost identical as expected between the two platforms.
The nested case however does not have a cholesky decomposition step.

### TODO:

- Spend a few more minutes and try to understand the nested code 
- Start filling in the MS given we are certain/confident we can match pez with a real example.


## Jan 17th 2018

All _key_ simulations are complete. We are trying to fit phylglmm to real data from a few studies (mainly the studies that used pez and gls). For the _dunes_ dataset used in Li and Ives (2017), we noticed our estimates were off and lower likelihood. 
We realized phyZ is _unordered_ therefore it does not match the data. After we ordered phyZ and cov alphabetically, phylo_lmm matches pez. 

We realized the results don't match in the nested model. 
We suspect the problem is in _sp:site_. We ran one of the models without the nested term and it matches with pez. 
We should test this more carefully in a stepwise additive term method to see if we can match the simple models.

From now on, write down exactly what we need to do before we do it and keep a honest journal.

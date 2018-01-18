## Journal

## Jan 17th 2018

All _key_ simulations are complete. We are trying to fit phylglmm to real data from a few studies (mainly the studies that used pez and gls). For the _dunes_ dataset used in Li and Ives (2017), we noticed our estimates were off and lower likelihood. We realized phyZ is _unordered_ therefore it does not match the data. After we ordered phyZ and cov alphabetically, phylo_lmm matches pez. 

We realized the results don't match in the nested model. We suspect the problem is in _sp:site_. We ran one of the models without the nested term and it matches with pez. We should test this more carefully in a stepwise additive term method to see if we can match the simple models.

From now on, write down exactly what we need to do before we do it and keep a honest journal.

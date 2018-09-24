## Journal

## Sept 18

pez is problematic. 
We need to do compbren (forcing evo at the same location) and scale it in order to "control" fitting simulated parameter and match lme4.

PS: intercept random effects are fine, but not random slope

## August 28

Two _full_ models, single site vs multiple site.

Single site model will have 
1+X|sp + residual

mismatch model:
gls(glmmTMB proxy?) 
phylolm

match models:
lme4
glmmtmb
mcmcglmm
brms

multiple sites:
1+X|sp (1+X|taxa?!?) + 1|site + 1|site:sp + residual
pez
lme4
glmmtmb
brms 

1. efficiency (time)
2. data size 
3. goodness of fit?


Discussion 
Are we doing too much? 
AIC for reduce models and what assumptions and when is it okay to "neglict" random slope ... 
Diagnostic test for phylogenetic signal?



## August 27

- cannot use compute.brlen("Grafen") to make ultrametric trees because the "scaling" trick will not work (det(vcv) will be zero)
- cannot scale non-ultrametric trees (MCMCglmm)
- scale simulated tree at the first step to avoid all these problems.

## July 27th 2018

- I don't like having both journal and notes.

Stuff to report

- We more or less have an idea/ ball park of what happens when one assumes a _pure_ evolution/ errorless model. 
We need to be careful and standardize phyZ or else we will run into trouble as number of species increase. 
At the moment, only glmmTMB can fit this example, but nlminb will break when we try to fit more than 300 species.

- For bayesian methods, brms can do all cases (but really slow). 
MCMCglmm is a little faster but it cannot do nested/compound symmetric/site:sp interaction case without digging into the guts (not doing it). 
The reason MCMCglmm fails because it does not allow ginverse with interactions.
We can always (should definitely list all cases on the chart) find a model where all models can fit _easily_, but I am not comfortable with making _bold_ comparisons because matching platforms for a _fair_ comparison is hard (we can always do it and say it is hard).



## May 15th 2018

The order for phyZ and dataframe is more robust.
Do we want to stick phyZ before lFormula?
If we are going to stick it right before lFormula, we can even fix the order there, instead of mkBlist.

Going to clean up the repo and rerun everything.
I want to purrr some of the fast fitting steps, but simulating the phylogenetic tree still takes a long time.

GLS is very fragile. We have to simulate a single site and cannot simulate multiple site and filter for single site. Filtering for a single site will run into singlar problems (I don't know why), maybe because the way we are matrix multiplying in the tree simulation code.

## Mar 8th 2018

Now there is a switch in phylo.to.Z species branch function that does the equilvalent scaling.

## Mar 7th 2018

I managed to figure out how to matrix standardization works.
The two step _division_ through me off because it was probably a stability thing.
What it is doing is trying to computate a matrix scalar such that the determinate = 1. 
Knowing this fact, we can compute the matrix scalar and post hoc transform our sd estimate.
I have confirmed this works.

## Mar 6th 2018

For the single site simulations, I am going to simulate the full model and try to fit it using phylolm and gls and combine everything to the full model plot.
The plotting code is very ugly, I will clean it when the first draft of the ms is complete.

Most of the rough ideas are in the intro and methods.
Time to separate the plots and import them in the ms.

## Mar 5th 2018

sp:site works, but site:sp doesn't due to the _extra_ transpose problem. 
Let's put it up on the issues and deal with it later after the ms.


Still working on ms. 
Multiple observation now works for non-interaction random effects.
Need to figure out how to do interaction random effects.

## Feb 13th 2018

The scaling factor/standardizing factor is not trivial. 
The "pez" way is also a non-trivial scaling, thus, it is hard to match.
MLi is going to put together the paper and run some extra test to see what happens when we ignore phylogenetic correlations.


## Feb 12th 2018

MLi hacked around lFormula and succeed in reproducing the same results. 
At this point, it is way to hacky and it is still very "data frame" input dependent, but at the very least, the interaction problem is solved. 
The new version is much more robust when it comes to figuring out correlated slope case because all the hacks are before the KR product with the model frame.

The left over todos are:
- write the paper
- scaling and standarding phyZ (species branch matrix)
- can we do unbalance data now?


## Feb 9th 2018

That is false. We made a terrible mistake just like Ives (I think) and that is assuming we got the order correct. 
This is evident in the correlated slope example. 
We should not do the modification at the current location.
We need to do it _way_ upstream at the place where they use species in the model matrix (I think).
At the very least, we know what is wrong and know how to approach the _new_ (same old) problem.

I don't know what happened, but everything works.

## Feb 7th 2018

Need to rerun correlated slope multiple site lme4 models over night.
The plots for lme4 correlated slope model looks bad (i.e. the estimates are off and unable to estimate the correlation correctly). 
MLi went back and found the problem. 
The problem is we cannot use the kronecker product the same way as pez. 
The problem is fixed after we revert back to our original version but now it gives us a hint on the math behind it.


All the simulations (except for interaction REs) are complete. 
MLi is going to pretty the plots up a little and start putting together the ms.

## Feb 5th 2018

MLi is going to revert back to the previous version of phyloglmm and clean up the code a bit. 
We should rerun the pez simulations and make sure the models match exactly (as in the setups).
We should nail down the species level phylogenetic random effect and continue with the paper. 
MLi will work on the algebra while the simulations are running again.

## Jan 30th 2018

MLi: I can appreciate why the pez developer for putting the nested RE as a seperate case.
Whether we are right or not, we are confusing the heck of ourselves.
I think I need to roll back and not fuss to get pez's answer and rethink about the math.

We need to check carefully in our own phyloglmm code and see how we are building the interaction Zt.
Right now it is a mess. 

BMB suggest check for maxit and optim warnings.

Added a standardizing switch.

One last test, does standardizing cov(phylo) give a better fit?
ANS: The result is the same, we can even apply an appropriate post-hoc transformation.

The current phyloglmm is using REML, we have to change it appropriately to match pez.
I don't know the phy.to.Z equivalent of standardizing cov(phylo).

MLi has a hypothesis that pez's default nested code is fishy because it doesn't have a cholesky step and was not apart of Zt (i.e. it did not matrix multiple after applying the kronecker product step). 
The first thing MLi tried was comment out the Zt step and it turned out exactly as MLi predicted (see compare_pez.Rout)

MLi successfully replicated Li and Ives 2017 result.
There were three problems:
- max iteration and reltol were bad
- they used maximum likelihood instead of REML
- standardizing the variance-covariance matrix does matter
Everything matched perfectly after fixing these minor issues. 

## Jan 24th 2018

dune_lme4.Rout is the initial comparison of lme4 and pez with the _dune_ dataset in Li and Ives (2017).
The results are different between platforms and with the results presented in Li and Ives (MLi thinks it might be just package version or minor stuff).
MLi diagnosed the problem and it is the nested case that is causing the difference between the two platforms.
The main difference is pez is not doing the cholesky decomposition in the nested RE.

compare.Rout separates nested and non-nested RE and compare the two platforms. 
The non-nested case matches pretty well as expected. 
For the nested case, we can see the results are different. 
MLi wrote a hacked version that does the cholesky decomposition and it worked.


MLi made a separate file ( hacked_nested.R ) for the pez nested hack.
We are not rely/use this hack for the paper.
Now we are going to clean up the codes and separate the appropriate parts. 

## Jan 23rd 2018

MLi is somewhat confident in moving forward. 
We would like to understand what pez is doing for the nested case, but it shouldn't stop us from adding what we can do for the ms for now.
Once we figure out what is going on with the nested pez case, we need to redo the simulation for that part.


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

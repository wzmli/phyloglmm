## random slopes model discussion

### Intro
- the standard P model (e.g. PIC) is a correlated-residual model, i.e. we assume the brain~body relationship is universal; some taxa have higher-than-expected brain sizes, these deviations are phylog. correlated

- but it's entirely plausible that the relationship itself should be evolving, thus that the brain~body slope varies in a phylogenetically correlated way

- there is an existing literature, people have tried to do this already (REFS)

- random slopes models are not always appropriate, but they're probably relevant over a wider range of scenarios than people are currently thinking about [REF Schielzeth (? and Nakagawa) paper on random-slopes models in behaviour]

-----

- in principle they can be fitted whenever ... (? what are the conditions on identifiability?)  In a regular (grouped) mixed model we can only fit random slopes (closely analogous to randomized-block design ANOVA) when we have multiple values of the covariate measured within each grouping variable; otherwise random slope variation is confounded with ???.

- it may be easier to think about the random-slopes model in a strictly hierarchical setting (i.e., estimating different slopes for each family) - the PGLMM collapses to a standard random-slopes model, which would only be identifiable if ...

- however, in the PGLMM context, as long as there's variation in the predictor (e.g. body mass) among tips, there will (almost always?) be variation among taxa at some level, so random-slopes models will (almost always) be *theoretically* identifiable

- however (#2): how much data do we need in order to practically estimate the random slopes? Are we making a mistake (cf Schielzeth) by ignoring random slopes? What should we do when we don't have enough data? (cf. Barr et al 2013 "keep it maximal"; Bates, Vashisth, etc etc who don't like that, possibility of Bayesian approaches, etc etc etc)

## tip variability and confounding

- random-intercept model only
- no repeated measures within tips
- if we could constrain residual variance to zero

- compare cases:

($\sigma^2_r = 0$ vs. $\sigma^2_r > 0$) *
(1 obs per tip vs. >1 obs per tip) *
(1|taxon included in model vs. not included in model) *
(LMM (residual variance included) vs GLMM (not included))

when we have (1|taxon) or residual we might be doing the
equivalent of Pagel's lambda (cf. Boettiger)


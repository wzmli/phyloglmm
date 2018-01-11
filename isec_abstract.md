## Abstract for ISEC 2018

Analyses of the traits and distributions of many species within a community must allow for similarities among species that are determined by their shared evolutionary history: this area of ecological/evolutionary statistics has been termed _phylogenetic comparative analysis_.
In particular, when we estimate the effects
of species traits on their probability of occurrence in
a particular site, we must allow for correlation due to
phylogenetic relationships.
The primary approaches to this problem,
phylogenetic generalized least square regression (PGLS) and
phylogenetic generalized linear mixed models (PGLMM), have
been widely and effectively used.
Nevertheless, these methods still have 
constraints in both computational efficiency and the
range of models they can handle that limit
the size and complexity of systems we can analyze.
Here we propose an alternative method
that incorporates phylogenetic correlations
as a Gaussian process,
modeling the evolutionary changes along each branch
of a phylogenetic tree as independent
(possibly multivariate) random Normal deviates.
This approach is an order of magnitude faster than
existing implementations, and allows us to model some
processes 
(such as the correlated evolution of multiple traits)
that have previously been neglected because they do not
fit simply into existing analytic frameworks.
We demonstrate the method with
simulated phylogenies and evolutionary models of varying
complexity, as well as real data from several previous studies.
Although our algorithmic approach is general and could be
implemented in a wide range of computational platforms,
we implemented our approach using the lme4 R package (the most widely used package for fitting mixed effect models). We compare our results against two existing R packages (gls in nlme and communityPGLMM in pez) to explore the limitations of different methods and quantify simulation accuracy and computational efficiency.





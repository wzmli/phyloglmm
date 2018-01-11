## Abstract for ISEC 2018

Phylogenetic regression via independent constract and further generalizations to generalized least sqaure regression (PGLS) and generalized linear mixed models (PGLMM) are widely used to incorportate phylogenetic relationships in species trait-environmental base applications/problems.
Despite the choices of modelling strategies avaliable, there are many limitations in the current state of art of statistical analysis that researchers often constraint to limit or neglect phylogenetic relationships, i.e. high dimensional systems, computation power and etc. 
Here we propose an alternative, yet straightforward but little-used method of constructing phylogenetic correlations directly from phylogenetic tree by summing the evolutionary changes that occurred on all of the branches in the phylogeny in its past.  
We apply relatively simple phylogenetic comparative approaches to data from simulated phylogenies that incorportate various complexity of phylogenetic signals in the model. 
We implemented our approach in lme4 R package, the most widely used package for fitting mixed effect models, and compare against two existing R packages (gls in nlme and communityPGLMM in pez) to explore the limitations of different methods and quantify simulation acccury and computational efficency.





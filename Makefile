### testing phylog in lme4-verse

### Hooks 
current: target
-include target.mk


##################################################################

# make files and directories

Sources = Makefile .gitignore README.md sub.mk LICENSE.md 
include sub.mk
# include $(ms)/perl.def

Makefile:

##################################################################

Sources += $(wildcard *.R)

######################################################################

### simulate phylogenetic tree

simulate_tree.Rout: parameters.R simulate_tree.R
	$(run-R)

### Fitting

fit_pez.Rout: parameters.R simulate_tree.Rout fit_pez.R
	$(run-R)

fit_lme4.Rout: parameters.R simulate_tree.Rout phyloglmm_setup.R phyloglmm.R
	$(run-R)






### Makestuff



## Change this name to download a new version of the makestuff directory
# Makefile: start.makestuff

-include $(ms)/git.mk
-include $(ms)/visual.mk
-include $(ms)/pandoc.mk

-include $(ms)/wrapR.mk
-include $(ms)/flextex.mk


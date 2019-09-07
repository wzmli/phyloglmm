## This is phyloglmm, a project directory with a makestuff _submodule_
## makestuff/makemod.Makefile

current: target
-include target.mk

# -include makestuff/perl.def

######################################################################

# Content

# ms.tex
# main.tex
# notes.md
# todo.md
# journal.md

ms.pdf: main.tex phyloglmm_ms.tex
	pdflatex phyloglmm_ms
	bibtex phyloglmm_ms
	pdflatex phyloglmm_ms

##################################################################

Sources += $(wildcard *.R)

### debug Morgan's example

debug.Rout: phyloglmm_setup.Rout ./debug_examp/worked_example_phylolmm.rds ./debug_examp/worked_example_phylolmm.R
	$(run-R)


#2, parameters.R
#3, phyloglmm_setup.R


######################################################################

### simulate phylogenetic tree


#parameters.Rout:

phyloglmm_setup.Rout:

simulate_tree.Rout: names.R parameters.R simulate_tree.R
	$(run-R)

phylo_re_plot.Rout: phylo_re_plot.R
	$(run-R)

phytree_var_plot.Rout: phytree_var_plot.R
	$(run-R)

pez_simulate_tree.Rout: pez_simulate_tree.R
	$(run-R)

simulate_poistree.Rout: parameters.R simulate_poistree.R
	$(run-R)

fit.pois.Rout: simulate_poistree.Rout new_phylo_setup.R phyloglmm_pois.R
	$(run-R)

######################################################################

example.Rout: example.R
	$(run-R)

### Single site 

fit.gls.%.Rout: names.R parameters.R simulate_tree.R fit_gls.R
	$(run-R)

fit.gls.ss.xlarge.104.Rout: fit_gls.R

fit.gls.ss.large.7777.Rout: fit_gls.R
fit.glmmTMB.ss.large.1.Rout:

### phylolm

fit.phylolm.%.Rout: names.R parameters.R simulate_tree.R fit_phylolm.R
	$(run-R)

fit.phylolm.ss.xlarge.1.Rout: fit_phylolm.R

### lme4 can fit single and multiple sites

fit.lme4.%.Rout: names.R parameters.R simulate_tree.R new_phylo_setup.R phyloglmm.R
	$(run-R)

fit.lme4pez.%.Rout: names.R parameters.R simulate_tree.R new_phylo_setup.R lme4pez.R
	$(run-R)

fit.lme4.ss.xlarge.2.Rout: phyloglmm.R
fit.lme4.ms.xlarge.12.Rout: phyloglmm.R
fit.lme4.ms.small.1.Rout:
fit.lme4.ms.xlarge.1.Rout: phyloglmm.R
### tmb

lme4_profile.Rout: lme4_profile.R
	$(run-R)

fit.glmmTMB.%.Rout: names.R parameters.R simulate_tree.R new_phylo_setup.R glmmTMBhacked.R fit_tmb.R
	$(run-R)

fit.lme4.ss.large.2.Rout: fit_tmb.R
fit.glmmTMB.ms.large.62.Rout: fit_tmb.R
fit.glmmTMB.ss.xlarge.2.Rout: fit_tmb.R

### pez can only fit multiple sites

fit.pez.ms.small.1.Rout: fit_pez.R

### Multiple sites compound symmetric case

fit_cs.lme4.%.Rout: names.R parameters.R simulate_tree.R phyloglmm_setup.Rout fit_cs_lme4.R
	$(run-R)

fit.pez.%.Rout: names.R parameters.R simulate_tree.R new_phylo_setup.R fit_pez.R
	$(run-R)

fit.phyr.%.Rout: names.R parameters.R simulate_tree.R fit_phyr.R
	$(run-R)

fit.phyr.ms.small.1.Rout: fit_phyr.R

fit.pez.ms.small.1.Rout: fit_pez.R

fit.lme4.ms.small.2.Rout: phyloglmm.R

compare.pez.%.Rout: names.R parameters.R simulate_tree.R new_phylo_setup.R compare_pez.R
	$(run-R)

compare.pez.ms.small.1.Rout: compare_pez.R

fit.pez.ms.med.1.Rout: fit_pez.R
fit.lme4.ms.med.1.Rout: phyloglmm.R

fit.brms.%.Rout: names.R parameters.R simulate_tree.R new_phylo_setup.R  fit_brms.R
	$(run-R)

fit.brms.ss.small.1.Rout: fit_brms.R
### Collect and plot results

collect_gls.Rout: collect_gls.R
	$(run-R)

collect.Rout: collect.R
	$(run-R)

plot.Rout: ./datadir/collect.RDS collect_brms.Rout plot.R
	$(run-R)

ssplot.Rout: ssplot.R
	$(run-R)

msplot.Rout: msplot.R
	$(run-R)

csplot.Rout: csplot.R
	$(run-R)

compare_pez_plot.Rout: compare_pez_plot.R
	$(run-R)

peztest.Rout: new_phylo_setup.R peztest.R
	$(run-R)

Ignore += outline.html


### Fitting using other platforms (NEED TO FIX/CLEAN)

fit.MCMCglmm.%.Rout: names.R parameters.R simulate_tree.R MCMCglmmhacked.R fit_MCMCglmm.R
	$(run-R)

fit.MCMCglmm.ss.small.1.Rout: fit_MCMCglmm.R

fit_glmmPQL.Rout: parameters.Rout simulate_tree.Rout fit_glmmPQL.R
	$(run-R)


fit_tmb.Rout: parameters.R simulate_tree.R phyloglmm_setup.R new_phylo_setup.R glmmTMBhacked.R fit_tmb.R
	$(run-R)

fit_gls.Rout: parameters.R simulate_tree.R fit_gls.R
	$(run-R)

fit_lme4.Rout: parameters.R simulate_tree.R new_phylo_setup.R phyloglmm.R
	$(run-R)

#### Li et al 2017 examples

read_data.Rout: hacked_code/data_clean/dune_traits_Z.txt hacked_code/0_pkg_func.R hacked_code/1-data.R
	$(run-R)

phylosig.Rout: read_data.Rout hacked_code/3-phylosig.R
	$(run-R)

forward_selection.Rout: phylosig.Rout get_RE.R hacked_code/0_pkg_func.R hacked_code/2-forward_selection_fixed_first.R
	$(run-R)

analyses.Rout: forward_selection.Rout hacked_code/4-analyses.R
	$(run-R)


### fitting dune's example using lme4

dune_lme4.Rout: phylosig.Rout get_RE.R hacked_code/0_pkg_func.R new_phylo_setup.R hacked_nested.R dune_lme4.R
	$(run-R)

compare_pez.Rout: phylosig.Rout get_RE.R hacked_code/0_pkg_func.R phyloglmm_setup.R hacked_nested.R compare_pez.R
	$(run-R)

compare_pez_lme4.Rout: phylosig.Rout get_RE.R hacked_code/0_pkg_func.R phyloglmm_setup.R hacked_nested.R compare_pez_lme4.R
	$(run-R)


newsim.Rout: newsim.R
	$(run-R)


retest.Rout: retest.R
	$(run-R)

######################################################################

## This ¶ can be deleted once makestuff is set up
msrepo = https://github.com/dushoff
makestuff/Makefile: makestuff
makestuff:
	git submodule add -b master $(msrepo)/makestuff

######################################################################

### Content

######################################################################

### Makestuff

Sources += Makefile makestuff

ms = makestuff
Makefile: makestuff/Makefile

makestuff/%.mk: makestuff/Makefile ;
makestuff/Makefile:
	git submodule update -i

-include makestuff/os.mk
-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk

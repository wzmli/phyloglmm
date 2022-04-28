## This is phyloglmm, a resting screens project directory
## makestuff/project.Makefile

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

## Mike rules
phyloglmm_ms.pdf: main.tex phyloglmm_ms.tex
	pdflatex phyloglmm_ms
	bibtex phyloglmm_ms
	pdflatex phyloglmm_ms

## JD rule
Sources += main.tex phyloglmm_ms.tex
phyloglmm_ms.pdf: main.tex phyloglmm_ms.tex

phyloglmm_ms_resubmit.pdf: main.tex phyloglmm_ms.tex plot.Rout

##################################################################

Sources += $(wildcard *.R)

######################################################################

### simulate phylogenetic tree

#parameters.Rout:

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

fit.pois.Rout: simulate_poistree.Rout phyloglmm_pois.R
	$(run-R)

######################################################################

example.Rout: example.R
	$(run-R)

### Single site 

### gls

fit.gls.%.Rout: names.R parameters.R simulate_tree.R fit_gls.R
	$(run-R)

### phylolm

fit.phylolm.ss.small.1.Rout: fit_phylolm.R
fit.phylolm.%.Rout: names.R parameters.R simulate_tree.R fit_phylolm.R
	$(run-R)

### lme4 can fit single and multiple sites


# datadir/lme4/lme4.ss.small.1.Rout: phyloglmm.R parameters.R
datadir/lme4/lme4.%.rds: names.R parameters.R simulate_tree.R phyloglmm.R
	$(run-R)

fit.lme4.ms.large.1.Rout: phyloglmm.R
fit.lme4.%.Rout: names.R parameters.R simulate_tree.R phyloglmm.R
	$(run-R)

fit.lme4pez.%.Rout: names.R parameters.R simulate_tree.R new_phylo_setup.R lme4pez.R
	$(run-R)

### tmb

fit.glmmTMB.ss.large.1.Rout: fit_tmb.R
fit.glmmTMB.ms.large.1.Rout: fit_tmb.R
datadir/glmmTMB/glmmTMB.%.rds: names.R parameters.R simulate_tree.R fit_tmb.R
	$(run-R)

fit.glmmTMB.%.Rout: names.R parameters.R simulate_tree.R fit_tmb.R
	$(run-R)

fit.lme4test.%.Rout: names.R parameters.R simulate_tree.R new_phylo_setup.R phylolme4.R
	$(run-R)

collect_lme4test.Rout: collect_lme4test.R

lme4_profile.Rout: lme4_profile.R
	$(run-R)

### pez can only fit multiple sites

fit.pez.%.Rout: names.R parameters.R simulate_tree.R fit_pez.R
	$(run-R)

### Multiple sites compound symmetric case

fit_cs.lme4.%.Rout: names.R parameters.R simulate_tree.R phyloglmm_setup.Rout fit_cs_lme4.R
	$(run-R)


### phyr
fit.phyr.%.Rout: names.R parameters.R simulate_tree.R fit_phyr.R
	$(run-R)

compare.pez.%.Rout: names.R parameters.R simulate_tree.R new_phylo_setup.R compare_pez.R
	$(run-R)

## brms

fit.brms.ss.med.1.Rout: 
fit.brms.%.Rout: names.R parameters.R simulate_tree.R fit_brms.R
	$(run-R)

### Collect and plot results

collect_gls.Rout: collect_gls.R
	$(run-R)

collect.Rout: collect.R
	$(run-R)

plot.Rout: ./datadir/collect.RDS plot.R
	$(run-R)

compare_pez_plot.Rout: compare_pez_plot.R
	$(run-R)

peztest.Rout: peztest.R
	$(run-R)

Ignore += outline.html


### Fitting using other platforms (NEED TO FIX/CLEAN)

fit.MCMCglmm.ss.large.1.Rout: fit_MCMCglmm.R
fit.MCMCglmm.%.Rout: names.R parameters.R simulate_tree.R MCMCglmmhacked.R fit_MCMCglmm.R
	$(run-R)

fit.lme4.ss.large.1.Rout: phyloglmm.R
fit.MCMCglmm.ss.large.1.Rout: fit_MCMCglmm.R


fit_glmmPQL.Rout: parameters.Rout simulate_tree.Rout fit_glmmPQL.R
	$(run-R)


fit.brms.ss.small.1.Rout: fit_brms.R
fit_tmb.Rout: parameters.R simulate_tree.R fit_tmb.R
	$(run-R)

fit.gls.ss.small.1.Rout: fit_gls.R
fit_gls.Rout: parameters.R simulate_tree.R fit_gls.R
	$(run-R)

fit.lme4.ss.small.1.Rout: phyloglmm.R
fit_lme4.Rout: parameters.R simulate_tree.R phyloglmm.R
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

dune_lme4.Rout: phylosig.Rout get_RE.R hacked_code/0_pkg_func.R hacked_nested.R dune_lme4.R
	$(run-R)

compare_pez.Rout: phylosig.Rout get_RE.R hacked_code/0_pkg_func.R hacked_nested.R compare_pez.R
	$(run-R)

compare_pez_lme4.Rout: phylosig.Rout get_RE.R hacked_code/0_pkg_func.R hacked_nested.R compare_pez_lme4.R
	$(run-R)


newsim.Rout: newsim.R
	$(run-R)


retest.Rout: retest.R
	$(run-R)

kronecker_sim.Rout: kronecker_sim.R
	$(run-R)

kronecker_fit.Rout: Kronecker_fit.R
	$(run-R)

######################################################################

## Interaction example
phylointeraction.Rout: phylointeraction.R
	$(run-R)


### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls $@

-include makestuff/os.mk
-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk
-include makestuff/texdeps.mk
-include makestuff/wrapR.mk

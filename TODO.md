## plots

* tweak plots, go over all

##  ASK MIKE

- significant changes:
   - interactions tweak in simulate_tree ...???
   - add form=~sp in gls/corBrownian???
- lme4: Reconsider solution for phylosp/phylonm!

* rerun all from scratch?
* run pez/phyr longer?
* adjust convergence rules to account for effective sample size
* make rules for making dirs appropriately?
* targets ??? parallelize runs? stop on failure?

see `run_all` ...

## problems

perl -wf makestuff/wrapR/Rtrim.pl fit.pez.ms.xlarge.2.wrapR.rout > fit.pez.ms.xlarge.2.Rout
Can't open fit.pez.ms.xlarge.2.wrapR.rout: No such file or directory at makestuff/wrapR/Rtrim.pl line 4.
Use of uninitialized value $f in substitution (s///) at makestuff/wrapR/Rtrim.pl line 6.
Use of uninitialized value $f in substitution (s///) at makestuff/wrapR/Rtrim.pl line 7.
Use of uninitialized value $f in print at makestuff/wrapR/Rtrim.pl line 8.
/bin/mv -f   fit.pez.ms.xlarge.2.Rlog ./.fit.pez.ms.xlarge.2.Rlog
/bin/mv: cannot stat 'fit.pez.ms.xlarge.2.Rlog': No such file or directory
make: *** [Makefile:136: fit.pez.ms.xlarge.2.Rout] Error 1




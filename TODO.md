## new

- (ML) where does phylotree.png come from? I would like to tweak it to be a little bit 'skinnier' top to bottom (and be reproducible).  Could re-make it myself but ...
- (ML) check notation in full multi-group model and explanation
- finish cleaning up notation?


## misc

- write cover letter
- finish responses
- read over full ms

## plots

* should figure 5 be scale="free" instead?


## OLD
	
### "ask Mike" (all resolved?)

- significant changes:
   - add form=~sp in gls/corBrownian???
- lme4: Reconsider solution for phylosp/phylonm!
- improvements for brms, MCMCglmm priors etc. ??
   - why are RE priors so wide ... ?
   - pay attention to divergences etc.?
- is gls collection stuff in `collect.R` correct? resid=0, phylosig = sigma?
- is pez much slower than previously? Why?

* move 'pez' to 'slow' category?
* switch loop order to interleave model types etc. ? (seed rep as outermost loop?)
* make rules for making dirs appropriately?
* set up for SHARCnet/furrr ??
* targets ??? parallelize runs? stop on failure?
* make rules for .rds file, not .Rout?

* speed up collect?
see `run_all` ...

## problems/incomplete runs

* pez.ms missing: 33, 34, 41, 42, 49, 50
* pez.xlarge missing: all? (memory?)





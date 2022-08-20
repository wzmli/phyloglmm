## new

- (ML) 500-species results appear to be missing in the single-group speed comparison. Can we get them back?
	- MLi: will do them thanks!.
- (ML) do we have any idea why glmmTMB/lme4 times are now very similar (were we messing something up before, or ??)
	- MLi: This we need chat about this. 
- (ML) where does phylotree.png come from? I would like to tweak it to be a little bit 'skinnier' top to bottom (and be reproducible).  Could re-make it myself but ...
	- MLi: did it by hand, let's put it as the agenda item.
- define $\rho$ (correlation) in model definitions (shows up first in figures)
- consider `scale=free` for Figure 5?

- (BMB) tweak alignment of vector-matrix multiplication?
- (BMB/ML) finish cleaning up notation?


## misc
- write cover letter 
- finish responses
- read over full ms

## plots

* should figure 5 be scale="free" instead?


##  ASK MIKE

- significant changes:
   - add form=~sp in gls/corBrownian??? (done)
		- MLi: for correctness yes! For our example specifically, it doesn't matter. I ran it with a couple examples with form=~sp and without and it just returns identical output.  
- lme4: Reconsider solution for phylosp/phylonm!
	- MLi: ??
- improvements for brms, MCMCglmm priors etc. ??
   - why are RE priors so wide ... ? 
		- MLi: confused? They look reason to me. 
   - pay attention to divergences etc.? 
- is gls collection stuff in `collect.R` correct? resid=0, phylosig = sigma?
	- MLi: gls returns one parameter. 
- is pez much slower than previously? Why?
	- MLi: My guess is I reran everything on my laptop instead of my overclocked tower. So the slower platforms are more noticible. 

* move 'pez' to 'slow' category?
* switch loop order to interleave model types etc. ? (seed rep as outermost loop?)
* make rules for making dirs appropriately?
* set up for SHARCnet/furrr ??
* targets ??? parallelize runs? stop on failure?
* make rules for .rds file, not .Rout?
	- MLi: Yes! Will do. 

* speed up collect?
see `run_all` ...

## problems/incomplete runs

* pez.ms missing: 33, 34, 41, 42, 49, 50
* pez.xlarge missing: all? (memory?)
	- MLi: no xlarge for pez because it is way too slow. I don't recall the exact time, it was taking more than 2 or 3 hours for a single fit and I stopped. 



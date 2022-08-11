## plots

* tweak plots, go over all


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



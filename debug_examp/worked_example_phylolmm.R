####################################
#### Worked example, phylo_lmm. ####
####        Reminder:           ####       
#### 1) names in Z matrix       ####
#### 2) Lists in lme4 fit       ####
####################################

library(lme4)
library(ape)

phylo_examp <- readRDS(input_files[[1]])
try_phy <- phylo_examp[[1]]
titercurves_reduced <- phylo_examp[[2]]

## transform from phylogeny to Z matrix
brl_mat <- compute.brlen(try_phy)
phyZ <- phylo.to.Z(brl_mat)

## fit model
lme4time <- system.time(
	lme4fit <- phylo_lmm(
	  log(Titer) ~ poly(Day, degree = 2, raw = TRUE) + Log_Dose +
	  + (1 | Citation)
	  + (1 | unique_line)
		+ (1 | Scientific_Name)
		+ (0 + Day | Scientific_Name) 
		+ (0 + Daysq | Scientific_Name)	  
#	  + (0 + Daysq | Scientific_Name) 
#		+ (1+X|sp)
		, data = titercurves_reduced
		, phylonm = c("Scientific_Name")
		, phylo = try_phy
		, phyloZ = phyZ
		, nsp = 45
		, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
	)
)

## Random effect structure screwy (list filled in incorrectly)
debug(phylo_lmm)

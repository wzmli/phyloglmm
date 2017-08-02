#### Fitting using glmmPQL

library(MASS)
library(ape)

# rownames(dat) <- as.character(dat$phylo)

dat$allGrp <- factor(1) ## dummy grouping var because glmmPQL needs a group ...

fit_glmmPQL <- glmmPQL(Y~X
	, random = ~1|sp 
	, data = dat
	, family = "gaussian"
	, correlation = corBrownian(phy=phy)
	, verbose = FALSE
	)

summary(fit_glmmPQL)

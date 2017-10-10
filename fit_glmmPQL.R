#### Fitting using glmmPQL
library(dplyr)
library(MASS)
library(ape)

dat <- (dat
	%>% rowwise()
	%>% mutate(phylo=paste("t",sp,sep="")
		, obs=phylo
		)
)

dat <- data.frame(dat)
rownames(dat) <- as.character(dat$phylo)

# rownames(dat) <- as.character(dat$phylo)

dat$allGrp <- factor(1) ## dummy grouping var because glmmPQL needs a group ...

fit_glmmPQL <- glmmPQL(Y~1
	, random = list(sp=~1)
	, data = dat
	, family = "gaussian"
	, correlation = corBrownian(1,phy=phy)
	, verbose = FALSE
	)

print(summary(fit_glmmPQL))
print(VarCorr(fit_glmmPQL))

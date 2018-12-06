library(lme4)

ff <- list.files("datadir/lme4/", pattern = "ms.small")


confint_dat <- data.frame(seed = numeric(length(ff))
  , B0_lower = NA
	, B0_upper = NA 
	, B1_lower = NA
	, B1_upper = NA
)

i = 101
while(i < (length(ff)+1)){
	print(i)
	lme4obj <- readRDS(paste("datadir/lme4/", ff[i], sep=""))
	dat <- confint(lme4obj[[1]], method = "profile"
		, parm= c("(Intercept)","X")
		)
	confint_dat[i, "seed"] <- ff[i]
	confint_dat[i,"B0_lower"] <- dat["(Intercept)",1]
	confint_dat[i,"B0_upper"] <- dat["(Intercept)",2]
	confint_dat[i,"B1_lower"] <- dat["X",1]
	confint_dat[i,"B1_upper"] <- dat["X",2]
	i <- i + 1
}

print(confint_dat)
# rdsave(confint_dat)
saveRDS(confint_dat, "datadir/lme4_ms_small_profile.RDS")

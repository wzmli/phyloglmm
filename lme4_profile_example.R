## lme4 profile example 

lme4obj <- readRDS("datadir/lme4/lme4.ms.small.78.rds")

## check for convergence
lme4obj[[1]]@optinfo$conv$opt

## check for NaN/1/-1/...
print(lme4obj[[1]])

lme4_profile <- confint(lme4obj[[1]]
  , method = "profile"
  , parm = c("(Intercept)","X")
    )

pp <- profile(lme4obj[[1]],
              which = c("(Intercept)","X"),
              verbose=20)
confint(pp)

sessionInfo()

library(brms)
library(Matrix)
library(dplyr)
library(cmdstanr)

set_cmdstan_path("cmdstan/")

dat <- (dat
        %>% rowwise()
        %>% mutate(phylo=paste("t",sp,sep="")
                   , obs=phylo
        )
        #	 %>% filter(site == 1)
)

dat <- data.frame(dat)

inv.phylo <- MCMCglmm:::inverseA(phy,nodes="TIPS",scale=FALSE)

brms_time1 <- proc.time()
  A <- solve(inv.phylo$Ainv)
  rownames(A) <- rownames(inv.phylo$Ainv)
  brms_fit <- brm(y_main~ X + (1+X|gr(sp,cov=A))
      , family= gaussian()
      , data2 = list(A = A)
      , prior = c(prior(normal(0,1), "b")
			, prior(normal(0,1), "Intercept")
			, prior(normal(10,1), "sd")
                  ##			, prior(normal(phyrho.B01,1), "cor(Intercept,X)")
                  ## 1. source of divergent transitions if we use normal(1,1) [~ 1000 divergent trans]
                  ## 2. they go away with normal(10,1) but that's nowhere near the true value (1)
                  ##   (and the prior drags the estimate up to 10) [0 divergent transitions]
                  ## 3. they go away with normal(1, 0.001) (i.e. centered at true value with
                  ##  small SD/high precision) but that's cheating
                  ## next steps, *if* we wanted to take them
                  ##    - can we get away with some 0.001 < sigma < 1 ?
                  ##    - would allowing a heavier tail help? e.g. studentt(10, 1, 1) ???
                  ##    - look at diagnostics to try to figure out where/what is driving the divergences
                  ##      (the most typical problem, which may not apply here, is a 'funnel' - we would
                  ##       look at the bivariate distribution of the sigma estimate and any of the RE values
                  ##       I'm not sure how that applies here because we don't have latent variables??)
                  ##
                , prior(normal(1, 0.1), "sigma")
      )
      , data=dat
			, iter = stan_nitt
			, chains = 2
			# , adapt_delta = 0.8
		)


## fit once, for compilation
## might want to time solve() too, but will be a drop in the bucket
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

## Help with priors 
## https://rdrr.io/cran/brms/man/set_prior.html

brms_time1 <- proc.time()
  brms_dummy <- brm(y_main~ X + (1+X|gr(sp, cov = A))
                , family= gaussian()
                , data2 = list(A = A)
                  ## , cov_ranef = list(sp = A)
                , prior = c(prior(normal(0,1))
                #   , prior(normal(0,1), coef = "X")
                #   , prior(normal(20,2), class = sd, coef = "Intercept")
                #   , prior(normal(10,2), class = sd, coef = "X")
                #   , prior(normal(1,0.1), "sigma")
#                , prior(student_t(2,0,20), "sd")
#                , prior(student_t(4,0,20), "sigma")
              # , prior(normal(0,1),"b")          ## work
              # , prior(normal(0,2),"Intercept") ## work
              , prior(normal(10,1),"sd")        ## work
              , prior(normal(1,0.001),"sigma")    ## work
                  )
      , data=dat
      , iter = 10
      , chains = 1
      , backend="cmdstanr"
        )
brms_time2 <- proc.time()
brms_compiletime <- brms_time2 - brms_time1

## time only the sampling component
brms_time1 <- proc.time()
brms_fit <- update(brms_dummy
                   , iter = stan_nitt
                   , adapt_delta = 0.95
                 , chains = 2)
brms_time2 <- proc.time()

print(summary(brms_fit))

print(brms_fit)

brms_time <- brms_time2 - brms_time1
print(brms_time)

print(summary(brms_fit))

brms_list <- list(brms_fit, brms_time, brms_compiletime)
saveRDS(brms_list,file=paste("datadir/brms/brms",numsite,size,tree_seed,"rds",sep="."))

#rdnosave()

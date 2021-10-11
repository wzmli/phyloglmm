library(brms)
library(Matrix)
library(dplyr)

dat <- (dat
        %>% rowwise()
        %>% mutate(phylo=paste("t",sp,sep="")
                   , obs=phylo
        )
        #	 %>% filter(site == 1)
)

dat <- data.frame(dat)

inv.phylo <- MCMCglmm:::inverseA(phy,nodes="TIPS",scale=FALSE)


## fit once, for compilation
## might want to time solve() too, but will be a drop in the bucket
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
brms_time1 <- proc.time()
  brms_dummy <- brm(y_main~ X + (1+X|gr(sp, cov = A))
                , family= gaussian()
                , data2 = list(A = A)
                  ## , cov_ranef = list(sp = A)
                , prior = c(prior(normal(0,1), "b")
                , prior(normal(0,1), "Intercept")
                , prior(student_t(2,0,20), "sd")
                , prior(student_t(4,0,20), "sigma")
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
                 , chains = 2)
brms_time2 <- proc.time()

print(summary(brms_fit))

brms_time <- brms_time2 - brms_time1
print(brms_time)

print(summary(brms_fit))

brms_list <- list(brms_fit, brms_time, brms_compiletime)
saveRDS(brms_list,file=paste("datadir/brms/brms",numsite,size,tree_seed,"rds",sep="."))

#rdnosave()

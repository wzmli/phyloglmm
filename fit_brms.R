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

inv.phylo <- MCMCglmm:::inverseA(phy,nodes="TIPS",scale=TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

brms_time <- system.time(
    brms_fit <- brm(Y~ X + (1+X|sp)
      , family= gaussian()
      , cov_ranef = list(sp = A)
      , prior = c(prior(normal(0,1), "b")
			, prior(normal(0,1), "Intercept")
			, prior(student_t(2,0,20), "sd")
         , prior(student_t(4,0,20), "sigma")         
      )
      , data=dat)
)

print(summary(brms_fit))

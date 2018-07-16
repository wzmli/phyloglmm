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

brms_time1 <- proc.time()
  A <- solve(inv.phylo$Ainv)
  rownames(A) <- rownames(inv.phylo$Ainv)
  brms_fit <- brm(new_y~ X + (1+X|sp) + (1|site) + (1|site:sp)
      , family= gaussian()
      , cov_ranef = list(sp = A)
      , prior = c(prior(normal(0,1), "b")
			, prior(normal(0,1), "Intercept")
			, prior(student_t(2,0,20), "sd")
         , prior(student_t(4,0,20), "sigma")
      )
      , data=dat
			, iter = nitt
			, chains = 2
		)
brms_time2 <- proc.time()

print(summary(brms_fit))

brms_time <- brms_time2 - brms_time1
print(summary(brms_fit))

brms_list <- list(brms_fit,brms_time)
saveRDS(brms_list,file=paste("datadir/brms",size,seed,"rds",sep="."))

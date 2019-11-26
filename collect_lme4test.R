library(tidyverse)

reppath <- "./lme4test/reps2/"
noreppath <- "./lme4test/noreps2/"

msrep <- list.files(path=reppath, pattern= "rds")
msnorep <- list.files(path=noreppath, pattern = "rds")

turndat <- function(pp,x){
  ll <- readRDS(paste0(pp,x))
  seed <- c(strsplit(x,"[.]"))[[1]][4]
  size <- c(strsplit(x,"[.]"))[[1]][3]
  tt <- (tidy(ll)
    %>% mutate(seed=seed, size=size)
  )
}

ll_reps <- lapply(msrep,function(x)turndat(reppath,x))
repdat <- bind_rows(ll_reps)

dd1 <- (repdat
  %>% filter(!is.na(group))
  %>% mutate(par = paste(group,term,sep="_")
             , type = "rep")
)



ll_noreps <- lapply(msnorep,function(x)turndat(noreppath,x))
norepdat <- bind_rows(ll_noreps)

dd2 <- (norepdat
       %>% filter(!is.na(group))
       %>% mutate(par = paste(group,term,sep="_")
                  , type = "norep")
)

dd <- rbind(dd1,dd2)



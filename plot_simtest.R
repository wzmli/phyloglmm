## plotting simulation results

library(dplyr)
library(ggplot2)
library(reshape2)

ss_small <- readRDS("results/lme4_ss_simtest.small.RDS")

ss_med <- readRDS("results/lme4_ss_simtest.med.RDS")

ss_large <- readRDS("results/lme4_ss_simtest.large.RDS")

ms <- readRDS("results/ms_simtest.RDS")

mscs <- readRDS("results/mscs_simtest.RDS")

ss_small$size <- "small"
ss_med$size <- "medium"
ss_large$size <- "large"

ss_dat <- rbind(ss_small,ss_med,ss_large)

gg_sssd <- (ggplot(ss_dat,aes(x=size,y=sd0vec_ss))
  + geom_boxplot()
  + ggtitle("Single site simulation (1|sp)")
  + theme_bw()
  + geom_hline(yintercept =2)
)
print(gg_sssd)

gg_ssres <- (ggplot(ss_dat, aes(x=size, y=residvec_ss))
  + geom_boxplot()
  + ggtitle("Single site simulation (residual)")
  + theme_bw()
  + geom_hline(yintercept = 10)
)

print(gg_ssres)


## multiple sites

ms_pez_s <- ms[[1]] %>% mutate(dataSize = "small", platform = "pez")
ms_pez_m <- ms[[2]] %>% mutate(dataSize = "medium", platform = "pez")
ms_lme4_s <- ms[[3]] %>% mutate(dataSize = "small", platform = "lme4")
ms_lme4_m <- ms[[4]] %>% mutate(dataSize = "medium", platform = "lme4")
ms_lme4_l <- ms[[5]] %>% mutate(dataSize = "large", platform = "lme4")

ms_pez <- rbind(ms_pez_s,ms_pez_m)
ms_lme4 <- rbind(ms_lme4_s,ms_lme4_m,ms_lme4_l)

mpez <- melt(ms_pez,by=c("dataSize","platform"))
print(ggplot(mpez,aes(x=dataSize,y=value))+geom_boxplot()+facet_wrap(~variable,scales = "free") + ggtitle("ms_pez"))

mlme4 <- melt(ms_lme4,by=c("dataSize","platform"))
print(ggplot(mlme4,aes(x=dataSize,y=value))+geom_boxplot()+facet_wrap(~variable,scales = "free") + ggtitle("ms_lme4"))

## cs 

mscs_pez_s <- mscs[[1]] %>% mutate(dataSize = "small", platform = "pez")
mscs_pez_m <- mscs[[2]] %>% mutate(dataSize = "medium", platform = "pez")
mscs_lme4_s <- mscs[[3]] %>% mutate(dataSize = "small", platform = "lme4")
mscs_lme4_m <- mscs[[4]] %>% mutate(dataSize = "medium", platform = "lme4")
mscs_lme4_l <- mscs[[5]] %>% mutate(dataSize = "large", platform = "lme4")

mscs_pez <- rbind(mscs_pez_s,mscs_pez_m)
mscs_lme4 <- rbind(mscs_lme4_s,mscs_lme4_m,mscs_lme4_l)

mcspez <- melt(mscs_pez,by=c("dataSize","platform"))
print(ggplot(mcspez,aes(x=dataSize,y=value))+geom_boxplot()+facet_wrap(~variable,scales = "free") + ggtitle("mscs_pez"))

mcslme4 <- melt(mscs_lme4,by=c("dataSize","platform"))
print(ggplot(mcslme4,aes(x=dataSize,y=value))+geom_boxplot()+facet_wrap(~variable,scales = "free") + ggtitle("mscs_lme4"))


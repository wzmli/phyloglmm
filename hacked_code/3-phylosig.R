############################################################################################
# Question #4: phylogenetic signal
library(ape)
library(lme4)
library(phylolm)
library(picante)
library(dplyr)

traitnames <- names(dune.traits2)[-1]
# traitnames <- traitnames[c(3,6:16)]

trait.phy <- data.frame(year = c(2012), trait = traitnames, 
                        spnum = 0,
                        signal = 0, LR = 0, pvalue = 0, K = 0, 
                        Kpvalue = array(0, c(length(traitnames),1)))

counter <- 1

dat <- filter(dune.veg2, sp%in% dune.traits2$sp)

## hack: join traits to dat 

dat <- left_join(dat,dune.traits2,by="sp")

print(dat)


#   if (full.data == 1) {
#     data.sp <- aggregate(presence ~ sp, data = dat, sum)
#     #hist(data.sp$presence, breaks=30)
#     common.sp <- data.sp$sp[data.sp$presence >= 3]
#     
#     dat <- dat[is.element(dat$sp, common.sp), ]
#     dat$sp <- droplevels(dat$sp)
#     show("Only common species")
#     species.flag <- "common"
#   } else {
#     show("All species")
#     species.flag <- "all"
#   }

dat$sp <- as.factor(dat$sp)
dat$site <- as.factor(dat$site)

dat$presence <- (dat$freq > 0)
dat$Y <- log(dat$freq + 1)



for (itrait in traitnames) {
  # line up traits
  print(itrait)
  sub_trait <- dune.traits2
  rownames(sub_trait) <- as.character(sub_trait$sp)
  sub_trait$sp <- as.character(sub_trait$sp)
  
  for(i in 2:dim(sub_trait)[2]){
    if(class(sub_trait[, i]) == "character"){
      sub_trait[, i] = as.factor(sub_trait[,i])
    }
  }
  
  sub_trait$X <- sub_trait[,names(sub_trait) == itrait]
  sub_trait <- sub_trait[is.element(sub_trait$sp, levels(dat$sp)),]
  sub_trait = na.omit(select(sub_trait, sp, X))
  # sub_trait <- sub_trait[!is.na(sub_trait$X),]
  sub_trait$sp <- as.factor(sub_trait$sp)
  
  # This needs to be done for Q1 "common spcies"!
  sub_dat <- dat[is.element(dat$sp, sub_trait$sp),]
  sub_dat$sp <- droplevels(sub_dat$sp)
  
  nspp <- nlevels(sub_dat$sp)
  nsite <- nlevels(sub_dat$site)
  
  phy <- dune.phylo2
  
  # sort data set to match phy
  sub_trait <- sub_trait[match(phy$tip.label,as.character(sub_trait$sp)),]   
  # cbind(phy$tip.label, as.character(sub_trait$sp))
  
  # set up X
  sub_trait$X <- as.numeric(sub_trait$X)
  if(length(unique(sub_trait$X)) > 2) sub_trait$X <- (sub_trait$X - mean(sub_trait$X))/sd(sub_trait$X)
  # sub_trait$X <- (sub_trait$X - mean(sub_trait$X))/sd(sub_trait$X)
  
  if(length(unique(sub_trait$X))>2){
    z <- phylolm(X ~ 1, phy = phy, model = "lambda", data = sub_trait)
    sig <- z$optpar
    # z <- gls(X ~ 1, data = sub_trait, correlation = corPagel(1,phy), method = "ML")
    # sig <- attr(z$apVar,"Pars")[1]
    z0 <- lm(X ~ 1, data = sub_trait)
    
    LR <- logLik(z)[[1]] - logLik(z0)[[1]]
    pvalue <- pchisq(2*(logLik(z)[[1]] - logLik(z0)[[1]]), df=1, lower.tail = F)/2
  } else{
    catX <- (sub_trait$X == sub_trait$X[1])
    z <- binaryPGLMM(catX ~ 1, data = sub_trait, phy = phy)
    sig <- z$s2
    LR <- NA
    pvalue <- z$P.H0.s2
  }
  
  blomK = phylosignal(sub_trait$X, phy, reps = 10000)
  
  trait.phy$trait[counter] <- itrait
  trait.phy$spnum[counter] <- nlevels(sub_dat$sp)
  trait.phy$signal[counter] <- sig
  trait.phy$LR[counter] <- LR
  trait.phy$pvalue[counter] <- pvalue
  trait.phy$K[counter] = blomK$K
  trait.phy$Kpvalue[counter] <- blomK$PIC.variance.P
  counter <- counter + 1
}



trait.phy$signal = round(trait.phy$signal, 5)
trait.phy$pvalue = round(trait.phy$pvalue, 5)
trait.phy$bothSig = trait.phy$pvalue <= 0.05 & trait.phy$Kpvalue <=0.05
arrange(trait.phy, year, bothSig)
# write.csv(select(trait.phy, trait, signal, pvalue, K, Kpvalue), file = "phylo_signal.csv")

# write.csv(trait.phy, file = "trait_phylo_signal.csv", row.names =  FALSE)

# get traits with strong phylo-signal by date
trait_phy_sig = filter(trait.phy, year == "2012") %>% 
  filter(pvalue <= 0.05) %>%  .[["trait"]]
trait_phy_sig
# get traits that sites are selecting for


## Everything below doesn't work, comment out for now
# traits_lme4 = trait_selection(veg.long = dune.veg2, 
#                                            trait = dune.traits2,
#                                            binary = FALSE, trans =  "log")# %>% 
# traits_lme4 %>%  mutate(lmer.Pr = round(lmer.Pr, 6)) %>% 
#   arrange(lmer.Pr) # part of table 3
# 
# traits_glme4 = trait_selection(veg.long = dune.veg2, 
#                                             trait = dune.traits2, 
#                                             binary = TRUE)
# # part of table 3S

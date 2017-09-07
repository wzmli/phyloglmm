options(stringsAsFactors = F)
library(vegan)
library(stringr) # work with strings
library(dplyr) # data manipulation
library(brranching) # for pyhlomatic, to get the phylogeny
library(tidyr) # data manipulation
library(reshape2) # long to wide table
library(ape) # plot phylogeny
library(pez) # for communityPGLMM
library(lme4)
library(phylolm) # to test phylo signal of traits
library(picante)

# this function can be used to test for phylogenetic patterns (attraction/repulsion) 
#  it will make log-transformation for abundance data and then clean
#  all data (phylogeny and veg data) to have the same species.
# veg.long: the data set in long format, it must has these columns:
  # sp: species names
  # site: site names
  # freq: species abundance at each site
# phylo: the phylogeny to be used
# binary: logical, abundance or presence/absence data?
# trans: transformation, "log" for log-transformation.
phylo_pattern = function(veg.long, phylo, binary = FALSE, trans = NULL){
  # transformation of freq
  if(!is.null(trans)){
    if (trans == "log") {
      veg.long$Y <- log(veg.long$freq + 1)
    } 
  }
  # to make sure phylo has the right species
  phy <- drop.tip(phylo, tip = phylo$tip.label[
    !phylo$tip.label %in% unique(veg.long$sp)])
  # to make sure veg.long has the right species
  veg.long = filter(veg.long, sp %in% phy$tip.label)
  
  # sp and site need to be factors  
  veg.long$sp = as.factor(veg.long$sp)
  veg.long$site = as.factor(veg.long$site)
  nspp <- nlevels(veg.long$sp)
  nsite <- nlevels(veg.long$site)
  
  # to standardize the phylo var-covar matrix  
  Vphy <- vcv(phy)
  Vphy <- Vphy[order(phy$tip.label),order(phy$tip.label)]
  Vphy <- Vphy / max(Vphy)
  Vphy <- Vphy / det(Vphy) ^ (1 / nspp)
  
  show(c(nlevels(veg.long$sp), Ntip(phy)))
  if(nlevels(veg.long$sp) != Ntip(phy)){
    stop("The vegetation data and the phylogeny have different number of species")
  }
  
  # specify random effects
  re.site <- list(1, site = veg.long$site, covar = diag(nsite))
  re.sp <- list(1, sp = veg.long$sp, covar = diag(nspp))
  re.sp.phy <- list(1, sp = veg.long$sp, covar = Vphy)
  # sp is nested within site
  re.nested.phy <- list(1, sp = veg.long$sp, covar = Vphy, site = veg.long$site)
  re.nested.rep <- list(1, sp = veg.long$sp, covar = solve(Vphy), site = veg.long$site)
  
  if(!dir.exists("rds")) dir.create("rds")
  # to hold results as each model runs a while
  
  if(binary == FALSE){
    # this is the full model, with attraction to test
    z <- communityPGLMM(
      Y ~ 1, data = veg.long, family = "gaussian", sp = veg.long$sp, site = veg.long$site,
      random.effects = list(re.sp, re.sp.phy, re.site, re.nested.phy), 
      REML = F, verbose = F, s2.init = c(0.98, 0.001, 0.001, 0.0065), 
      reltol = 10^-5, maxit = 40)
    cat("sigma_c2_attraction = ", z$s2n, "\n")
    cat("attraction converge code = ", z$convcode, "\n")
    saveRDS(z, file = "rds/Q1_pattern_log_attraction.rds")
    
    # this is the model without nested term
    z0 <- communityPGLMM(
      Y ~ 1, data = veg.long, family = "gaussian", sp = veg.long$sp, site = veg.long$site,
      random.effects = list(re.sp, re.site, re.sp.phy), 
      REML = F, verbose = F, s2.init = c(0.98, 0.001, 0.001), 
      reltol = 10^-5, maxit = 40)
    saveRDS(z0, file = "rds/Q1_pattern_log_null.rds")
    cat("model without nested term converge code = ", z0$convcode, "\n")
    cat("p attraction = ", pchisq(2 * (z$logLik - z0$logLik), df = 1, lower.tail = F) / 2, "\n")
    
    # this is the full model, with repulsion to test
    z.rep <- communityPGLMM(
      Y ~ 1, data = veg.long, family = "gaussian", sp = veg.long$sp, site = veg.long$site, 
      random.effects = list(re.sp, re.sp.phy, re.site, re.nested.rep),
      REML = F, verbose = F, s2.init = c(0.98, 0.001, 0.028, 0.001), 
      reltol = 10^-5, maxit = 40)
    saveRDS(z.rep, file = "rds/Q1_pattern_log_repulsion.rds")
    cat("sigma_c2_repulsion = ", z.rep$s2n, "\n")
    cat("repulsion model converge code = ", z.rep$convcode, "\n")
    cat("p repulsion = ", pchisq(2 * (z.rep$logLik - z0$logLik), df = 1, lower.tail = F) / 2, "\n")
    
    sigma_output = data.frame(models = c("attraction", "repulsion", "no_nested"),
                              sp_num = rep(nspp, 3), 
                              sigma_sp = c(z$s2r[1], z.rep$s2r[1], z0$s2r[1]), 
                              sigma_sp_phy = c(z$s2r[2], z.rep$s2r[2], z0$s2r[2]), 
                              sigma_nested = c(z$s2n[1], z.rep$s2n[1], NA), 
                              sigma_site = c(z$s2r[3], z.rep$s2r[3], z0$s2r[3]), 
                              sigma_resid = c(z$s2resid, z.rep$s2resid, z0$s2resid),
                              p_nested = c(pchisq(2 * (z$logLik - z0$logLik), df = 1, lower.tail = F) / 2,
                                           pchisq(2 * (z.rep$logLik - z0$logLik), df = 1, lower.tail = F) / 2, 
                                           NA),
                              AIC = c(z$AIC, z.rep$AIC, z0$AIC))
  }
  
  if(binary == TRUE){
    veg.long$presence = veg.long$freq > 0
    # this is the full model, with attraction to test
    z <- communityPGLMM(
      presence ~ 1, data = veg.long, family = "binomial",
      sp = veg.long$sp, site = veg.long$site,
      random.effects = list(re.sp, re.sp.phy, re.site, re.nested.phy),
      REML = F, verbose = F, s2.init = c(2.8, 0.001, 0.01, 0.045),
      reltol = 10^-5, tol.pql = 10^-4)
    saveRDS(z, file = "rds/Q1_pattern_binary_attraction.rds")
    cat("sigma_c2_attract = ", z$s2n, "\n")
    cat("attraction converge code = ", z$convcode, "\n")
    w = communityPGLMM.binary.LRT(z, re.number = 4) # to test significance of nested term
    cat("p_attract = ", w$Pr, "\n")
    
    # this is the full model, with repulsion to test
    z.rep <- communityPGLMM(
      presence ~ 1, data = veg.long, family = "binomial",
      sp = veg.long$sp, site = veg.long$site,
      random.effects = list(re.sp, re.sp.phy, re.site, re.nested.rep),
      REML = F, verbose = F, s2.init = c(3, 0.001, 0.19, 0.001),
      reltol = 10^-5, tol.pql = 10^-4)
    saveRDS(z.rep, file = "rds/Q1_pattern_binary_repulsion.rds")
    cat("sigma_c2_repulsion = ", z.rep$s2n, "\n")
    cat("repulsion converge code = ", z.rep$convcode, "\n")
    w2 = communityPGLMM.binary.LRT(z.rep, re.number = 4)
    cat("p_repulsion = ", w2$Pr, "\n")
    
    # this is the model without nested term
    z0 = communityPGLMM(
      presence ~ 1, data = veg.long, family = "binomial",
      sp = veg.long$sp, site = veg.long$site, 
      random.effects = list(re.sp, re.sp.phy, re.site),
      REML = F, verbose = F, s2.init = c(3, 0.001, 0.19),
      reltol = 10^-5, tol.pql = 10^-4)
    saveRDS(z0, file = "rds/Q1_pattern_binary_null.rds")
    
    sigma_output = data.frame(models = c("attraction", "repulsion", "no_nested"),
                              spp_num = rep(nspp, 3),
                              sigma_sp = c(z$s2r[1], z.rep$s2r[1], z0$s2r[1]), 
                              sigma_sp_phy = c(z$s2r[2], z.rep$s2r[2], z0$s2r[2]), 
                              sigma_nested = c(z$s2n[1], z.rep$s2n[1], NA), 
                              sigma_site = c(z$s2r[3], z.rep$s2r[3], z0$s2r[3]),
                              p_nested = c(w$Pr, w2$Pr, NA))
  }  
  sigma_output
}

SigAIC <- function(mod, Penalty=qchisq(0.95,1)){   
  # function to calculate SigAIC value of the model
  # from Jamil et al. 2013 JVS
  LL <- logLik(mod) 
  ll <- as.numeric(LL)
  df <- attr(LL, "df")
  as.numeric(-2*ll + Penalty*df)
}

Model_select_tier1R<-function(start.model, block, sig.aic = FALSE){
  #  model selection tier 1, function to select the Random effects (traits)
  # from Jamil et al. 2013 JVS
  # terms = "0"
  Sigaic <- NA
  M <- vector("list", length = length(block))    
    formcc = as.character(formula(start.model))
  for(j in 1:length(block))  {
    print(block[j])
    if(any(str_detect(string = formcc, pattern = block[j]))){ # the trait already as fixed term
    M[[j]] <- update(start.model, as.formula(paste(". ~ . + (0 + ", block[j], "|site)")))
    } else { # trait not in the fixed term yet
    M[[j]] <- update(start.model, as.formula(paste(". ~ . +", block[j],  " + (0 + ", block[j], "|site)")))
    }
  }      
  if(sig.aic == TRUE) {
    LL1 <- unlist(lapply(M, SigAIC))
  } else {
    LL1 <- unlist(lapply(M, AIC))
  }
  k <- order(LL1, decreasing = FALSE)[1]
  cat("sigAIC after including term", block[k], LL1[k], "\n")
  # terms <-paste(terms, block[k], sep = "+")            
  list(B = block[-k], term = block[k], new.model = M[[k]], SigAic = LL1[k]) 
} 

Model_select_tier2F <- function(start.model, block, sig.aic = FALSE){
  #  model selection tier 2, function to select the trait in main effects
  M = list() 
  for(j in 1:length(block)){
    M[[j]] <-  update(start.model, as.formula(paste(". ~ . + ", block[j])))
  }
  if(sig.aic == TRUE) {
    LL1 <- unlist(lapply(M, SigAIC))
  } else {
    LL1 <- unlist(lapply(M, AIC))
  } 
  k <- order(LL1,decreasing=FALSE)[1] 
  cat("sigAIC after including term", block[k], LL1[k], "\n")
  list(B = block[-k], new.model = M[[k]], term = block[k], SigAic = LL1[k])
}

# this function will test how much of phylogenetic variation decreased after
# including functional traits in the model. You can specify traits as random terms
# by setting trait.re. This function will return the proportion of decreases. 
# It will also create a folder named as rds to hold resulted models, which can
# be used to extract details about models.
  # veg.long: veg data, same requirment as that in the phylo_pattern function
  # phylo: phylogeny
  # trait: trait by species matrix, traits as columns
  # trait.re: a vector of trait names for random terms. Default will use all traits in the trait matrix
  # binary: abundance or presence/absence?
  # trans: set to be "log" for abundance data
phylo_explained_by_multi_traits_re_sel = function(veg.long, phylo = pb.phylo, trait = pb.trait, 
                                                  trait.re = NULL,
                                                  binary = FALSE, trans = NULL){
  # transform frequency data
  veg.long$presence <- as.numeric(veg.long$freq > 0)
  
  if(!is.null(trans)){
    if(trans == "log") {veg.long$Y <- log(veg.long$freq + 1)}
  }
  
  trait = na.omit(trait)
  # make sure vegetation data has the right species
  veg.long = filter(veg.long, sp %in% trait$sp)
  
  # change characters into factors
  for(i in 2:dim(trait)[2]){
    if(class(trait[, i]) == "character"){
      trait[, i] = as.factor(trait[,i])
    }
  }
  
  ntraits = dim(trait)[2] - 1 # the first column is the species column: sp
  if(!is.null(trait.re)){
    which.re = which(names(trait)[-1] %in% trait.re)
  }
  # X1, X2 ... to make the following model set ups easier
  names(trait)[2:dim(trait)[2]] = paste0("X", 1:ntraits)
  for(i in 2:dim(trait)[2]){ # scale traits
    if(class(trait[, i]) == "numeric"){
      trait[,i] = (trait[,i] - mean(trait[,i]))/sd(trait[,i])
    }
  }
  
  # merge vegetation data with traits data
  dat = na.omit(left_join(veg.long, trait, by = "sp"))
  # sp and site need to be factors
  dat$sp <- as.factor(dat$sp)
  dat$site <- as.factor(dat$site)
  
  nspp <- length(unique(dat$sp))
  nsite <- nlevels(dat$site)
  
  # trim and standardize the phylogeny
  phy <- drop.tip(phylo, tip = phylo$tip.label[!phylo$tip.label %in% unique(dat$sp)])
  Vphy <- vcv(phy)
  Vphy <- Vphy[order(phy$tip.label),order(phy$tip.label)]
  Vphy <- Vphy/max(Vphy)
  Vphy <- Vphy/det(Vphy)^(1/nspp)
  
  show(c(nspp, Ntip(phy)))
  if(nspp != Ntip(phy)){
    stop("The vegetation data and the phylogeny have different number of species")
  }
  
  # random effect for site
  re.site <- list(1, site = dat$site, covar = diag(nsite))
  re.sp <- list(1, sp = dat$sp, covar = diag(nspp))
  re.sp.phy <- list(1, sp = dat$sp, covar = Vphy)
  
  re.nested.phy <- list(1, sp = dat$sp, covar = Vphy, site = dat$site)
  
  if(is.null(trait.re)){
    re.traits = vector("list", length = ntraits)
    for(j in 1:ntraits){
      re.traits[[j]] = list(unname(unlist(dat[paste0("X", j)])), site = dat$site, covar = diag(nsite))
    }
    re.all = vector("list", length = 4 + ntraits)
    re.all[[1]] = re.sp
    re.all[[2]] = re.sp.phy
    re.all[3:(2 + ntraits)] =  re.traits
    re.all[[(2 + ntraits + 1)]] = re.site
    re.all[[(2 + ntraits + 2)]] = re.nested.phy
  } else {
    ntraits = length(which.re)
    re.traits = vector("list", length = ntraits)
    for(j in seq_along(which.re)){
      re.traits[[j]] = list(unname(unlist(dat[paste0("X", which.re[j])])), site = dat$site, covar = diag(nsite))
    }
    re.all = vector("list", length = 4 + ntraits)
    re.all[[1]] = re.sp
    re.all[[2]] = re.sp.phy
    re.all[3:(2 + ntraits)] =  re.traits
    re.all[[(2 + ntraits + 1)]] = re.site
    re.all[[(2 + ntraits + 2)]] = re.nested.phy
  }
  str(re.all)
  
  if(!dir.exists("rds")) dir.create("rds")
  # to hold results as each model runs a while
  if(binary == FALSE){
    dat2 = select(dat, Y, starts_with("X"))
    # full model, with traits and phylo nested term
    z <- communityPGLMM(as.formula("Y ~ 1 + .") , data = dat2,
                        family = "gaussian", sp = dat$sp, site = dat$site, 
                        random.effects = re.all, REML = F, verbose = F,
                        reltol = 10^-5, maxit = 40,
                        s2.init = c(1.5, rep(0.01, (length(re.all)-1))))
    saveRDS(z, "rds/trait_multi_reg_log_full.rds")
    cat("s2_attract_with_traits = ", z$s2n[1], "\n")
    cat("trait+nested_phylo convcode = ", z$convcode, "\n")
    
    # model without traits, but with phylo nested term
    z0 <- communityPGLMM(as.formula("Y ~ 1 + .") , data = dat2,
                         family = "gaussian", sp = dat$sp, site = dat$site, 
                         random.effects = list(re.sp, re.sp.phy, re.site, re.nested.phy), 
                         REML = F, verbose = F, 
                         reltol = 10^-5, maxit = 40,
                         s2.init = c(1.5, .01, .01, .01)) # remove re.traits
    saveRDS(z0, "rds/trait_multi_reg_log_no_traits.rds")
    cat("s2_attract_without_traits = ", z0$s2n[1], "\n")
    cat("no_traits convcode = ", z0$convcode, "\n")
    
    # model without phylo nested term, but with traits
    z1 <- communityPGLMM(as.formula("Y ~ 1 + .") , data = dat2,
                         family = "gaussian", sp = dat$sp, site = dat$site, 
                         random.effects = re.all[-length(re.all)], REML = F, verbose = F, 
                         reltol = 10^-5, maxit = 40,
                         s2.init = c(1.5, rep(0.01, (length(re.all)-2)))) # remove re.nested.phy
    saveRDS(z1, "rds/trait_multi_reg_log_no_nested_phy.rds")
    cat("no_nested_phy convcode = ", z1$convcode, "\n")
    
    # model without traits and without phylo nested term
    z00 <- communityPGLMM(as.formula("Y ~ 1 + .") , data = dat2,
                          family = "gaussian", sp = dat$sp, site = dat$site, 
                          random.effects = list(re.sp, re.sp.phy, re.site), 
                          REML = F, verbose = F, 
                          reltol = 10^-5, maxit = 40,
                          s2.init = c(1.5, .01, .01)) # remove re.traits and re.nested.phy
    saveRDS(z00, "rds/trait_multi_reg_log_no_traits_no_nested.rds")
    
    convergent = data.frame(model = c("trait+nested_phylo", "no_traits", "no_nested_phy", "no_traits_nested_phylo"),
                            convcode = c(z$convcode, z0$convcode, z1$convcode, z00$convcode))
    write.csv(file = "Q4_log_multi_traits_convegent.csv", convergent)
    cat("no_traits_nested_phylo convcode = ", z00$convcode, "\n")
    
    q <- ntraits
    p_Xs <- 0
    for(k in 1:q)
      p_Xs <- p_Xs + 2^(-q) * choose(q,k) * 
      pchisq(2 * (z$logLik - z0$logLik), df = q, lower.tail = F)
    
    out_multi_reg_trait = t(data.frame(s2_attract_with_traits = z$s2n[1],
                                       p_attract_with_traits = pchisq(2 * (z$logLik - z1$logLik), df = 1, lower.tail = F)/2,
                                       s2_attract_without_traits = z0$s2n,
                                       p_attract_without_traits = pchisq(2 * (z0$logLik - z00$logLik), df = 1, lower.tail = F)/2,
                                       p_traits = p_Xs,
                                       s2_attract_desc = (z0$s2n -  z$s2n[1])/z0$s2n))
    write.csv(out_multi_reg_trait, file = "Q4_log_multi_traits_var_explained.csv")
  }
  
  if(binary == TRUE){
    dat2 = select(dat, presence, starts_with("X"))
    # full model, with traits and phylo nested term
    z <- communityPGLMM(as.formula("presence ~ 1 + .") , data = dat2,
                        family = "binomial", sp = dat$sp, site = dat$site, 
                        random.effects = re.all, REML = F, verbose = F, 
                        s2.init = c(1.5, rep(0.01, (length(re.all)-1))))
    saveRDS(z, "rds/binary_trait_multi_reg_full.rds")
    # to test the significance of phylo nested term with traits included
    w_nested_phy = communityPGLMM.binary.LRT(z, re.number = length(re.all))
    # to test the significance of traits
    w_traits = communityPGLMM.binary.LRT(z, re.number = c(3:(2 + ntraits)))
    cat("s2_attract_with_traits = ", z$s2n[1], "\n")
    cat("w_nested_phy = ", w_nested_phy$Pr, "\n")
    cat("w_traits = ", w_traits$Pr, "\n")
    
    #     z1 <- communityPGLMM(as.formula("presence ~ .") , data = dat2,
    #                          family = "binomial", sp = dat$sp, site = dat$site, 
    #                          random.effects = re.all[-length(re.all)], REML = F, verbose = F, 
    #                          s2.init = c(1.5, rep(0.01, (length(re.all)-2)))) # remove re.nested.phy
    #     saveRDS(z1, "binary_trait_multi_reg_no_nested_phy.rds")
    #   
    
    # model without traits, but with phylo nested term  
    z0 <- communityPGLMM(as.formula("presence ~ 1 + .") , data = dat2,
                         family = "binomial", sp = dat$sp, site = dat$site, 
                         random.effects = list(re.sp, re.sp.phy, re.site, re.nested.phy), 
                         REML = F, verbose = F, 
                         s2.init = c(1.5, .01, .01, .01)) # remove re.traits
    saveRDS(z0, "rds/binary_trait_multi_reg_no_traits.rds")
    cat("s2_attract_without_traits = ", z0$s2n[1], "\n")
    w_nested_phy_notraits = communityPGLMM.binary.LRT(z0, re.number = 4)
    
    convergent = data.frame(model = c("trait+nested_phylo", "no_traits"),
                            convcode = c(z$convcode, z0$convcode))
    write.csv(convergent, file = "Q4_binary_multi_traits_convegent.csv")
    
    q <- ntraits
    p_Xs <- 0
    for(k in 1:q)
      p_Xs <- p_Xs + 2^(-q) * choose(q,k) * pchisq(2 * w_traits$LR, df = q, lower.tail = F)
    
    #     
    #     z00 <- communityPGLMM(as.formula("presence ~ .") , data = dat2,
    #                           family = "binomial", sp = dat$sp, site = dat$site, 
    #                           random.effects = list(re.sp, re.sp.phy, re.site), 
    #                           REML = F, verbose = F, 
    #                           s2.init = c(1.5, .01, .01)) # remove re.traits and re.nested.phy
    #     saveRDS(z00, "binary_trait_multi_reg_no_traits_no_nested.rds")
    
    out_multi_reg_trait = t(data.frame(s2_attract_with_traits = z$s2n[1],
                                       p_attract_with_traits = w_nested_phy$Pr,
                                       s2_attract_without_traits = z0$s2n[1],
                                       p_attract_without_traits = w_nested_phy_notraits$Pr,
                                       p_traits = p_Xs,
                                       s2_attract_desc = (z0$s2n[1] - z$s2n[1])/z0$s2n[1]))
    write.csv(out_multi_reg_trait, file = "Q4_binary_multi_traits_var_explained.csv")
  }  
  return(out_multi_reg_trait)
}

# this function can be used to search for additional traits that can further reduce
# phylogenetic variation, after including traits selected by lmer.
selection = function(veg, trait, phylo, binary = FALSE, 
                     added.traits = NULL, fixed.terms = NULL){
  # removed already selected traits
  if(!is.null(added.traits)) {
    which.added = which(names(trait) %in% added.traits)
    trait2 = trait[-which.added]
  } else {
    trait2 = trait
  }
  
  if(!is.null(added.traits)){
    output_1 = data.frame(traits = names(trait2)[-1], 
                          s2_nested_w_t = NA, s2_nested_wo_t = NA, 
                          aic_w_t = NA, aic_wo_t = NA) %>% 
      mutate(traits = paste(paste(added.traits, collapse = "+"), traits, sep = "+"))
  } else{
    output_1 = data.frame(traits = names(trait2)[-1], 
                          s2_nested_w_t = NA, s2_nested_wo_t = NA, 
                          aic_w_t = NA, aic_wo_t = NA)
  }
  
  
  for (i in 2:dim(trait2)[2]){
    if(any(is.na(trait2[, i]))) cat(sum(is.na(trait2[, i])), 
                                    "sp do not have value for", 
                                    names(trait2)[i], "removed")
    
    if(!is.null(added.traits)) {
      # if(!is.null(fixed.terms)){
      #   mt = setdiff(fixed.terms, added.traits)
      #   if(names(trait2)[i] %in% mt){
      #     dat = left_join(veg, trait[, c(1, which.added)], by = "sp") %>% 
      #       left_join(trait2[, c(1, which(names(trait2) %in% mt))], by = "sp")
      #   } else{
      #     dat = left_join(veg, trait[, c(1, which.added)], by = "sp") %>% 
      #       left_join(trait2[, c(1, which(names(trait2) %in% mt))], by = "sp") %>% 
      #       left_join(trait2[, c(1, i)], by = "sp")
      #   }
      # } else{
      dat = left_join(veg, trait[, c(1, which.added)], by = "sp") %>%
        left_join(trait2[, c(1, i)], by = "sp")
      } else {
      dat = left_join(veg, trait2[, c(1, i)], by = "sp")
    }
    dat = na.omit(dat) 
    
    dat$sp = as.factor(dat$sp)
    dat$site = as.factor(dat$site)
    nspp = nlevels(dat$sp)
    nsite = nlevels(dat$site)
    
    # sub_trait = na.omit(trait[, c(1, i)])
    
    phy = drop.tip(phylo, tip = phylo$tip.label[
      !phylo$tip.label %in% unique(dat$sp)])
    Vphy = vcv(phy)
    Vphy = Vphy[order(phy$tip.label),order(phy$tip.label)]
    Vphy = Vphy/max(Vphy)
    Vphy = Vphy/exp(determinant(Vphy)$modulus/nspp)[1]
    
    show(c(nlevels(dat$sp), Ntip(phy)))
    
    re.site = list(1, site = dat$site, covar = diag(nsite))
    re.sp = list(1, sp = dat$sp, covar = diag(nspp))
    re.sp.phy = list(1, sp = dat$sp, covar = Vphy)
    re.nested.phy = list(1, sp = dat$sp, covar = Vphy, site = dat$site)
    # re.nested.rep = list(1, sp = dat$sp, covar = solve(Vphy), site = dat$site) 
    
    names(dat)[names(dat) == names(trait2)[i]] = "X"
    re.trait.to.add = list(as.matrix(dat["X"]), site = dat$site, covar = diag(nsite))
    
    if(!is.null(added.traits)){
      re.trait = vector("list", length = 1 + length(added.traits))
      for(j in 1:length(added.traits)){
        re.trait[[j]] = list(unname(unlist(dat[added.traits[j]])), site = dat$site, covar = diag(nsite))
      }
      re.trait[[1 + length(added.traits)]] = re.trait.to.add
      n.re.trait = 1 + length(added.traits)
    } else {
      re.trait = re.trait.to.add
      n.re.trait = 1
    }
    
    re.all = vector("list", length = 4 + n.re.trait)
    re.all[[1]] = re.sp
    re.all[[2]] = re.sp.phy
    if(n.re.trait == 1){
      re.all[[3]] =  re.trait
      re.all[[4]] = re.site
      re.all[[5]] = re.nested.phy
    } else{
      re.all[3:(2 + n.re.trait)] =  re.trait
      re.all[[(2 + n.re.trait + 1)]] = re.site
      re.all[[(2 + n.re.trait + 2)]] = re.nested.phy
    }
    str(re.all)
    
    if(!is.null(fixed.terms)){
      mt = setdiff(fixed.terms, added.traits)
      dat = left_join(dat, trait[c("sp", mt)], by = "sp")
      dat$sp = as.factor(dat$sp)
      dat$site= as.factor(dat$site)
      }
    
    if(binary == FALSE){
      dat$Y = log(dat$freq + 1)
      # full model
      if(!is.null(added.traits)){
        if(!is.null(fixed.terms)){
          if(names(trait2)[i] %in% mt){ # the random term added already in the fixed terms
            fm = paste("Y ~ 1", paste(fixed.terms, collapse = "+"), sep = "+")
          } else { # added random term not in the fixed terms, add it back
          fm = paste("Y ~ 1", paste(fixed.terms, collapse = "+"), "X", sep = "+")
          }
        } else { # no fixed terms specified, random terms only, put them back to fixed terms
        fm = paste("Y ~ 1", paste(added.traits, collapse = "+"), "X", sep = "+")
        }
      } else {# neither random nor fixed terms specified
        fm = "Y ~ 1 + X"
      }
      print(fm)
      z <- communityPGLMM(formula = as.formula(fm), data = dat, family = "gaussian", 
                          sp = dat$sp, site = dat$site, 
                          random.effects = re.all, 
                          REML = F, verbose = F, 
                          s2.init = c(1.5, rep(0.01, (length(re.all)-1))), 
                          reltol = 10^-5, maxit = 40)
      # saveRDS(z, file = paste0("select/", output_1$traits[i]), "_z.rds")
      # print(z$convcode)
      
      # no trait
      z0 <- communityPGLMM(as.formula(fm), data = dat, family = "gaussian", 
                           sp = dat$sp, site = dat$site, 
                           random.effects = re.all[c(1,2, length(re.all)-1, length(re.all))], 
                           REML = F, verbose = F, s2.init = c(1.5, .01, .1, .0064),
                           reltol = 10^-5, maxit = 40)
      # saveRDS(z0, file = paste0("select/", output_1$traits[i]), "_z0.rds")
      # print(z0$convcode)
    }
    
    if(binary == TRUE){
      dat$presence = as.numeric(dat$freq > 0)
      # full model
      if(!is.null(added.traits)){
        if(!is.null(fixed.terms)){
          if(names(trait2)[i] %in% mt){
            fm = paste("presence ~ 1", paste(fixed.terms, collapse = "+"), sep = "+")
          } else {
            fm = paste("presence ~ 1", paste(fixed.terms, collapse = "+"), "X", sep = "+")
          }
        } else {
          fm = paste("presence ~ 1", paste(added.traits, collapse = "+"), "X", sep = "+")
        }
      } else {
        fm = "presence ~ 1 + X"
      }
      print(fm)

      z <- communityPGLMM(formula = as.formula(fm), data = dat, family = "binomial", 
                          sp = dat$sp, site = dat$site, 
                          random.effects = re.all, 
                          REML = F, verbose = F, 
                          s2.init = c(1.5, rep(0.01, (length(re.all)-1))), 
                          reltol = 10^-5, maxit = 40)
      # saveRDS(z, file = paste0("select/", output_1$traits[i]), "_z_binary.rds")
      
      # no trait
      z0 <- communityPGLMM(as.formula(fm), data = dat, family = "binomial", 
                           sp = dat$sp, site = dat$site, 
                           random.effects = re.all[c(1,2, length(re.all)-1, length(re.all))], 
                           REML = F, verbose = F, s2.init = c(1.5, .01, .1, .0064),
                           reltol = 10^-5, maxit = 40)
      # saveRDS(z0, file = paste0("select/", output_1$traits[i]), "_z0_binary.rds")
    }
    
    output_1[i-1, 2] = z$s2n[1]
    output_1[i-1, 3] = z0$s2n[1]
    # output_1[i-1, 4] = z$AIC
    # output_1[i-1, 5] = z0$AIC
    # write.csv(output_1, file = "output_1_forward_selection.csv")
  }
  output_1 = mutate(output_1, prop = (s2_nested_wo_t - s2_nested_w_t)/s2_nested_wo_t) %>% 
    arrange(desc(prop)) # get the trait that decreases nested c^2 the most
  output_1
  return(output_1)
}

# this functional can be used to search environmental variables that closely related
# species have similar responses to. This may provide insight about future traits
# to measure.
# dat: long format for both veg and env data
# veg.long: used only if dat is NULL
# envi: environmental variables as columns, sites as rows
phylo_signal_slopes_envi = function(dat = NULL, veg.long = NULL, phylo, envi = NULL,  
                                    binary = FALSE){
  # remove singleton sp and include envi variables and scale them.
  if(is.null(dat)){
    veg.long = filter(veg.long, sp %in% phylo$tip.label)
    dat = left_join(veg.long, envi, by = "site")
  } else{
    dat = filter(dat, sp %in% phylo$tip.label)
  }
  
  dat$sp <- as.factor(dat$sp)
  dat$site <- as.factor(dat$site)
  
  nspp <- nlevels(dat$sp)
  nsite <- nlevels(dat$site)
  
  # transform frequency data
  dat$presence <- as.numeric(dat$freq > 0)
  dat$Y <- log(dat$freq + 1)
  
  # the covar matrix to be standardized to have determinant  of 1
  phy <- drop.tip(phylo, tip=phylo$tip.label[
    !phylo$tip.label %in% unique(dat$sp)])
  Vphy <- vcv(phylo)
  Vphy <- Vphy[order(phy$tip.label), order(phy$tip.label)]
  Vphy <- Vphy/max(Vphy)
  Vphy <- Vphy/det(Vphy)^(1/nspp)
  
  show(c(nspp, Ntip(phy)))
  if(nspp != Ntip(phy)){
    stop("The vegetation data and the phylogeny have different number of species")
  }
  
  # output
  envir <- data.frame(var=names(dat)[4:(dim(dat)[2] - 2)], 
                      lmer.logLik=NA, lmer.logLik0=NA, lmer.X.Pr=NA, 
                      PGLMM.intercept.star.s2=NA, PGLMM.intercept.phy.s2=NA, 
                      PGLMM.slope.star.s2=NA, PGLMM.slope.phy.s2=NA, 
                      PGLMM.resid.s2=NA, PGLMM.slopephy.logLik=NA,
                      PGLMM.slopephy.logLik0=NA, PGLMM.slope.phy.Pr=NA, 
                      PGLMM.slope.star.logLik=NA, PGLMM.slope.star.logLik0=NA, 
                      PGLMM.slope.star.Pr=NA, PGLMM.intercept.phy.Pr = NA,
                      PGLMM.B = NA, PGLMM.B.Pr = NA, PGLMM.B.se = NA)
  levels(envir$var) <- names(dat)[4:(dim(dat)[2] - 2)] # presence, Y at the end
  
  # start with the fourth column
  for (i in 4:(dim(dat)[2] - 2)) {
    print(i)
    if(binary == FALSE){
      # lmer
      z.1 <- lmer(Y ~ 1 + (1 | sp) + as.matrix(dat[, i]) + (0 + as.matrix(dat[, i]) | sp), 
                  data = dat, REML=F)
      z.0 <- lmer(Y ~ 1 + (1 | sp) + as.matrix(dat[, i]), data = dat, REML=F)
      envir[i-3,1] <- names(dat)[i]
      envir[i-3,2:4] <- c(logLik(z.1), logLik(z.0), 
                          pchisq(2 * (logLik(z.1) - logLik(z.0)), df=1, lower.tail=F)/2)
    }
    
    if(binary == TRUE){
      # lmer
      z.1 <- glmer(presence ~ 1 + (1 | sp) + as.matrix(dat[, i]) + (0 + as.matrix(dat[, i])| sp), 
                   data = dat, family = binomial)
      z.0 <- glmer(presence ~ 1 + (1 | sp) + as.matrix(dat[, i]), 
                   data = dat, family = binomial)
      envir[i-3,1] <- names(dat)[i]
      envir[i-3,2:4] <- c(logLik(z.1), logLik(z.0), 
                          pchisq(2 * (logLik(z.1) - logLik(z.0)), df=1, lower.tail=F)/2)
    }
    # re.1.site <- list(1, site = dat$site, covar = diag(nsite))
    re.1.star <- list(1, sp = dat$sp, covar = diag(nspp))
    re.1.phy <- list(1, sp = dat$sp, covar = Vphy)
    
    re.X.star <- list(as.matrix(dat[, i]), sp = dat$sp, covar = diag(nspp))
    re.X.phy <- list(as.matrix(dat[, i]), sp = dat$sp, covar = Vphy)
    
    # the goal is to test the significance of re.X.phy (phylo signal of slopes in abund~envi)
    if(binary == FALSE){
      z <- communityPGLMM(Y ~ 1 + as.matrix(dat[, i]), data = dat, family = "gaussian",
                          sp = dat$sp, site = dat$site, 
                          random.effects = list(re.1.star, re.1.phy, re.X.star, re.X.phy), 
                          REML = F, verbose = F, s2.init=.1)
      print(z$convcode)
      
      z0 <- communityPGLMM(Y ~ 1 + as.matrix(dat[, i]), data = dat, family = "gaussian", 
                           sp = dat$sp, site = dat$site, 
                           random.effects = list(re.1.star, re.1.phy, re.X.star), 
                           REML = F, verbose = F, s2.init=z$ss[c(1,2,3,5)]^2)
      print(z0$convcode)
      
      z00 <- communityPGLMM(Y ~ 1 + as.matrix(dat[, i]), data = dat, family = "gaussian", 
                            sp = dat$sp, site = dat$site, 
                            random.effects = list(re.1.star, re.1.phy), 
                            REML = F, verbose = F, s2.init=z$ss[c(1,2,3)]^2)
      print(z00$convcode)
      
      z000 <- communityPGLMM(Y ~ 1 + as.matrix(dat[, i]), data = dat, family = "gaussian", 
                             sp = dat$sp, site = dat$site, 
                             random.effects = list(re.1.star), 
                             REML = F, verbose = F, s2.init=z$ss[c(1,3)]^2)
      print(z000$convcode)
      
      envir[i-3,5:19] <- c(z$ss^2, 
                           z$logLik, 
                           z0$logLik, 
                           pchisq(2 * (z$logLik - z0$logLik), df=1, lower.tail=F)/2, 
                           z0$logLik, 
                           z00$logLik, 
                           pchisq(2 * (z0$logLik - z00$logLik), df=1, lower.tail=F)/2,
                           pchisq(2 * (z00$logLik - z000$logLik), df=1, lower.tail=F)/2,
                           z$B[2], z$B.pvalue[2], z$B.se[2])  
    }
    
    if(binary == TRUE){
      z <- communityPGLMM(presence ~ 1 + as.matrix(dat[, i]), data = dat, family = "binomial", 
                          sp = dat$sp, site = dat$site, 
                          random.effects = list(re.1.star, re.1.phy, re.X.star, re.X.phy), 
                          REML = F, verbose = F, s2.init = .1)
      print(z$convcode)
      w0 = communityPGLMM.binary.LRT(z, re.number = 4)
      z0 <- communityPGLMM(presence ~ 1 + as.matrix(dat[, i]), data = dat, family = "binomial", 
                           sp = dat$sp, site = dat$site, 
                           random.effects = list(re.1.star, re.1.phy, re.X.star), 
                           REML = F, verbose = F, s2.init = .1)
      print(z0$convcode)
      w00 = communityPGLMM.binary.LRT(z0, re.number = 3)
      #         z00 <- communityPGLMM(presence ~ 1 + dat[, i], data = dat, family = "binomial", 
      #                               sp = dat$sp, site = dat$site, 
      #                               random.effects = list(re.1.star, re.X.star), 
      #                               REML = F, verbose = F, s2.init=z$ss[c(1,3,5)]^2)
      envir[i-3,5:19] <- c(z$ss^2, NA, # no s2.resid
                           NA, # no logLik output
                           w0$LR, 
                           w0$Pr, 
                           NA, # no logLik output
                           w00$LR, 
                           w00$Pr,
                           NA, # intercept.phy not really interesting
                           z$B[2], z$B.pvalue[2], z$B.se[2])
      names(envir)[names(envir) == "PGLMM.slopephy.logLik0"] = "slopephy.logLik.minus.loglik0"
      names(envir)[names(envir) == "PGLMM.s2phy.logLik0"] = "s2phy.logLik.minus.loglik0"
    }
    write.csv(envir, file = "envir.csv")
  }
  
  envir$lmer.X.Pr = round(envir$lmer.X.Pr, 5)
  envir$PGLMM.slope.star.Pr = round(envir$PGLMM.slope.star.Pr, 5)
  envir$PGLMM.slope.phy.Pr = round(envir$PGLMM.slope.phy.Pr, 5)
  envir$PGLMM.intercept.phy.Pr = round(envir$PGLMM.intercept.phy.Pr, 5)
  envir$PGLMM.B.Pr = round(envir$PGLMM.B.Pr, 5)
  select(envir, var, lmer.X.Pr, PGLMM.slope.star.Pr, PGLMM.slope.phy.Pr, PGLMM.B, PGLMM.B.Pr) %>% 
    rename(slope.var.lmer.Pr = lmer.X.Pr,
           slope.var.indep.pglmm.Pr = PGLMM.slope.star.Pr, 
           slope.var.phylo.pglmm.Pr = PGLMM.slope.phy.Pr, 
           slope.estimated.pglmm = PGLMM.B, 
           slope.estimated.pglmm.Pr = PGLMM.B.Pr)
}

# this function to test whether trait-abundance relationships varied from site to site
trait_selection = function(veg.long, trait = pb.trait, 
                           binary = FALSE, trans = NULL){
  # transform frequency data
  veg.long$presence <- as.numeric(veg.long$freq > 0)
  if(!is.null(trans)){
    if(trans == "log") {veg.long$Y <- log(veg.long$freq + 1)}
  }
  veg.long$sp <- as.factor(veg.long$sp)
  veg.long$site <- as.factor(veg.long$site)
  
  for(i in 2:dim(trait)[2]){
    if(class(trait[, i]) == "character"){
      trait[, i] = as.factor(trait[,i])
    }
  }
  
  counter <- 1
  traitnames = names(trait)[-1]
  output <- data.frame(trait = traitnames, 
                       sp_richness = NA, dep = c("freq"), 
                       lmer.logLik=NA, lmer.logLik0=NA, lmer.Pr=NA)
  
  for (itrait in traitnames) {
    # line up traits
    sub_trait <- trait[is.element(trait$sp, veg.long$sp),]
    sub_trait$X <- sub_trait[, names(sub_trait) == itrait]
    sub_trait <- sub_trait[!is.na(sub_trait$X),] # rm missing values
    
    # This needs to be done for Q1 "common spcies"!
    sub_dat <- veg.long[is.element(veg.long$sp, sub_trait$sp),] # have same sp in traits data
    sub_dat$sp = droplevels(sub_dat$sp)
    sub_dat$site = droplevels(sub_dat$site)
    
    nspp <- nlevels(sub_dat$sp)
    nsite <- nlevels(sub_dat$site)
    
    # set up Vtrait
    sub_trait$X <- as.numeric(sub_trait$X)
    sub_trait$X <- (sub_trait$X - mean(sub_trait$X))/sd(sub_trait$X)
    Vtrait <- outer(sub_trait$X, sub_trait$X)
    
    sub_dat$X <- sub_trait$X
    show(c(nlevels(sub_dat$sp), dim(Vtrait)))
    
    if(binary == FALSE){
      # lmer
      z.1 <- lmer(Y ~ 1 + (1 | sp) + sub_dat$X + (0 + sub_dat$X | site), 
                  data = sub_dat, REML=F)
      z.0 <- lmer(Y ~ 1 + (1 | sp) + sub_dat$X, data = sub_dat, REML=F)
      lme4_results <- c(logLik(z.1), logLik(z.0), 
                        pchisq(2 * (logLik(z.1) - logLik(z.0)), 
                               df=1, lower.tail=F)/2)
      output$trait[counter] <- itrait  
      output$sp_richness[counter] <- nlevels(sub_dat$sp)
      output$dep[counter] <- "freq"
      output$lmer.logLik[counter] <- lme4_results[1]
      output$lmer.logLik0[counter] <- lme4_results[2]
      output$lmer.Pr[counter] <- lme4_results[3]
    }
    
    if(binary == TRUE){
      # lmer
      z.1 <- glmer(presence ~ 1 + (1 | sp) + sub_dat$X + (0 + sub_dat$X | site), 
                   data = sub_dat, family = binomial)
      z.0 <- glmer(presence ~ 1 + (1 | sp) + sub_dat$X, 
                   data = sub_dat, family = binomial)
      lme4_results <- c(logLik(z.1), logLik(z.0), 
                        pchisq(2 * (logLik(z.1) - logLik(z.0)), 
                               df=1, lower.tail=F)/2)
      output$trait[counter] <- itrait  
      output$sp_richness[counter] <- nlevels(sub_dat$sp)
      output$dep[counter] <- "presence"
      output$lmer.logLik[counter] <- lme4_results[1]
      output$lmer.logLik0[counter] <- lme4_results[2]
      output$lmer.Pr[counter] <- lme4_results[3]
    }
    
    counter <- counter + 1
    # show(c(output$p_trait[counter],output$s2_trait[counter],output$s2_attract[counter],output$s2_attract0[counter]))
  }
  output
}

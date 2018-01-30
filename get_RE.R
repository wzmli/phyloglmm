get_RE <- function(veg.long, phylo = pb.phylo, trait = pb.trait, 
                   trait.re = NULL, stand = TRUE,
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
  
  if(stand){
 Vphy <- Vphy/max(Vphy)                    ## Don't need?
 Vphy <- Vphy/det(Vphy)^(1/nspp)           ## Don't need?
  }
  show(c(nspp, Ntip(phy)))
  if(nspp != Ntip(phy)){
    stop("The vegetation data and the phylogeny have different number of species")
  }
  
  # random effect for site
  re.site <- list(1, site = dat$site, covar = diag(nsite))
  re.sp <- list(1, sp = dat$sp, covar = diag(nspp))
  re.sp.phy <- list(1, sp = dat$sp, covar = Vphy)
  
  re.nested.phy <- list(1, sp = dat$sp, covar = Vphy, site = dat$site)
  
  return(list(re.site,re.sp,re.sp.phy,re.nested.phy))
}

get_phylo <- function(veg.long, phylo = pb.phylo, trait = pb.trait, 
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
  return(phy)
  }

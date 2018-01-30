library(lme4)
library(pez)

str(dat)
reltol <- 0.000000001
maxit <- 100000000

## abundance data ----
start.model = lmer(Y ~ 1 + (1|sp) + (1|site), data = dat, REML = FALSE)
SigAIC(start.model)
AIC(start.model)

# search for fixed terms ------
i = 0
modellist = list()
(Aic1 =  AIC(start.model))
Aic0 = Aic1 + 1
block <- as.vector(dune.traits2) ## hack?
block <- names(dune.traits2)[-1]
fix.termsF = NULL
AicF = AIC(start.model)
while (Aic1 < Aic0){
  i = i + 1
  Aic0 = Aic1 # the start model's AIC to be Aic0
  cat("Initial Aic", Aic1, "\n")
  Model.sel <- Model_select_tier2F(start.model, block)
  selected.model <- Model.sel$new.model
  Aic1 = AIC(selected.model) # the new model's AIC to be Aic1
  cat("After Aic", Aic1, "\n")
  block <- Model.sel$B # remove the tested trait from the trait list
  terms <- Model.sel$terms # the tested trait
  # save model and terms
  if (Aic1 < Aic0){ # if the new model has lower AIC
    start.model = selected.model # make the new model as the base model
    final.modelF = selected.model # and possible final model
    fix.termsF = c(fix.termsF, Model.sel$term) # update fixed terms to be included
  } else { # if the new model has higher AIC
    final.modelF = start.model # we get the final model
  }
  modellist[[i]] = selected.model
  AicF = c(AicF, Aic1)
  print(c(Model.sel$term, Aic1))
}
final.modelF

start.model = final.modelF

## search for random terms  -------
i = 0
modellist = list()
(Aic1 =  AIC(start.model))
Aic0 = Aic1 + 1
block <- as.vector(dune.traits2) ## hack?
block <- names(dune.traits2)[c(2,6)]
random.terms = NULL
AicR = AIC(start.model)

while (Aic1 < Aic0){
  i = i + 1
  Aic0 = Aic1 # Aic of start model
  Model.sel <- Model_select_tier1R(start.model, block)
  start.model <- Model.sel$new.model # new model to be the new start model
  Aic1 = AIC(start.model) # AIC of new model
  block <- Model.sel$B # update trait list by removing the tested one
  terms <- Model.sel$terms 
  # save model and terms
  if (Aic1 < Aic0){ # if the new model has lower AIC
    final.modelR = start.model # the new model will be a possible final model
    random.terms = c(random.terms, Model.sel$term) # all random terms
  }
  modellist[[i]] = start.model # save new model
  AicR = c(AicR, Aic1)
  print(c(Model.sel$term, Aic1))
}
# with log.sla as random effect
# with log.sla and annual as fixed terms

lmer(Y ~ 1 + log.sla + annual + (1|sp) + (1|site) + (0+log.sla|site), data = dat, REML = F)


re.sla = list(unname(unlist(dat["log.sla"])), site = dat$site, covar = diag(nsite))

REs <- get_RE(veg.long = dune.veg2
  , trait = dune.traits2[c(1, 2)]
  , trait.re = c("log.sla")
  , phylo = dune.phylo2
  , trans = "log")

## get_RE returns a list of random effects in this order:
## re.site
## re.sp
## re.phy
## re.nested.phy

re.site <- REs[[1]]
re.sp <- REs[[2]]
re.sp.phy <- REs[[3]]
re.nested.phy <- REs[[4]]

z_traitsRE = communityPGLMM(formula = "Y ~ 1 + log.sla + annual", 
                            data = dat, family = "gaussian", 
                            sp = dat$sp, site = dat$site, 
                            random.effects = list(re.sp, re.sp.phy, re.site, 
                                                  re.sla, re.nested.phy), 
                            REML = F, verbose = F, 
                            s2.init = c(1.5, rep(0.01, 4)), 
                            reltol = reltol, maxit = maxit)
z_traitsRE$AIC
z_traitsRE$s2r %>% round(6)
z_traitsRE$s2resid
z_traitsRE$s2n
z_traitsRE$B
z_traitsRE$B.pvalue


## table 2A, without traits as Random effects
z_no_traitsRE = communityPGLMM(formula = "Y ~ 1 + log.sla + annual",
                               data = dat, family = "gaussian", 
                               sp = dat$sp, site = dat$site, 
                               random.effects = list(re.sp, re.sp.phy, re.site, 
                                                     re.nested.phy), 
                               REML = F, verbose = F, 
                               s2.init = c(1.5, rep(0.01, 3)), 
                               reltol = reltol, maxit = maxit)

z_no_traitsRE$s2r %>% round(6)
z_no_traitsRE$s2resid
z_no_traitsRE$s2n
z_no_traitsRE$B
z_no_traitsRE$B.pvalue

(z_no_traitsRE$s2n - z_traitsRE$s2n)/z_no_traitsRE$s2n # 0.1909708

# any additional traits that can further reduce the s2n? ----
sel_1_pglmm = selection(veg = dune.veg2, trait = dune.traits2, phylo = dune.phylo2, 
                        binary = F, added.traits = c("log.sla"),
                        fixed.terms = c("log.sla", "annual"))
#                   traits s2_nested_w_t s2_nested_wo_t aic_w_t aic_wo_t      prop
# 1 log.sla+log.seed.mass   0.006199256    0.008292037      NA       NA 0.2523844
# 2      log.sla+log.ldmc   0.006612347    0.008231912      NA       NA 0.1967422
# 3        log.sla+annual   0.006643365    0.008211569      NA       NA 0.1909750
# 4    log.sla+log.height   0.006692554    0.008179769      NA       NA 0.1818163

sel_2_pglmm = selection(veg = dune.veg2, trait = dune.traits2, phylo = dune.phylo2, 
                        binary = F, added.traits = c("log.sla", "log.seed.mass"),
                        fixed.terms = c("log.sla", "annual", "log.seed.mass"))
# as log.seed.mass is selected as random term, we should put it in the fixed term also
#                             traits s2_nested_w_t s2_nested_wo_t aic_w_t aic_wo_t      prop
# 1     log.sla+log.seed.mass+annual   0.006199340    0.008292037      NA       NA 0.2523743
# 2 log.sla+log.seed.mass+log.height   0.006333981    0.008255810      NA       NA 0.2327851
# 3   log.sla+log.seed.mass+log.ldmc   0.006426451    0.008293603      NA       NA 0.2251316

## the final selected traits will be ----
# log.sla, annual, log.seed.mass for fixed terms
# log.sla, log.seed.mass for random terms
# so that the phylo-sginal decreased the most (~25%)

z_final_traitsRE = phylo_explained_by_multi_traits_re_sel(veg.long = dune.veg2,
                                                          trait = dune.traits2[c(1, 2, 5, 6)], 
                                                          trait.re = c("log.sla", "log.seed.mass"),
                                                          phylo = dune.phylo2, 
                                                          trans = "log")
## extract details to fill Table 2B ----
# rds files are removed as they are large.
z = readRDS("rds/trait_multi_reg_log_full.rds")
z$s2n # sigma^2_c,  0.006199256
z$s2resid # 0.3332155
round(z$s2r, 6)
# intercept, intercept_phy, sla, seed mass, site
# 0.036447 0.000364 0.015503 0.013434 0.003638
# coef of fixed terms
z$formula
z$B
z$B.pvalue
z$AIC

z0 = readRDS("rds/trait_multi_reg_log_no_traits.rds")
z0$AIC
z0$s2n
round(z0$s2r, 6)
z0$s2resid
z0$B
z0$B.pvalue

(z0$s2n - z$s2n)/z0$s2n # 0.2523825

## presence/absence data ----
dat$presence = as.numeric(dat$freq > 0)
start.model = glmer(presence ~ 1 + (1|sp) + (1|site), data = dat, family = binomial())
SigAIC(start.model)
AIC(start.model)

i = 0
modellist = list()    
(SigAIC1 =  AIC(start.model)) 
SigAIC0 = SigAIC1 + 1
block = as.vector(traits)
fix.terms = NULL
Sig.AICR = AIC(start.model)

while (SigAIC1 < SigAIC0){
  i = i + 1
  SigAIC0 = SigAIC1
  cat("Initial sigAIC", SigAIC1, "\n")
  Model.sel <- Model_select_tier2F(start.model, block, sig.aic = F) 
  selected.model <- Model.sel$new.model
  SigAIC1 = AIC(selected.model)
  cat("After sigAIC", SigAIC1, "\n")
  block <- Model.sel$B
  terms <- Model.sel$terms
  # save model and terms
  if (SigAIC1 < SigAIC0){
    start.model = selected.model
    final.modelF = selected.model
    fix.termsF = c(fix.terms, Model.sel$term)
  } else {
    final.modelF = start.model
  }
  modellist[[i]] = selected.model
  Sig.AICR = c(Sig.AICR, SigAIC1)
  print(c(Model.sel$term, SigAIC1))
} 


final.modelF
start.model = final.modelF

i = 0
modellist = list()    
Sig.AIC = NULL
(SigAIC1 =  AIC(start.model)) 
SigAIC0 = SigAIC1 + 1
block = as.vector(traits)
random.terms = NULL
Sig.AICR = AIC(start.model)

while (SigAIC1 < SigAIC0){
  i = i + 1
  SigAIC0 = SigAIC1
  Model.sel <- Model_select_tier1R(start.model, block, sig.aic = F) 
  start.model <- Model.sel$new.model
  SigAIC1 = AIC(start.model)
  block <- Model.sel$B
  terms <- Model.sel$terms
  # save model and terms
  if (SigAIC1 < SigAIC0){
    final.modelR = start.model
    random.termsF = c(random.terms, Model.sel$term)
  }
  modellist[[i]] = start.model
  random.terms = c(random.terms, Model.sel$term)
  Sig.AICR = c(Sig.AICR, SigAIC1)
  print(c(Model.sel$term, SigAIC1))
} 


# how much do there traits reduce phylo variation?
dune_multi_traits_binary_selected = phylo_explained_by_multi_traits_re_sel(veg.long = dune.veg2, 
                                                                           trait = dune.traits2[c(1, 2, 6)], 
                                                                           trait.re = "log.sla", 
                                                                           phylo = dune.phylo2, 
                                                                           binary = T)
# s2_attract_with_traits    1.089702e-02
# p_attract_with_traits     3.595031e-01
# s2_attract_without_traits 3.963401e-02
# p_attract_without_traits  8.755730e-02
# p_traits                  2.269264e-05
# s2_attract_desc           7.250589e-01

# keep reaching
sel_1_binary = selection(veg = dune.veg2, trait = dune.traits2, phylo = dune.phylo2, 
                         binary = T, added.traits = c("log.sla", "annual"))
# traits s2_nested_w_t s2_nested_wo_t aic_w_t aic_wo_t      prop
# 1 log.sla+annual+log.seed.mass   0.002175528     0.04011705      NA       NA 0.9457705
# 2    log.sla+annual+log.height   0.009915030     0.03962318      NA       NA 0.7497669
# 3      log.sla+annual+log.ldmc   0.010414697     0.03897093      NA       NA 0.7327573
sel_1_binary2 = selection(veg = dune.veg2, trait = dune.traits2, phylo = dune.phylo2, 
                          binary = T, added.traits = c("log.sla", "annual", "log.seed.mass"))
# traits s2_nested_w_t s2_nested_wo_t aic_w_t aic_wo_t      prop
# 1 log.sla+annual+log.seed.mass+log.height  0.0008283252     0.03979028      NA       NA 0.9791827
# 2   log.sla+annual+log.seed.mass+log.ldmc  0.0021192854     0.03930879      NA       NA 0.9460862

sel_1_binary3 = selection(veg = dune.veg2, trait = dune.traits2, phylo = dune.phylo2, 
                          binary = T, added.traits = c("log.sla", "annual", "log.seed.mass", "log.height"))
# traits s2_nested_w_t s2_nested_wo_t aic_w_t aic_wo_t      prop
# 1 log.sla+annual+log.seed.mass+log.height+log.ldmc   0.001198059  0.03927794   NA    NA 0.9694979


dune_multi_traits_binary2 = phylo_explained_by_multi_traits(veg.long = dune.veg2, 
                                                            trait = dune.traits2[, -4], 
                                                            phylo = dune.phylo2, 
                                                            binary = T)
# s2_attract_with_traits    0.0008713322
# p_attract_with_traits     0.4886422099
# s2_attract_without_traits 0.0398037571
# p_attract_without_traits  0.0874398906
# p_traits                  0.0002811005
# s2_attract_desc           0.9781092978

## extract details to fill Table S2 ----
# rds files are removed as they are large.
# z.binary = readRDS("rds/binary_trait_multi_reg_full.rds")
z.binary$s2n
round(z.binary$s2r, 6)
z.binary

# z0.binary = readRDS("rds/binary_trait_multi_reg_no_traits.rds")
z0.binary$s2n
z0.binary$s2r %>% round(6)
z0.binary
(z0.binary$s2n - z.binary$s2n)/z0.binary$s2n

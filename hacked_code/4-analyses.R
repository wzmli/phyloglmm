
##################### analyses ---------------------
#### Patterns (Table 1) ----
dune_pattern = phylo_pattern(veg.long = dune.veg2, 
                              phylo = dune.phylo2, 
                              trans = "log")
dune_pattern
#       models sp_num   sigma_sp sigma_sp_phy sigma_nested  sigma_site sigma_resid    p_nested      AIC
# 1 attraction     28 0.08008604 0.0009584519 7.805741e-03 0.002654685   0.3565561 0.005321897 1122.177
# 2  repulsion     28 0.09006573 0.0012628749 2.644303e-13 0.002905533   0.3908993 0.500000000 1128.701
# 3  no_nested     28 0.09006584 0.0029055416           NA 0.001262860   0.3908993          NA 1126.701

dune_pattern_binary = phylo_pattern(veg.long = dune.veg2, binary = TRUE,
                                 phylo = dune.phylo2)
dune_pattern_binary
#       models spp_num  sigma_sp sigma_sp_phy sigma_nested  sigma_site  p_nested
# 1 attraction      28 0.9074453 1.966551e-11 3.213948e-02 0.005790093 0.1465297
# 2  repulsion      28 0.9704567 2.126978e-11 2.034572e-12 0.017977236 0.5000000
# 3  no_nested      28 0.9704636 5.041391e-12           NA 0.017977424        NA

#### Select traits to check % of phylogenetic variation they can reduce ----
# Table 2
# see the 2-forward_selection_fixed_first.R

#### phylogentic signal of traits (Table 3) ----
# see 3-phylosig.R

#### environmental variables (Table 4) ----
dune_envi_log = phylo_signal_slopes_envi(dat = select(dune.all, -(log.sla:annual)) %>% 
                                  mutate(Management = as.numeric(Management)), 
                                phylo = dune.phylo2, binary = FALSE)
dune_envi_log
#          var slope.var.lmer.Pr slope.var.indep.pglmm.Pr slope.var.phylo.pglmm.Pr slope.estimated.pglmm slope.estimated.pglmm.Pr
# 1     log.A1           0.00075                  0.00075                  0.50000           -0.03366327                  0.37422
# 2   Moisture           0.00000                  0.00000                  0.50000           -0.02085463                  0.47082
# 3 Management           0.22106                  0.27570                  0.18610           -0.05816512                  0.05716
# 4        Use           0.09931                  0.10810                  0.34976           -0.01072823                  0.79869
# 5     Manure           0.00000                  0.00000                  0.02511            0.01578035                  0.72443

dune_envi_binary = phylo_signal_slopes_envi(dat = select(dune.all, -(log.sla:annual)) %>% 
                                  mutate(Management = as.numeric(Management)), 
                                phylo = dune.phylo2, 
                                binary = T)
dune_envi_binary
#          var slope.var.lmer.Pr slope.var.indep.pglmm.Pr slope.var.phylo.pglmm.Pr slope.estimated.pglmm slope.estimated.pglmm.Pr
# 1     log.A1           0.00208                  0.00494                  0.50000          -0.094440018                  0.49830
# 2   Moisture           0.00000                  0.00000                  0.50000          -0.105790450                  0.30745
# 3 Management           0.48536                  0.50000                  0.34757          -0.195326895                  0.05114
# 4        Use           0.25296                  0.19076                  0.50000          -0.084944856                  0.53017
# 5     Manure           0.00002                  0.00013                  0.05940          -0.005358705                  0.97140


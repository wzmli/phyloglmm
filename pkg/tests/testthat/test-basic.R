library("ape")
library("testthat")
library("phyloglmm")

phyloZ <- phylo.to.Z(garamszegi_phy)
datG <- transform(garamszegi_simple,
                  obs = factor(seq(nrow(garamszegi_simple))),
                  sp = phylo)

test_that("basic phylo.to.Z", {
  expect_s4_class(phyloZ, "Matrix")
  expect_equal(nrow(phyloZ), ape::Ntip(garamszegi_phy))
  expect_equal(ncol(phyloZ), ape::Ntip(garamszegi_phy) + ape::Nnode(garamszegi_phy) - 1)
})

test_that("basic lme4 LMM", {
  phylo_lmm_fit <- phylo_lmm(phen~cofactor+(1|sp)
                           , data=datG
                           , phylonm = "sp"
                           , phylo = phylo
                           , phyloZ = phyloZ
                           , REML = TRUE
                           , control = lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
                             )
  expect_equal(c(lme4::VarCorr(phylo_lmm_fit)$sp), 207.02517)
})

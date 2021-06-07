library("ape")
library("testthat")
library("phyloglmm")

phyloZ <- phylo.to.Z(garamszegi_phy)

augment_data <- function(x) {
  transform(x,
            obs = factor(seq(nrow(x))),
            sp = factor(phylo)
            )
}

simp_VC <- function(x) {
  v <- lme4::VarCorr(x)
  if (inherits(x, "glmmTMB")) v <- v$cond
  unlist(lapply(v, c))
}

datG <- augment_data(garamszegi_simple)
datP <- augment_data(garamszegi_pois)

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


test_that("phylo_glmm", {
  phylo_glmm_fit <- phylo_glmm(phen_pois~cofactor+(1|sp)+(1|obs)
                             , data=datP
                             , phylonm = "sp"
                             , family = poisson
                             , phylo = phylo
                             , phyloZ=phyloZ
                             , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                               )
  phylo_glmmTMB_fit <- phylo_glmmTMB(phen_pois~cofactor+(1|sp)+(1|obs)
                             , data=datP
                             , phylonm = "sp"
                             , family = poisson
                             , phyloZ=phyloZ
                               )
## Warning message:
## In condReStruc[i] <- `*vtmp*` :
##   number of items to replace is not a multiple of replacement length

  expect_equal(simp_VC(phylo_glmm_fit), simp_VC(phylo_glmmTMB_fit),
               tolerance = 1e-3)
})

sf <- function(x) system.file("test_data", x, package="phyloglmm", mustWork=TRUE)

dd <- read.table(sf("GH6_data.txt"), header=TRUE)
pp <- ape::read.nexus(sf("GH6_tree.nex"))
phyloZ <- phylo.to.Z(pp)
cc <- lmerControl(check.nobs.vs.nlev="ignore",
                  check.nobs.vs.nRE="ignore")

test_that("check tips vs data species", {
  expect_error(phylo_lmm_fit <- phylo_lmm(v~bio1+(1|sp)
                         , data=dd,
                         , phylonm = "sp"
                         , phyloZ = phyloZ
                         , control = cc),
               "in data but not phyloZ")
})

library("ape")
library("testthat")
library("phyloglmm")

cc <- lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
phyloZ0 <- phylo.to.Z(garamszegi_phy)

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

phyloZ <- phyloZ0[levels(datG$sp), ]

test_that("basic phylo.to.Z", {
  expect_s4_class(phyloZ, "Matrix")
  expect_equal(nrow(phyloZ), ape::Ntip(garamszegi_phy))
  expect_equal(ncol(phyloZ), ape::Ntip(garamszegi_phy) + ape::Nnode(garamszegi_phy) - 1)
})

phylo_lmm_fit <- phylo_lmm(phen~cofactor+(1|sp)
                         , data=datG
                         , phylonm = "sp"
                         , phylo = garamszegi_phy
                         , REML = TRUE
                         , control = cc
                           )

test_that("basic lme4 LMM", {
  expect_equal(c(lme4::VarCorr(phylo_lmm_fit)$sp), 207.02517)
})

test_that("fit with external phyloZ", {
  phylo_lmm_fitZ <- phylo_lmm(phen~cofactor+(1|sp)
                           , data=datG
                           , phylonm = "sp"
                           , phyloZ = phyloZ
                           , REML = TRUE
                           , control = cc
                             )
  ## too tedious to check identity with all these environments & calls
  ## floating around ...
  expect_identical(simp_VC(phylo_lmm_fit), simp_VC(phylo_lmm_fitZ))
})

test_that("detect phyloZ/sp mismatch", {
  expect_error(phylo_lmm(phen~cofactor+(1|sp)
          , data=datG
          , phylonm = "sp"
          , phyloZ = phyloZ0
          , REML = TRUE
          , control = cc
            ),
            "row names of phyloZ must match")
})

test_that("require sp to be a factor", {
  datG$sp2 <- as.character(datG)
  expect_error(phylo_lmm(phen~cofactor+(1|sp)
          , data=datG
          , phylonm = "sp2"
          , phyloZ = phyloZ0
          , REML = TRUE
          , control = cc
            ),
            "species identifier in data must be a factor")
})


test_that("phylo_glmm", {
  phylo_glmm_fit <- phylo_glmm(phen_pois~cofactor+(1|sp)+(1|obs)
                             , data=datP
                             , phylonm = "sp"
                             , family = poisson
                             , phylo = garamszegi_phy
                             , control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")
                               )

  phylo_glmmTMB_fit <- phylo_glmmTMB(phen_pois~cofactor+(1|sp)+(1|obs)
                             , data=datP
                             , phylonm = "sp"
                             , family = poisson
                             , phylo = garamszegi_phy
                               )

  expect_equal(simp_VC(phylo_glmm_fit), simp_VC(phylo_glmmTMB_fit),
               tolerance = 1e-3)
})

set.seed(101)
tt <- ape::rtree(20)
dd <- data.frame(y = rnorm(20), x = rnorm(20), sp = paste0("t", 1:20))
dd[20, "sp"] <- "T20"
dd$sp <- factor(dd$sp)
pZ <- phylo.to.Z(tt)
cc <- lmerControl(check.nobs.vs.nlev="ignore",
                  check.nobs.vs.nRE="ignore")

test_that("check tips vs data species", {
  expect_error(phylo_lmm(y ~ x + (1|sp),
                         data = dd,
                         phylonm = "sp",
                         phyloZ = pZ,
                         REML = FALSE,
                         control = cc),
               "in phyloZ but not data:.*in data but not phyloZ:")
})

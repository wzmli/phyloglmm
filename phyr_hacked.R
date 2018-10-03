prep_dat_pglmm_hacked <- function (formula, data, tree, repulsion = FALSE, prep.re.effects = TRUE, 
          family = "gaussian", prep.s2.lme4 = FALSE, tree_site = NULL, 
          bayes = FALSE) 
{
  if (!all(c("sp", "site") %in% names(data))) {
    stop("The data frame should have a column named as 'sp' and a column named as 'site'.")
  }
  data = dplyr::arrange(as.data.frame(data), site, sp)
  data$sp = as.factor(data$sp)
  sp = data$sp
  data$site = as.factor(data$site)
  site = data$site
  spl = levels(sp)
  sitel = levels(site)
  nspp = nlevels(sp)
  nsite = nlevels(site)
  fm = unique(lme4::findbars(formula))
  if (prep.re.effects) {
    if (is.null(fm)) 
      stop("No random terms specified, use lm or glm instead")
    if (any(grepl("sp__", fm))) {
      if (class(tree) == "phylo") {
        if (length(setdiff(spl, tree$tip.label))) 
          stop("Some species not in the phylogeny, please either drop these species or update the phylogeny")
        if (length(setdiff(tree$tip.label, spl))) {
          warning("Drop species from the phylogeny that are not in the data", 
                  immediate. = TRUE)
          tree = ape::drop.tip(tree, setdiff(tree$tip.label, 
                                             spl))
        }
        Vphy <- ape::vcv(tree)
#         Vphy <- Vphy/max(Vphy)
#         Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/nspp)
#         Vphy = Vphy[spl, spl]
      }
      if (inherits(tree, c("matrix", "Matrix"))) {
        if (length(setdiff(spl, row.names(tree)))) 
          stop("Some species not in the cov matrix, please either drop these species or update the matrix")
        if (length(setdiff(row.names(tree), spl))) {
          warning("Drop species from the cov matrix that are not in the data", 
                  immediate. = TRUE)
        }
        tree = tree[spl, spl]
        if ((det(tree) - 1) > 1e-04) {
          warning("The cov matrix is not standarized, we will do this now...", 
                  immediate. = TRUE)
          tree <- tree/max(tree)
          tree <- tree/exp(determinant(tree)$modulus[1]/nrow(tree))
          if ((det(tree) - 1) > 1e-04) 
            warning("Failed to standarized the cov matrix", 
                    immediate. = TRUE)
        }
        Vphy = tree
      }
    }
    else {
      Vphy = NULL
    }
    if (any(grepl("site__", fm))) {
      if (is.null(tree_site)) 
        stop("tree_site not specified")
      if (class(tree_site) == "phylo") {
        if (length(setdiff(sitel, tree_site$tip.label))) 
          stop("Some species not in the phylogeny tree_site, please either drop these species or update the phylogeny")
        if (length(setdiff(tree_site$tip.label, sitel))) {
          warning("Drop species from the phylogeny tree_site that are not in the data", 
                  immediate. = TRUE)
          tree = ape::drop.tip(tree_site, setdiff(tree_site$tip.label, 
                                                  sitel))
        }
        Vphy_site <- ape::vcv(tree_site)
        Vphy_site <- Vphy_site/max(Vphy_site)
        Vphy_site <- Vphy_site/exp(determinant(Vphy_site)$modulus[1]/nsite)
        Vphy_site = Vphy_site[sitel, sitel]
      }
      if (inherits(tree_site, c("matrix", "Matrix"))) {
        if (length(setdiff(sitel, row.names(tree_site)))) 
          stop("Some species not in the cov matrix tree_site, please either drop these species or update the matrix tree_site")
        if (length(setdiff(row.names(tree_site), sitel))) {
          warning("Drop species from the cov matrix that are not in the data", 
                  immediate. = TRUE)
        }
        tree_site = tree_site[sitel, sitel]
        if ((det(tree_site) - 1) > 1e-04) {
          warning("The cov matrix is not standarized, we will do this now...", 
                  immediate. = TRUE)
          tree_site <- tree_site/max(tree_site)
          tree_site <- tree_site/exp(determinant(tree_site)$modulus[1]/nrow(tree_site))
          if ((det(tree_site) - 1) > 1e-04) 
            warning("Failed to standarized the cov matrix", 
                    immediate. = TRUE)
        }
        Vphy_site = tree_site
      }
    }
    else {
      Vphy_site = NULL
    }
    if (nrow(data) != nspp * nsite) {
      message("the dataframe may have been removed for NAs as its number of row is not nspp * nsite \n\n              we will recreate the full data frame for you.")
      data_all = dplyr::arrange(expand.grid(site = sitel, 
                                            sp = spl), site, sp)
      data_all = dplyr::left_join(data_all, data, by = c("site", 
                                                         "sp"))
      nna.ind = which(!is.na(data_all[, as.character(formula)[2]]))
      if (nrow(data) != length(nna.ind)) 
        stop("something wrong with NAs")
    }
    if (all(grepl("@", fm) == FALSE)) {
      n_repulsion = 1
    }
    else {
      n_repulsion = sum(sapply(fm[grepl("@", fm)], function(x) {
        xx = strsplit(as.character(x)[3], "@")[[1]]
        sum(grepl("__", xx))
      }))
    }
    if (length(repulsion) == 1) 
      repulsion = rep(repulsion, n_repulsion)
    if (length(repulsion) != n_repulsion) 
      stop("the number of repulsion terms specified is not right, please double check")
    nested_repul_i = 1
    random.effects = lapply(fm, function(x) {
      x2 = as.character(x)
      x2 = gsub(pattern = "^0 ?[+] ?", replacement = "", 
                x2)
      if (grepl("[+]", x2[2])) 
        stop("(x1 + x2|g) form of random terms are not allowed yet, pleast split it")
      if (x2[2] == "1") {
        if (!grepl("[@]", x2[3])) {
          if (grepl("__$", x2[3])) {
            coln = gsub("__$", "", x2[3])
            if (coln %nin% c("sp", "site")) 
              stop("group variable with phylogenetic var-covar matrix must be named as either sp or site")
            d = data[, coln]
            xout_nonphy = list(1, d, covar = diag(nlevels(d)))
            names(xout_nonphy)[2] = coln
            if (coln == "sp") {
              xout_phy = list(1, d, covar = Vphy)
            }
            if (coln == "site") {
              xout_phy = list(1, d, covar = Vphy_site)
            }
            names(xout_phy)[2] = coln
            xout = list(xout_nonphy, xout_phy)
          }
          else {
            d = data[, x2[3]]
            xout = list(1, d, covar = diag(length(unique(d))))
            names(xout)[2] = x2[3]
            xout = list(xout)
          }
        }
        else {
          sp_or_site = strsplit(x2[3], split = "@")[[1]]
          if (!grepl("__", x2[3])) {
            message("Nested term without specify phylogeny, use identity matrix instead")
            xout = list(1, sp = data[, sp_or_site[1]], 
                        covar = diag(length(unique(data[, sp_or_site[1]]))), 
                        site = data[, sp_or_site[2]])
            names(xout)[c(2, 4)] = sp_or_site
            xout = list(xout)
          }
          else {
            if (sp_or_site[1] == "sp__" & !grepl("__", 
                                                 sp_or_site[2])) {
              if (bayes) {
                if (repulsion[nested_repul_i]) {
                  xout = list(1, sp, covar = solve(Vphy), 
                              data[, sp_or_site[2]])
                }
                else {
                  xout = list(1, sp, covar = Vphy, data[, 
                                                        sp_or_site[2]])
                }
              }
              else {
                n_dim = length(unique(data[, sp_or_site[2]]))
                if (repulsion[nested_repul_i]) {
                  xout = as(kronecker(diag(n_dim), solve(Vphy)), 
                            "dgCMatrix")
                }
                else {
                  xout = as(kronecker(diag(n_dim), Vphy), 
                            "dgCMatrix")
                }
                xout = list(xout)
              }
              nested_repul_i <<- nested_repul_i + 1
            }
            if (sp_or_site[1] == "sp" & sp_or_site[2] == 
                "site__") {
              if (repulsion[nested_repul_i]) {
                xout = as(kronecker(solve(Vphy_site), 
                                    diag(nspp)), "dgCMatrix")
              }
              else {
                xout = as(kronecker(Vphy_site, diag(nspp)), 
                          "dgCMatrix")
              }
              xout = list(xout)
              nested_repul_i <<- nested_repul_i + 1
            }
            if (sp_or_site[1] == "sp__" & sp_or_site[2] == 
                "site__") {
              if (repulsion[nested_repul_i]) {
                Vphy2 = solve(Vphy)
              }
              else {
                Vphy2 = Vphy
              }
              nested_repul_i <<- nested_repul_i + 1
              if (repulsion[nested_repul_i]) {
                Vphy_site2 = solve(Vphy_site)
              }
              else {
                Vphy_site2 = Vphy_site
              }
              nested_repul_i <<- nested_repul_i + 1
              xout = as(kronecker(Vphy_site2, Vphy2), 
                        "dgCMatrix")
              xout = list(xout)
            }
            if (nrow(data) != nspp * nsite) 
              xout[[1]] = xout[[1]][nna.ind, nna.ind]
            xout = list(xout)
          }
        }
      }
      else {
        if (grepl("@", x2[3])) 
          stop("sorry, random terms for slopes cannot be nested")
        if (grepl("__$", x2[3])) {
          coln = gsub("__$", "", x2[3])
          if (coln %nin% c("sp", "site")) 
            stop("group variable with phylogenetic var-covar matrix must be named as either sp or site")
          d = data[, coln]
          xout_nonphy = list(data[, x2[2]], d, covar = diag(nlevels(d)))
          names(xout_nonphy)[2] = coln
          xout_phy = list(data[, x2[2]], d, covar = Vphy)
          names(xout_phy)[2] = coln
          xout = list(xout_nonphy, xout_phy)
        }
        else {
          d = data[, x2[3]]
          xout = list(data[, x2[2]], d, covar = diag(dplyr::n_distinct(d)))
          names(xout)[2] = x2[3]
          xout = list(xout)
        }
      }
      xout
    })
    random.effects = unlist(random.effects, recursive = FALSE)
    names(random.effects) = unlist(sapply(fm, function(x) {
      x2 = as.character(x)
      x3 = paste0(x2[2], x2[1], x2[3])
      if (grepl("__$", x2[3]) & !grepl("@", x2[3])) {
        x4 = gsub("__$", "", x3)
        return(c(x4, x3))
      }
      x3
    }), recursive = T)
    if (prep.s2.lme4) {
      s2_init = numeric(length(random.effects))
      names(s2_init) = names(random.effects)
      dv = sqrt(var(lm(formula = lme4::nobars(formula), 
                       data = data)$residuals)/length(s2_init))
      s2_init[grep(pattern = "__", x = names(s2_init))] = dv
      s2_init[grep(pattern = "@", x = names(s2_init))] = dv
      no__ = gsub(pattern = "__", replacement = "", x = formula)
      if (grepl(pattern = "@", no__[3])) {
        no__[3] = gsub(pattern = "[+] *[(]1 *[|] *[sitep]{2,4}@[sitep]{2,4}[)]", 
                       replacement = "", x = no__[3])
      }
      fm_no__ = as.formula(paste0(no__[2], no__[1], no__[3]))
      if (family == "gaussian") {
        itheta = lme4::lFormula(fm_no__, data)$reTrms$theta
      }
      else {
        itheta = lme4::glFormula(fm_no__, data)$reTrms$theta
      }
      lt = grep("__|@", names(s2_init), invert = TRUE)
      if (length(lt) != length(itheta)) {
        warning("automated s2.init failed, set to NULL", 
                immediate. = TRUE)
        s2_init = NULL
      }
      else {
        s2_init[lt] = itheta
      }
    }
    else {
      s2_init = NULL
    }
  }
  else {
    random.effects = NA
    s2_init = NULL
  }
  formula = lme4::nobars(formula)
  return(list(formula = formula, data = data, sp = sp, site = site, 
              random.effects = random.effects, s2_init = s2_init, tree = tree, 
              tree_site = tree_site, Vphy = Vphy, Vphy_site = Vphy_site))
}

phyr::communityPGLMM
function (formula, data = NULL, family = "gaussian", tree = NULL, 
          tree_site = NULL, repulsion = FALSE, sp, site, random.effects = NULL, 
          REML = TRUE, bayes = FALSE, s2.init = NULL, B.init = NULL, 
          reltol = 10^-6, maxit = 500, tol.pql = 10^-6, maxit.pql = 200, 
          verbose = FALSE, ML.init = TRUE, marginal.summ = "mean", 
          calc.DIC = FALSE, default.prior = "inla.default", cpp = TRUE, 
          optimizer = c("nelder-mead-nlopt", "bobyqa", "Nelder-Mead", 
                        "subplex"), prep.s2.lme4 = FALSE) 
{
  optimizer = match.arg(optimizer)
  if ((family %nin% c("gaussian", "binomial")) & (bayes == 
                                                  FALSE)) {
    stop("\nSorry, but only binomial (binary) and gaussian options are available for\n         communityPGLMM at this time")
  }
  if (bayes) {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("To run communityPGLMM with bayes = TRUE, you need to install the packages 'INLA'.  \n           Please run in your R terminal:\n           install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
    }
    if ((family %nin% c("gaussian", "binomial", "poisson"))) {
      stop("\nSorry, but only binomial (binary), poisson (count), and gaussian options \n           are available for Bayesian communityPGLMM at this time")
    }
  }
  fm_original = formula
  prep_re = if (is.null(random.effects)) 
    TRUE
  else FALSE
  if (prep_re) {
    dat_prepared = prep_dat_pglmm(formula, data, tree, repulsion, 
                                  prep_re, family, prep.s2.lme4, tree_site, bayes)
    formula = dat_prepared$formula
    data = dat_prepared$data
    sp = dat_prepared$sp
    site = dat_prepared$site
    random.effects = dat_prepared$random.effects
    tree = dat_prepared$tree
    tree_site = dat_prepared$tree_site
  }
  else {
    formula = lme4::nobars(formula)
    if (missing(sp)) 
      sp = as.factor(data$sp)
    if (missing(site)) 
      site = as.factor(data$site)
  }
  if (prep.s2.lme4) 
    s2.init = dat_prepared$s2_init
  if (bayes & ML.init & (family %in% c("binomial", "gaussian"))) {
    if (family == "gaussian") {
      ML.init.z <- phyr:::communityPGLMM.gaussian(formula = formula, 
                                           data = data, sp = sp, site = site, random.effects = random.effects, 
                                           REML = REML, s2.init = s2.init, B.init = B.init, 
                                           reltol = reltol, maxit = maxit, verbose = verbose, 
                                           cpp = cpp, optimizer = optimizer)
      s2.init <- c(ML.init.z$s2r, ML.init.z$s2n, ML.init.z$s2resid)
      B.init <- ML.init.z$B[, 1, drop = TRUE]
    }
    if (family == "binomial") {
      if (is.null(s2.init)) 
        s2.init <- 0.25
      ML.init.z <- communityPGLMM.binary(formula = formula, 
                                         data = data, sp = sp, site = site, random.effects = random.effects, 
                                         REML = REML, s2.init = s2.init, B.init = B.init, 
                                         reltol = reltol, maxit = maxit, tol.pql = tol.pql, 
                                         maxit.pql = maxit.pql, verbose = verbose, cpp = cpp, 
                                         optimizer = optimizer)
      s2.init <- c(ML.init.z$s2r, ML.init.z$s2n)
      B.init <- ML.init.z$B[, 1, drop = TRUE]
    }
  }
  if (bayes & ML.init & (family %nin% c("binomial", "gaussian"))) {
    warning("ML.init option is only available for binomial and gaussian families. You will have to \n            specify initial values manually if you think the default are problematic.")
  }
  if (bayes) {
    z <- communityPGLMM.bayes(formula = formula, data = data, 
                              family = family, sp = sp, site = site, random.effects = random.effects, 
                              s2.init = s2.init, B.init = B.init, verbose = verbose, 
                              REML = REML, marginal.summ = marginal.summ, calc.DIC = calc.DIC, 
                              default.prior = default.prior)
  }
  else {
    if (family == "gaussian") {
      z <- phyr:::communityPGLMM.gaussian(formula = formula, data = data, 
                                   sp = sp, site = site, random.effects = random.effects, 
                                   REML = REML, s2.init = s2.init, B.init = B.init, 
                                   reltol = reltol, maxit = maxit, verbose = verbose, 
                                   cpp = cpp, optimizer = optimizer)
    }
    if (family == "binomial") {
      if (is.null(s2.init)) 
        s2.init <- 0.25
      z <- communityPGLMM.binary(formula = formula, data = data, 
                                 sp = sp, site = site, random.effects = random.effects, 
                                 REML = REML, s2.init = s2.init, B.init = B.init, 
                                 reltol = reltol, maxit = maxit, tol.pql = tol.pql, 
                                 maxit.pql = maxit.pql, verbose = verbose, cpp = cpp, 
                                 optimizer = optimizer)
    }
  }
  z$formula_original = fm_original
  z$tree = tree
  z$tree_site = tree_site
  if (!is.null(names(random.effects))) {
    re.names = names(random.effects)[c(which(sapply(random.effects, 
                                                    length) %nin% c(1, 4)), which(sapply(random.effects, 
                                                                                         length) %in% c(1, 4)))]
    if (family == "gaussian") 
      re.names <- c(re.names, "residual")
    names(z$ss) = re.names
  }
  return(z)
}


phyr_hacked <- function (formula, data = NULL, family = "gaussian", tree = NULL, 
          tree_site = NULL, repulsion = FALSE, sp, site, random.effects = NULL, 
          REML = TRUE, bayes = FALSE, s2.init = NULL, B.init = NULL, 
          reltol = 10^-6, maxit = 500, tol.pql = 10^-6, maxit.pql = 200, 
          verbose = FALSE, ML.init = TRUE, marginal.summ = "mean", 
          calc.DIC = FALSE, default.prior = "inla.default", cpp = TRUE, 
          optimizer = c("nelder-mead-nlopt", "bobyqa", "Nelder-Mead", 
                        "subplex"), prep.s2.lme4 = FALSE) 
{
  optimizer = match.arg(optimizer)
  if ((family %nin% c("gaussian", "binomial")) & (bayes == 
                                                  FALSE)) {
    stop("\nSorry, but only binomial (binary) and gaussian options are available for\n         communityPGLMM at this time")
  }
  if (bayes) {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("To run communityPGLMM with bayes = TRUE, you need to install the packages 'INLA'.  \n           Please run in your R terminal:\n           install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
    }
    if ((family %nin% c("gaussian", "binomial", "poisson"))) {
      stop("\nSorry, but only binomial (binary), poisson (count), and gaussian options \n           are available for Bayesian communityPGLMM at this time")
    }
  }
  fm_original = formula
  prep_re = if (is.null(random.effects)) 
    TRUE
  else FALSE
  if (prep_re) {
    dat_prepared = prep_dat_pglmm_hacked(formula, data, tree, repulsion, 
                                  prep_re, family, prep.s2.lme4, tree_site, bayes)
    formula = dat_prepared$formula
    data = dat_prepared$data
    sp = dat_prepared$sp
    site = dat_prepared$site
    random.effects = dat_prepared$random.effects
    tree = dat_prepared$tree
    tree_site = dat_prepared$tree_site
  }
  else {
    formula = lme4::nobars(formula)
    if (missing(sp)) 
      sp = as.factor(data$sp)
    if (missing(site)) 
      site = as.factor(data$site)
  }
  if (prep.s2.lme4) 
    s2.init = dat_prepared$s2_init
  if (bayes & ML.init & (family %in% c("binomial", "gaussian"))) {
    if (family == "gaussian") {
      ML.init.z <- phyr:::communityPGLMM.gaussian(formula = formula, 
                                           data = data, sp = sp, site = site, random.effects = random.effects, 
                                           REML = REML, s2.init = s2.init, B.init = B.init, 
                                           reltol = reltol, maxit = maxit, verbose = verbose, 
                                           cpp = cpp, optimizer = optimizer)
      s2.init <- c(ML.init.z$s2r, ML.init.z$s2n, ML.init.z$s2resid)
      B.init <- ML.init.z$B[, 1, drop = TRUE]
    }
    if (family == "binomial") {
      if (is.null(s2.init)) 
        s2.init <- 0.25
      ML.init.z <- communityPGLMM.binary(formula = formula, 
                                         data = data, sp = sp, site = site, random.effects = random.effects, 
                                         REML = REML, s2.init = s2.init, B.init = B.init, 
                                         reltol = reltol, maxit = maxit, tol.pql = tol.pql, 
                                         maxit.pql = maxit.pql, verbose = verbose, cpp = cpp, 
                                         optimizer = optimizer)
      s2.init <- c(ML.init.z$s2r, ML.init.z$s2n)
      B.init <- ML.init.z$B[, 1, drop = TRUE]
    }
  }
  if (bayes & ML.init & (family %nin% c("binomial", "gaussian"))) {
    warning("ML.init option is only available for binomial and gaussian families. You will have to \n            specify initial values manually if you think the default are problematic.")
  }
  if (bayes) {
    z <- communityPGLMM.bayes(formula = formula, data = data, 
                              family = family, sp = sp, site = site, random.effects = random.effects, 
                              s2.init = s2.init, B.init = B.init, verbose = verbose, 
                              REML = REML, marginal.summ = marginal.summ, calc.DIC = calc.DIC, 
                              default.prior = default.prior)
  }
  else {
    if (family == "gaussian") {
      z <- phyr:::communityPGLMM.gaussian(formula = formula, data = data, 
                                   sp = sp, site = site, random.effects = random.effects, 
                                   REML = REML, s2.init = s2.init, B.init = B.init, 
                                   reltol = reltol, maxit = maxit, verbose = verbose, 
                                   cpp = cpp, optimizer = optimizer)
    }
    if (family == "binomial") {
      if (is.null(s2.init)) 
        s2.init <- 0.25
      z <- communityPGLMM.binary(formula = formula, data = data, 
                                 sp = sp, site = site, random.effects = random.effects, 
                                 REML = REML, s2.init = s2.init, B.init = B.init, 
                                 reltol = reltol, maxit = maxit, tol.pql = tol.pql, 
                                 maxit.pql = maxit.pql, verbose = verbose, cpp = cpp, 
                                 optimizer = optimizer)
    }
  }
  z$formula_original = fm_original
  z$tree = tree
  z$tree_site = tree_site
  if (!is.null(names(random.effects))) {
    re.names = names(random.effects)[c(which(sapply(random.effects, 
                                                    length) %nin% c(1, 4)), which(sapply(random.effects, 
                                                                                         length) %in% c(1, 4)))]
    if (family == "gaussian") 
      re.names <- c(re.names, "residual")
    names(z$ss) = re.names
  }
  return(z)
}
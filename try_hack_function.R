hack_function(
    orig_fn = glmmTMB::glmmTMB,
    string_target_start =
        c(
            "call <- mf <- mc <- match.call()",
            "mkTMBStruc(formula, ziformula, dispformula,",
            "sparseX = sparseX"
          ),
    string_target_end =
        c(
            "call <- mf <- mc <- match.call()",
            "mkTMBStruc(formula, ziformula, dispformula,",
            "control = control)"
        ),
    string_target_fixed  = FALSE,
    string_target_protect = TRUE,
    action = c("append", "replace","replace"),
    string_replace = list(
        c("phyloZ <- get_phyloZ(phylo, phyloZ, data[[phylonm]])"),
        c("mkTMBStrucphylo(formula, ziformula, dispformula,"),
        c("sparseX = sparseX, phylonm = phylonm, phyloZ = phyloZ)",
          "TMBStruc$control <- lapply(control, eval, envir = TMBStruc)")
    ), ## string_replace
    add_args = list(phylo = NULL, phyloZ = NULL, phylonm = NULL)
)

teststr <- 'if (REML) randomArg <- c(randomArg, "beta")'
grep("if *\\(REML\\) *randomArg *<- *", teststr)

hack_function(
    orig_fn = glmmTMB:::mkTMBStruc,
    string_target_start =
        c(
            "condList  <- getXReTrms",
            'randomArg <- c(randomArg, "beta")'
        ),
    string_target_end =
        c(
            "contrasts = contrasts, sparse = sparseX[[\"cond\"]])",
            'randomArg <- c(randomArg, "beta")'
        ),
    string_target_fixed  = FALSE,
    string_target_protect = TRUE,
    action = c("replace","append"),
    string_replace = list(
        c("condListphylo <- getXReTrmsphylo(formula, mf, fr, contrasts = contrasts, phylonm = phylonm, phyloZ = phyloZ)"),
        c("n.edge <- ncol(phyloZ)",
          'n.site <- length(unique(fr[["site"]]))',
          'relength <- length(condReStruc)',
          'if (relength == 1) {',
          'REname <- unlist(strsplit(names(condReStruc), " "))',
          'rightbar <- REname[length(REname)]',
          'if (rightbar %in% phylonm) {',
          'condReStruc[[1]]$blockReps <- n.edge',
          'data.tmb$terms[[1]]$blockReps <- n.edge',
          '}',
          '}',
          'if (relength > 1) {',
          'for (i in seq_along(condReStruc)) {',
          'REname <- unlist(strsplit(names(condReStruc)[i], " "))',
          'rightbar <- REname[length(REname)]',
          'if (rightbar %in% phylonm) {',
          'condReStruc[[i]]$blockReps <- n.edge',
          'data.tmb$terms[[i]]$blockReps <- n.edge',
          'if (rightbar == "sp:site") {',
          'data.tmb$terms[[i]]$blockReps <- n.edge * n.site ## FIXME: pull out number of site from somewhere',
          '}',
          '}',
          '}',
          '}',
          'condList$Z <- t(condListphylo$reTrms$Zt)',
          'data.tmb$Z <- t(condListphylo$reTrms$Zt)',
          'parameters$b <- rep(0, ncol(data.tmb$Z))'
          )
    ), ## string_replace
    add_args = list(phylo = NULL, phyloZ = NULL, phylonm = NULL)
)

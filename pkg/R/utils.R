##' construct a Z matrix from a phylogeny
##' @name phylo_machinery
##' @param r a \code{phylo} object
##' @param stand standardize edge lengths by determinant of phylogenetic covariance matrix
##' @export
##' @importFrom Matrix t diag
##' @importFrom ape vcv
## FIXME:: rename to phylo_to_Z ?
phylo.to.Z <- function(r, stand = FALSE) {
  ntip <- length(r$tip.label)
  Zid <- Matrix::Matrix(0.0, ncol = length(r$edge.length), nrow = ntip)
  nodes <- (ntip + 1):max(r$edge)
  root <- nodes[!(nodes %in% r$edge[, 2])]
  for (i in 1:ntip) {
    cn <- i ## current node
    while (cn != root) {
      ce <- which(r$edge[, 2] == cn) ## find current edge
      Zid[i, ce] <- 1 ## set Zid to 1
      cn <- r$edge[ce, 1] ## find previous node
    }
  }
  tZid <- t(Zid)
  Z <- t(sqrt(r$edge.length) * tZid)
  if (stand) {
    V <- ape::vcv(r)
    ## V <- V/max(V)
    sig <- exp(as.numeric(determinant(V)["modulus"]) / ntip)
    ## sig <- det(V)^(1/ntip)
    ## all.equal(Z/sqrt(sig),
    ##  t(sqrt(r$edge.length / sig) * tZid), tolerance=2e-16)
    Z <- Z/sqrt(sig)
  }
  rownames(Z) <- r$tip.label
  colnames(Z) <- 1:length(r$edge.length)
  return(Z)
}

check_phylo_names <- function(phyloZ, data_sp) {
  in_phy_nin_data <- setdiff(rownames(phyloZ), unique(data_sp))
  in_data_nin_phy <- setdiff(unique(data_sp), rownames(phyloZ))
  if (length(in_phy_nin_data) > 0 ||
      length(in_data_nin_phy) > 0) {
    msg <- "mismatch between tip names in phylogeny/Z-matrix and names in data:"
    if (length(in_phy_nin_data) > 0) {
      msg <- c(msg,
               "\nin phyloZ but not data: ",
               paste(in_phy_nin_data, collapse = ", "))
    }
    if (length(in_data_nin_phy) > 0) {
      msg <- c(msg,
               "\nin data but not phyloZ: ",
               paste(in_data_nin_phy, collapse = ", "))
    }
    do.call(stop, as.list(msg))
  }
  ## ok
  return(NULL)
}


## check consistency etc. for phyloZ vs species data
get_phyloZ <- function(phylo = NULL, phyloZ = NULL, data_sp) {
  if (is.null(phylo) && is.null(phyloZ)) {
    stop("must provide either phylo or phyloZ")
  }
  if (!is.null(phylo) && !is.null(phyloZ)) {
    warning("both phylo and phyloZ provided: overwriting phyloZ")
  }
  if (!is.factor(data_sp)) {
    stop("species identifier in data must be a factor")
  }
  if (is.null(phylo)) {
    check_phylo_names(phyloZ, data_sp)
    if (!identical(rownames(phyloZ), levels(data_sp))) {
      stop("row names of phyloZ must match values and order of levels of the species identifier")
      ## FIXME: identify mismatches?
    }
    return(phyloZ)
  }
  phyloZ <- phylo.to.Z(phylo)
  check_phylo_names(phyloZ, data_sp)
  phyloZ <- phyloZ[levels(data_sp), ] ## match order
  return(phyloZ)
}

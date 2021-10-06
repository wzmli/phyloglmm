install_pkgs <- function(fn="pkglist") {
  all_pkgs <- scan(fn, what=character(1), quiet=TRUE)
  ## allow/remove comments?
  cran_pkgs <- grep("/", all_pkgs, invert=TRUE, value = TRUE)
  remote_pkgs <- grep("/", all_pkgs, value=TRUE)
  sapply(cran_pkgs, crancache::install_packages)
  sapply(remote_pkgs, remotes::install_github)
}

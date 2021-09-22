tt <- tempdir()
download.file("http://mpcm-evolution.com/OPM/Chapter11_OPM/data.zip",
              dest=file.path(tt, "data.zip"))
old.dir <- setwd(tt)
untar("data.zip")
garamszegi_phy <- ape::read.nexus("phylo.nex")
ff <- list.files(pattern="data_.*\\.txt")
for (f in ff) {
  x <- read.table(f, header=TRUE)
  n <- gsub("\\.txt$","",
            gsub("data_","garamszegi_",f))
  assign(n, x)
}
setwd(old.dir)
save(list = ls(pattern="garamszegi"), file="garamszegi.rda")
unlink(tt)

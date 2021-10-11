targetname <- unlist(strsplit(rtargetname,"[.]"))

platform <- targetname[2]
numsite <- targetname[3]
size <- targetname[4]
tree_seed <- as.numeric(targetname[5])

nspp <- switch(size,
               small = 25,
               med = 50,
               large = 100,
               xlarge = 500,
               stop("unknown size",nspp))

nsite <- switch(numsite,
                ss = 1,
                ms = 20,
                stop("unknown numsite", nsite))

rep = 10

datadir <- "datadir/"

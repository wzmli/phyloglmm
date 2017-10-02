targetname <- unlist(strsplit(rtargetname,"[.]"))

platform <- targetname[2]
seed <- as.numeric(targetname[3])

datadir <- "datadir/"

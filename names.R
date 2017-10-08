targetname <- unlist(strsplit(rtargetname,"[.]"))

platform <- targetname[2]
size <- targetname[3]
seed <- as.numeric(targetname[4])

datadir <- "datadir/"

targetname <- unlist(strsplit(rtargetname,"[.]"))

platform <- targetname[2]
size <- targetname[3]
seed <- as.numeric(targetname[4])

if(size == "small"){nspp <- 25}
if(size == "med"){nspp <- 50}
if(size == "large"){nspp <- 100}
if(size == "k"){nspp <- 1000}
datadir <- "datadir/"

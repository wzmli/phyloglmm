targetname <- unlist(strsplit(rtargetname,"[.]"))

platform <- targetname[2]
numsite <- targetname[3]
size <- targetname[4]
seed <- as.numeric(targetname[5])

if(size == "small"){nspp <- 25}
if(size == "med"){nspp <- 50}
if(size == "large"){nspp <- 100}
if(size == "xlarge"){nspp <- 600}


if(numsite == "ss"){nsite <- 1}
if(numsite == "ms"){nsite <- 10}

datadir <- "datadir/"

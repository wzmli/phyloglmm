targetname <- unlist(strsplit(rtargetname,"[.]"))

platform <- targetname[2]
size <- targetname[3]
seed <- as.numeric(targetname[4])

if(size == "small"){nspp <- 20}
if(size == "med"){nspp <- 100}
if(size == "large"){nspp <- 200}
if(size == "test"){nspp <- 50}
if(size == "xlarge"){nspp <- 500}
if(size == "k"){nspp <- 1000}
datadir <- "datadir/"

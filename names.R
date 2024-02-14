library(shellpipes)

pipeStar()

targetname <- unlist(strsplit(pipeStar(),"[.]"))

print(targetname)

numsite <- targetname[1]
size <- targetname[2]
tree_seed <- as.numeric(targetname[3])

size_vec <- c(small=25, med=50, large = 100,
              mlarge = 250, xlarge = 500)

nspp <- size_vec[[size]]
if (is.na(nspp)) stop("unknown size")

nsite <- switch(numsite
	, ss = 1
	, ms = 20
	, mms = 20
	, stop("unknown numsite", nsite)
)

rep = 10

datadir <- "datadir/"

nrep <- 1

if(exists("numsite")){
   if(numsite == "mms"){
      nrep <- 3
   }
}


saveEnvironment()

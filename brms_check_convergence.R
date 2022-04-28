brms_path
tt <- brmsss_res
for(i in 1:length(tt)){
	aa <- readRDS(paste(brms_path,tt[i],sep=""))
	a <- check_divergences(aa[[1]][["fit"]])
	if(!is.null(a)){
	   print(a)
	   print(tt[i])
	}
}

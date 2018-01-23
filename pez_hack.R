Z.i <- matrix(0,nrow=560,ncol=560)

ss <- function(x,i){
  rr <- (x-1)*28 + i
  cc <- (i-1)*20 + x
  return(c(rr,cc))
}

for(w in 1:20){
  for(r in 1:28){
    Z.i[ss(w,r)[1],ss(w,r)[2]] <- 1
  }
}

Zt <- matrix(0,nrow=560,ncol=560)

St.lengths[i] <- 560

St <- matrix(1,nrow=1,ncol=560)

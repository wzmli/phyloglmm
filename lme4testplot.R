library(ggplot2)

gg <- (ggplot(dd, aes(y=estimate,x=size,group=interaction(size,type),fill=type))
  + geom_violin()
  + facet_wrap(~par,scale="free_y")
  + theme_bw()
)
gg

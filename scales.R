library(tidyverse)
library(cowplot)
set.seed(101)
dd <- data_frame(x=rnorm(1000),
                 expx = exp(x),
                 plogisx = plogis(x))
ddt <- dd %>% gather(type,value)                         
ff <- function(x,trans) {
    dd0 <- ddt %>% filter(type==x)  ## select column , as data frame
    ggplot(dd0,aes(x=1,y=value))+geom_violin()+
        scale_y_continuous(trans=trans)+
        facet_wrap(~type)
}
mm <- map2(c("x","expx","plogisx"),c("identity","log","logit"),ff)

do.call(plot_grid,c(mm,list(nrow=1)))

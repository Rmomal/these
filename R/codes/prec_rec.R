library(tidyverse)
library(gridExtra)
library(grid)
precrec <- readRDS("~/these/pack1/R/Simu/PLN/erdos/d/precrec.rds")

p1<-ggplot(precrec[precrec$var==10,],aes(rec,prec,colour=method))+
  geom_smooth()+
  theme_bw()

p2<-ggplot(precrec[precrec$var==16,],aes(rec,prec,colour=method))+
  geom_smooth()+
  theme_bw()

p3<-ggplot(precrec[precrec$var==20,],aes(rec,prec,colour=method))+
  geom_smooth()+
  theme_bw()

p4<-ggplot(precrec[precrec$var==30,],aes(rec,prec,colour=method))+
  geom_smooth()+
  theme_bw()

grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 2, nrow = 2)

ggplot(precrec[ precrec$var==16 &precrec$param==6,],aes(rec,prec,colour=method))+
  geom_line()+
  theme_bw()

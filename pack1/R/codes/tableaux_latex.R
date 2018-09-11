
library(xtable)
variable<-"n"
path<-paste0("/Users/raphaellemomal/simulations/Simu/PLNcov/")
tab_res<-data.frame(var=parameters[[variable]])
for(type in c("erdos","tree","scale-free","cluster")){
  file<-paste0(path,type,"/",variable,"/auc.rds")
  tab<-data.frame(readRDS(file))
  vec<-tab %>% group_by(var) %>% summarise((med=median(treeggm-spiecResid)))
  tab_res<- cbind(tab_res, 100*vec[,2])
  colnames(tab_res)[ncol(tab_res)]<-type
}

xtable(tab_res)


for(variable in c("r")){
  tab_res<-data.frame(var=parameters[[variable]])
  file<-paste0(path,"cluster/",variable,"/auc.rds")
  tab<-data.frame(readRDS(file))
  vec<-tab %>% group_by(var) %>% summarise((med=median(treeggm-spiecResid)))
  tab_res<- cbind(tab_res, 100*vec[,2])
  colnames(tab_res)[ncol(tab_res)]<-type
  print(tab_res)
}


tab_res<-data.frame(prob=parameters[["prob"]])
file<-paste0(path,"erdos/prob/auc.rds")
tab<-data.frame(readRDS(file))
vec<-tab %>% group_by(var) %>% summarise((med=median(treeggm-spiecResid)))
tab_res<- cbind(tab_res, 100*vec[,2])
colnames(tab_res)[ncol(tab_res)]<-"erdos"
print(tab_res)

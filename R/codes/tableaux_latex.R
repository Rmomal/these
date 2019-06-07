
library(xtable)

parameters<-list(c(seq(10,30,2)),c(seq(10,120,10)), c(seq(0.5,5,0.5)/20),
                 c(seq(1,30,5)),c(seq(0.1,0.4,0.05)))
names(parameters)<-c("d","n","prob","r","dens")

variable<-"n"
path<-paste0("/Users/raphaellemomal/simulations/Simu/PLNcov/")
tab_res<-data.frame(var=parameters[[variable]])
for(type in c("erdos","tree","scale-free","cluster")){
  file<-paste0(path,type,"/",variable,"/auc.rds")
  tab<-data.frame(readRDS(file))
  vec<-tab %>% group_by(var) %>% summarise((med=median(EMmarg-spiecResid)))
  tab_res<- cbind(tab_res, 100*vec[,2])
  colnames(tab_res)[ncol(tab_res)]<-type
}


print(xtable(tab_res), include.rownames=FALSE)


for(variable in c("r")){
  tab_res<-data.frame(var=parameters[[variable]])
  file<-paste0(path,"cluster/",variable,"/auc.rds")
  tab<-data.frame(readRDS(file))
  vec<-tab %>% group_by(var) %>% summarise((med=median(EMCond-spiecResid)))
  tab_res<- cbind(tab_res, 100*vec[,2])
  colnames(tab_res)[ncol(tab_res)]<-type
  print(tab_res, include.rownames=FALSE)
}


tab_res<-data.frame(prob=parameters[["prob"]])
file<-paste0(path,"erdos/prob/auc.rds")
tab<-data.frame(readRDS(file))
vec<-tab %>% group_by(var) %>% summarise((med=median(EMmarg-spiecResid)))
tab_res<- cbind(tab_res, 100*vec[,2])
colnames(tab_res)[ncol(tab_res)]<-"erdos"
print(tab_res, include.rownames=FALSE)



datatable=data.frame(names=baran95$species.names, numbers=1:33)
print(xtable(datatable), include.rownames=FALSE)


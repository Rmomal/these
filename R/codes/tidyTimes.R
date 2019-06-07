library(tidyverse)
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
variable<-c("n","d")
methods<-c("EM_features","OneStep_features")
types<-c("cluster","erdos","scale-free")

#n<50
var=variable[1]

list_dftimes<-lapply(methods,function(met){
  lapply(types, function(typ){
    as_tibble(readRDS(paste0(path,typ,"/",var,"/",met,".rds"))) %>%
      filter(as.numeric(valeur_param)>=50) %>% mutate(time=as.numeric(timeFit),
                                                          method=met, type=typ)
  })

})
# bind types
list_dftimes<-lapply(list_dftimes, function(x){
  do.call(rbind,x)
})
dftimes<-  do.call(rbind,list_dftimes)
dftimes %>% group_by(method) %>% summarise(med=round(median(time),2),sd=paste0("(",round(sd(time),2),")"))


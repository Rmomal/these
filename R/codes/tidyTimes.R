library(tidyverse)
library(xtable)
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
variable<-c("n","d")
methods<-c("SpiecEasi","gCoda","ecoCopula","MRFcov","MInt","EMtree")
types<-c("cluster","erdos","scale-free")
diffs=c("easy","hard")


files=c()
for (method in methods){
 for (type in types){
   for(difficulty in diffs)
    files=c(files,paste0("/Users/raphaellemomal/simulations/Simu/PLN.2.0/TPFN/Results/",
                                 method,"_",type,"_TFPN_",difficulty,".rds"))
 }
}



dataTimes<-do.call(rbind, lapply(files,function(x){
  as_tibble(readRDS(x)) %>% 
    dplyr::select(method,type, times, difficulty)
}))


tab_erdoscluster=dataTimes %>% filter(type!="scale-free") %>% group_by(method,difficulty) %>% 
  summarise(med=round(median(times),2),sd=paste0("(",round(sd(times),2),")")) %>% 
  mutate(res=paste0(med,sd)) %>% select(-med,-sd) %>% 
  spread(method,res)

print(xtable(tab_erdoscluster), include.rownames=FALSE)

tab_scale=dataTimes %>% filter(type=="scale-free") %>% group_by(method,difficulty) %>% 
  summarise(med=round(median(times),2),sd=paste0("(",round(sd(times),2),")")) %>% 
  mutate(res=paste0(med,sd)) %>% select(-med,-sd) %>% 
  spread(method,res)

print(xtable(tab_scale), include.rownames=FALSE)
###############
Bigdata=do.call(rbind, lapply(files,function(x){
  as_tibble(readRDS(x)) 
}))
tab1=Bigdata %>% group_by(method,difficulty,type) %>% summarise(medFDR=round(median(FDR, na.rm=TRUE),2),
                                                                sd=paste0("(",round(sd(FDR,na.rm=TRUE),2),")")) %>% 
  mutate(res=paste0(medFDR,sd)) %>% select(-medFDR,-sd) %>% spread(method,res)

tab2=Bigdata %>% group_by(method,difficulty,type) %>% summarise(empty=sum(which(FP==0 & TP==0))/length(FP)) %>% 
  spread(method,empty)
tab3=Bigdata %>% group_by(method,difficulty,type) %>% summarise(meddens=round(median((TP+FP)/(TP+FN), na.rm=TRUE),2),
                                                                sd=paste0("(",round(sd((TP+FP)/(TP+FN),na.rm=TRUE),2),")")) %>% 
  mutate(res=paste0(meddens,sd)) %>% select(-meddens,-sd) %>%  spread(method,res)
print(xtable(tab1,include.rownames=FALSE))
print(xtable(tab3,include.rownames=FALSE))    

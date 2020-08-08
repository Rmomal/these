library(tidyverse)
library(xtable)
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
variable<-c("n","d")
methods<-c("SpiecEasi","gCoda","ecoCo.AIC","MRFcov1","MInt","EMtree")
types<-c("cluster","erdos","scale-free")
diffs=c("easy","hard")


# files=c()
# for (method in methods){
#   for (type in types){
#     for(difficulty in diffs){
#       resultsname=""
#       if(type=="cluster") resultsname="BadC_"
#       if(type=="erdos" && method%in%c("MRFcov1", "ecoCo.AIC","EMtree")) resultsname="Updated_"
#       if(type=="scale-free") resultsname="Updated_"
#       files=c(files,paste0(path,"TPFN/Updated Results/gather_update/",resultsname,
#                            method,"_",type,"_TFPN_",difficulty,".rds"))
#     }
#   }
# }
types<-c("scale-free")
files=c()
for (method in methods){
  for (type in types){
    for(difficulty in diffs){
      for(resultsname in  c("Updated_","BadSF_")){
        files=c(files,paste0(path,"TPFN/Updated Results/gather_update/",resultsname,
                             method,"_",type,"_TFPN_",difficulty,".rds"))
      }
    }
  }
}



dataTimes<-do.call(rbind, lapply(files,function(x){
  resultsname=strsplit(strsplit(x,"/")[[1]][10],"_")[[1]][1]
  as_tibble(readRDS(x)) %>% 
    dplyr::select(method,type, times, difficulty) %>% mutate(result= resultsname)
})) %>% filter(times<200)


tab_update_bad_SF=dataTimes %>%  group_by(method,result,difficulty) %>% 
  summarise(med=round(median(times),2),sd=paste0("(",round(sd(times),2),")")) %>% 
  mutate(res=paste0(med,sd)) %>% dplyr::select(-med,-sd) %>% 
  spread(method,res)

print(xtable(tab_update_bad_SF), include.rownames=FALSE)



tab_erdoscluster=dataTimes %>%filter(type!="scale-free") %>%  group_by(method,difficulty) %>% 
  summarise(med=round(median(times),2),sd=paste0("(",round(sd(times),2),")")) %>% 
  mutate(res=paste0(med,sd)) %>% dplyr::select(-med,-sd) %>% 
  spread(method,res)

print(xtable(tab_erdoscluster), include.rownames=FALSE)

tab_scale=dataTimes %>% filter(type=="scale-free") %>% group_by(method,difficulty) %>% 
  summarise(med=round(median(times),2),sd=paste0("(",round(sd(times),2),")")) %>% 
  mutate(res=paste0(med,sd)) %>% dplyr::select(-med,-sd) %>% 
  spread(method,res)


tab_erdos=dataTimes %>% filter(type=="erdos") %>% group_by(method,difficulty) %>% 
  summarise(med=round(median(times),2),sd=paste0("(",round(sd(times),2),")")) %>% 
  mutate(res=paste0(med,sd)) %>% dplyr::select(-med,-sd) %>% 
  spread(method,res)


tab_cluster=dataTimes %>% filter(type=="cluster") %>% group_by(method,difficulty) %>% 
  summarise(med=round(median(times),2),sd=paste0("(",round(sd(times),2),")")) %>% 
  mutate(res=paste0(med,sd)) %>% dplyr::select(-med,-sd) %>% 
  spread(method,res)

print(xtable(tab_scale), include.rownames=FALSE)
###############
Bigdata=do.call(rbind, lapply(files,function(x){
  resultsname=strsplit(strsplit(x,"/")[[1]][10],"_")[[1]][1]
  as_tibble(readRDS(x)) %>%   mutate(result= resultsname)
}))
tab1=Bigdata %>%filter(result=="BadSF") %>%  group_by(method,difficulty,type) %>% summarise(medFDR=round(median(FDR, na.rm=TRUE),2),
                                                                sd=paste0("(",round(sd(FDR,na.rm=TRUE),2),")")) %>% 
  mutate(res=paste0(medFDR,sd)) %>% dplyr::select(-medFDR,-sd) %>% spread(method,res)

tab2=Bigdata%>%filter(result=="BadSF") %>%  group_by(method,difficulty,type) %>% summarise(empty=sum(which(FP==0 & TP==0))/length(FP)) %>% 
  spread(method,empty)
tab3=Bigdata %>%filter(result=="BadSF") %>% group_by(method,difficulty,type) %>% summarise(meddens=round(median((TP+FP)/(TP+FN), na.rm=TRUE),2),
                                                                sd=paste0("(",round(sd((TP+FP)/(TP+FN),na.rm=TRUE),2),")")) %>% 
  mutate(res=paste0(meddens,sd)) %>% dplyr::select(-meddens,-sd) %>%  spread(method,res)
print(xtable(tab1,include.rownames=FALSE))
print(xtable(tab2,include.rownames=FALSE))
print(xtable(tab3,include.rownames=FALSE))    

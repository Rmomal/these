
type<-"erdos"
variable<-"d"
EM_features <- readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLN.2.0/",type,
                              "/",variable,"/EM_features.rds"))


EM_features %>% 
  ggplot(aes(valeur_param,time))+
  theme_minimal()+
  geom_boxplot()
EM_features %>% 
  filter(valeur_param%in%c("26","28","30")) %>% 
  ggplot(aes(time,alpha,color=valeur_param))+
  theme_minimal()+
  geom_point()

files<-c()
for(type in c("erdos","cluster","scale-free","tree")){
  cparam<-switch(type,"erdos"=c("n","d","prob"),"tree"=c("n","d"),
                 "cluster"=c("n","d","dens","r"),"scale-free"=c("n","d"))
  for(variable in cparam){
    for(method in c("EM","OneStep","gcodaResid","spiecResid"))
      if(variable=="n" & method=="gcodaResid"){
        print("hourra")
        files<-c(files,paste0(path,type,"/",variable,"/",method,"_features2.rds"),
                 paste0(path,type,"/",variable,"/",method,"_features.rds")) 
      }else{
        print("hourri")
        files<-c(files,paste0(path,type,"/",variable,"/",method,"_features.rds"))
      }
  }
}

Fat<-cbind(readRDS(files[1])[,c("time","valeur_param")],method="EMcv",type="Erdös",variable="n")
sapply(files[-1], function(x){
  elmts<-unlist(strsplit(x,"/"))
  variables<-elmts[length(elmts)-1]
  variables<-switch(variables,"n"="n","d"="p","prob"="edge probability","dens"="density","r"="ratio")
  types<-elmts[length(elmts)-2]
  method<-elmts[length(elmts)]
  #browser()
  method<-switch(method,"EM_features.rds"="EMcv","OneStep_features.rds"="EM1",
                 "gcodaResid_features.rds"="Gcoda","gcodaResid_features2.rds"="Gcoda",
                 "spiecResid_features.rds"="SpiecEasi")
  
  types<-switch(types,"erdos"="Erdös","cluster"="Cluster","scale-free"="Scale-free","tree"="Tree")
  Fat<<-rbind(Fat,cbind(readRDS(x)[,c("time","valeur_param")],method=method,type=types,variable=variables))
})     
saveRDS(Fat,paste0(path,"Fat_times.rds"))      

indices<-which(Fat$variable=="n" & (Fat$valeur_param%in%c("10")))
Fat<-Fat[-indices,]
Fat %>% as_tibble() %>% 
  mutate(valeur_param=as.numeric(valeur_param)) %>% 
  filter(method%in%c("EM1","EMcv","SpiecEasi") ) %>% 
  ggplot(aes(x=as.factor(valeur_param),y=time,color=method))+
  geom_boxplot()+
  theme_minimal()+
  scale_color_manual("",values=c("steelblue4","steelblue2","orange2","indianred2"))+
  facet_wrap(variable~type,scales="free")+
  labs(x="",y="Running times (s)")+
  theme(strip.text  = element_text(size=12))

#### tableau des médianes
library(xtable)


variable<-"n"
tab_res<-Fat %>% as_tibble() %>% 
  filter(variable%in%c("edge probability"))%>%
  mutate(valeur_param=as.numeric(valeur_param)) %>% 
 # mutate(class=cut(valeur_param,breaks=c(18,30))) %>% 
  group_by(method,variable) %>% 
  summarise(med=round(mean(time),2)) %>% 
  #unite(melt, c(type,method)) %>% 
  spread(variable,med)
tab_res<-tab_res[,c(1,3,2)]
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






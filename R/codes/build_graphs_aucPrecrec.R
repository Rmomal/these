############
#   AUC
############
library(ROCR)
#source("fonctions.R") # notemment build_crit

loop_auc<-function(method, path){
  for(type in c("cluster","erdos")){
    cparam<-switch(type,"scale-free"=c("d"),"erdos"=c("n","d","prob"), "cluster"=c("d","n","dens","r"))
    for(variable in cparam){

      length<-length(parameters[[variable]])
      df<-data.frame(valeur=rep(parameters[[variable]],Bgraph),graph=rep(1:Bgraph,each=length),
                     meth=rep(0,Bgraph*length))
      #fill df
      sapply(parameters[[variable]],function(x){
        print(paste0("type: ", type," // var: ", variable," // valeur :",x))
        lapply(1:Bgraph,function(nbgraph){
          path2<-paste0(path,type,"/",variable)
          df$meth[which(df$valeur==x & df$graph==nbgraph)]<<-build_crit(path2,nbgraph,x, variable,
                                                                        method=method,crit="auc")
        })
      })
      colnames(df)[ncol(df)]<-method
      #browser()
      #save df in previous auc file
      if(method=="EMCond"){
        saveRDS(df,paste0(path,type,"/",variable,"/auc.rds"))
      }else{
        auc <- readRDS(paste0(path,type,"/",variable,"/auc.rds"))
        auc<-left_join(auc,df,by=c("valeur","graph"))
        saveRDS(auc,paste0(path,type,"/",variable,"/auc.rds"))
      }
    }
  }
}
#======#
# RUN
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
path_nonfav<-"/Users/raphaellemomal/simulations/Simu/PLN_nonfav/"
method<-c("EMCond", "OneCond", "spiecResid","gcodaResid")
T1<-Sys.time()
sapply(method, function(x) loop_auc(x,path=path_nonfav))
T2<-Sys.time()
difftime(T2, T1)
Fat<-cbind(readRDS("/Users/raphaellemomal/simulations/Simu/PLN_nonfav/cluster/d/test/auc.rds"),
           type="cluster",variable="density")
# Fat AUC
# files<-c("/Users/raphaellemomal/simulations/Simu/PLN.2.0/scale-free/n/auc.rds",
#         "/Users/raphaellemomal/simulations/Simu/PLN.2.0/scale-free/d/auc.rds",
#         "/Users/raphaellemomal/simulations/Simu/PLN.2.0/scale-free/d/auc_n30.rds")

files<-c()
for(type in c("erdos","cluster")){
  cparam<-switch(type,"scale-free"=c("d"),"erdos"=c("n","d","prob"), "cluster"=c("d","n","dens","r"))
  for(variable in cparam){
    files<-c(files,paste0(path_nonfav,type,"/",variable,"/auc.rds"))
  }
}
# dens/ prob :
files<-list("/Users/raphaellemomal/simulations/Simu/PLN.2.0/erdos/prob/auc.rds",
         "/Users/raphaellemomal/simulations/Simu/PLN_nonfav/cluster/dens/auc.rds")

Fat<-do.call(rbind,lapply(files, function(x){
  elmts<-unlist(strsplit(x,"/"))
  variables<-elmts[length(elmts)-1]
  #if(variables=="dens") browser()
  variables<-switch(variables,"n"="n","d"="p","prob"="edge probability","dens"="density","r"="ratio")
  types<-elmts[length(elmts)-2]
  types<-switch(types,"erdos"="Erdös","cluster"="Cluster","scale-free"="Scale-free","tree"="Tree")
  if(elmts[length(elmts)]=="auc_n30.rds") variables="dn30"
  cbind(readRDS(x)[,c("valeur", "graph", "EMCond", "OneCond", "spiecResid", "gcodaResid")],type=types,variable=variables)
}))


saveRDS(Fat,paste0(path,"Fat_auc.rds"))


############
# Graphics
############

Fat_fav<-readRDS(paste0(path,"Fat_auc.rds"))
Fat_nonfav<-readRDS(paste0(path_nonfav ,"Fat_auc.rds"))
Fat_nonfav<-Fat_nonfav %>%
  mutate(type=paste0("Dense ",type))
Fat<-rbind(Fat_fav,Fat_nonfav)
graph_auc<-function(type,variable,lines=FALSE){
  file<-paste0(path,type,"/",variable,"/auc.rds")

  # tab<-data.frame(readRDS(file))
  # lignes<-which(is.na(tab[,4]))
  # if (length(lignes)!=0) tab<-tab[-lignes,]
  # browser()
  tab<-Fat %>%  as_tibble() %>%select(-OneCond) %>% filter(variable%in%c("density","edge probability")) %>%
   # mutate(variable=fct_recode(Fat$variable, "probability"="p","density"="dens")) %>% #filter(variable%in%c("n","p","ratio") & type=="Scale-free") %>%
    #filter(variable%in%c("edge probability","density","ratio")) %>%
   gather( method,value,EMCond,gcodaResid,spiecResid) %>%  group_by(type,variable,valeur,method) %>%
   # gather( method,value,EMCond, OneCond,gcodaResid,spiecResid) %>%  group_by(valeur,method) %>% filter(!is.na(value)) %>%
    summarise(mns=median(value),inf=quantile(value,0.25),sup=quantile(value,0.75)) #%>%
    #group_by(type,variable)
   #%>%
   # mutate(xmin=min(inf))
   # tab<-summarise(group_by(tab,valeur,method),mns=median(value),inf=quantile(value,0.25),sup=quantile(value,0.75))

  xmin=min(tab$inf)
  p<-ggplot(tab, aes(y=mns,x=valeur,shape=method,color=method))+
    geom_errorbar(aes(ymin=inf, ymax=sup), width=0,position=position_dodge(0))+
    geom_point(size=2.5)+
    labs(y="AUC",x="")+
    scale_shape_manual(values=c(16,15,17),breaks=c("EMCond","gcodaResid","spiecResid"),
                       labels=c("EMtree","gCoda","SpiecEasi"))+
    scale_color_manual(values=c("steelblue4","orange2","indianred2") ,breaks=c("EMCond","gcodaResid","spiecResid"),
                       labels=c("EMtree","gCoda","SpiecEasi"))+
    scale_y_continuous(limits = c(xmin,1))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())
   p<-p+geom_line(size=0.2)
  if(xmin<0.5) p<-p+geom_abline(slope=0,intercept=0.5,color="black", linetype="dashed")
  p1<-p+facet_wrap(type~variable, scales = "free_x",nrow=1)+
    theme(strip.text  = element_text(size=12))+ theme(panel.spacing = unit(1, "lines"))
    p1<-p1+geom_vline(xintercept = c(0.1,0.25), linetype="dashed", color="forestgreen", size=0.5)

  ggsave("panel_dens_seuils.png", plot=p1, width=6, height=3.5)
  #p+facet_grid(type~variable,scales="free")
}

grid_arrange_shared_legend(p_fav,p_nonfav,nrow=2,ncol=1,heights=c(5,1))
type<-"cluster"
variable<-"n"
graph_auc(type,variable, lines=TRUE)


makegraph_auc<-function(type, variable){
  pdf(paste0(path,"images/aucs",type,"_",variable,x,".pdf"),
      width=6,
      height=4,onefile=TRUE)
  print(graph_auc(type,variable,x))
  dev.off()
}

for(type in c("erdos","cluster","scale-free","tree")){
  cparam<-switch(type,"erdos"=c("n","prob"),"tree"=c("n"),
                 "cluster"=c("n","dens","r"),"scale-free"=c("n"))
  for(variable in cparam){
    makegraph_auc(type, variable)
  }
}


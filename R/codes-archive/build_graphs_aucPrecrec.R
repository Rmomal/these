############
#   AUC
############
library(ROCR)
#source("fonctions.R") # notemment build_crit

loop_auc<-function(method, path){
  for(type in c("erdos","cluster")){
    cparam<-switch(type,"cluster"=c("n","d","r"),"erdos"=c("n","d"))
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
      #save df in previous auc file
      if(method=="EMmarg"){
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
method<-c("EMmarg","EMCond", "OneCond","OneMarg","gcodaResid","spiecResid")
T1<-Sys.time()
sapply(method, function(x) loop_auc(x,path=path))
T2<-Sys.time()
difftime(T2, T1) 

# Fat AUC
files<-c()
for(type in c("erdos","cluster")){
  cparam<-switch(type,"cluster"=c("n","d","r"),"erdos"=c("n","d"))
  for(variable in cparam){
    files<-c(files,paste0(path,type,"/",variable,"/auc.rds"))
  }
}

Fat<-cbind(readRDS(files[1]),type="Erdös",variable="n")
sapply(files[-1], function(x){
  elmts<-unlist(strsplit(x,"/"))
  variables<-elmts[length(elmts)-1]
  variables<-switch(variables,"n"="n","d"="p","prob"="edge probability","dens"="density","r"="ratio")
  types<-elmts[length(elmts)-2]
  types<-switch(types,"erdos"="Erdös","cluster"="Cluster","scale-free"="Scale-free","tree"="Tree")
  Fat<<-rbind(Fat,cbind(readRDS(x),type=types,variable=variables))
})     


saveRDS(Fat,paste0(path,"Fat_auc.rds"))      


############
#  PrecRec
############

#####
#precrec simple
#####
loop_precrec<-function(path, methods){
  for(type in c("erdos","tree","cluster","scale-free")){
    cparam<-switch(type,"erdos"=c("d","n","prob"),"tree"=c("d","n"),
                   "cluster"=c("d","n","dens","r"),"scale-free"=c("d","n"))
    for(variable in cparam){
      print(paste0(type," // ", variable," // ",parameters[[variable]]))
      sapply(parameters[[variable]],function(x){
        precrec<-data.frame(prec=double(),rec=double(),method=character(),valeur=double(),graph=integer(),
                            stringsAsFactors=FALSE)
        path2<-paste0(path,type,"/",variable)
        sapply(1:Bgraph, function(nbgraph){
          sapply(methods, function(method){
            
            eval_crit<-build_crit(path2, nbgraph,x,variable,method, crit="precrec")
            
            tmp<-data.frame(eval_crit[[1]],eval_crit[[2]],method,x,nbgraph)
            colnames(tmp)<-c("prec","rec","method","valeur","graph")
            precrec<<- rbind(precrec,tmp)
          })
        }) 
     #   browser()
        saveRDS(precrec,paste0(path2,"/precrec/precrecSimple","_",x,".rds"))
      })
    }
  }
}
#======#
# RUN
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
method<-c("EMmarg", "EMCond", "OneMarg", "OneCond","gcodaResid","spiecResid")
T1<-Sys.time()
loop_precrec(path, method)
T2<-Sys.time()
difftime(T2, T1) 


#####
# precrec pool
#####
loop_precrecPool<-function(path,methods){
  for(type in c("erdos","tree","cluster","scale-free")){
    cparam<-switch(type,"erdos"=c("d","n","prob"),"tree"=c("d","n"),
                   "cluster"=c("d","n","dens","r"),"scale-free"=c("d","n"))
    for(variable in cparam){
      #fill df
      sapply(parameters[[variable]],function(x){
        print(paste0("type: ", type," // var: ", variable," // valeur :",x))
        path2<-paste0(path,type,"/",variable)
        fatListe_methodes<-lapply(methods,function(method){
          #c("treeggm","one_step","glasso","oracle","spiecResid")
          grosseListe<-lapply(1:Bgraph,function(nbgraph){
            build_crit(path2,nbgraph,x,variable,method=method,crit="precrecPool")
          })
          #les scores pour chaque croisement type*variable*valeur*method sont transformés et 
          # concatenes sur tous les graphes generes. On obtient deux grands vecteurs de prédictions et d'observations
          vec_pred<-do.call(c,lapply(grosseListe,function(x) x[[1]]))
          vec_obs<-do.call(c,lapply(grosseListe,function(x) x[[2]]))
          prediction<-prediction(vec_pred,vec_obs)
          # sur lesquels on calcule enfin les stat prec et rec
          ROC_precision<-performance(prediction,"prec")
          ROC_recal<-performance(prediction,"rec")
          precrec<-data.frame(ROC_precision@y.values,ROC_recal@y.values,prediction@cutoffs ,method)
          colnames(precrec)<-c("prec","rec","cut","method")
          return(precrec)
        })
        #les frames pour les diffferentes methodes sont concatenes en vue des plots
        fatfatListe<-as_tibble(do.call(rbind,fatListe_methodes)) #42e3 lignes
        saveRDS(fatfatListe,paste0(path2,"/precrec_pool/precrecPool_",x,".rds"))
      })
    }
  }
}
#======#
# RUN
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
method<-c("EMmarg", "EMCond", "OneMarg", "OneCond","gcodaResid","spiecResid")
T1<-Sys.time()
loop_precrecPool(path, method)
T2<-Sys.time()
difftime(T2, T1) 


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
  tab<-Fat_fav
  # tab<-data.frame(readRDS(file))
  # lignes<-which(is.na(tab[,4]))
  # if (length(lignes)!=0) tab<-tab[-lignes,]
  # browser()
  tab<-tab %>%  as_tibble() %>% filter(variable%in%c("edge probability","density","ratio")) %>% 
    gather( method,value,EMCond,OneCond,gcodaResid,spiecResid) %>%  group_by(type,variable,valeur,method) %>% 
    summarise(mns=median(value),inf=quantile(value,0.25),sup=quantile(value,0.75)) %>% 
    group_by(type,variable) 
  #%>% 
   # mutate(xmin=min(inf))
  #  tab<-summarise(group_by(tab,valeur,method),mns=median(value),inf=quantile(value,0.25),sup=quantile(value,0.75))
  
  xmin=min(tab$inf)
  p<-ggplot(tab, aes(y=mns,x=valeur,shape=method,color=method))+
    geom_errorbar(aes(ymin=inf, ymax=sup), width=0,position=position_dodge(0))+
    geom_point(size=2.5)+
    labs(y="AUC",x="")+
    scale_shape_manual(values=c(16,15,4,17), breaks=c("EMCond","OneCond","gcodaResid","spiecResid"),
                       labels=c("EMcv","EM1","gCoda","SpiecEasi"))+
    scale_color_manual(values=c("steelblue4","orange2","steelblue4","indianred2") ,breaks=c("EMCond","OneCond","gcodaResid","spiecResid"),
                       labels=c("EMcv","EM1","gCoda","SpiecEasi"))+
    scale_y_continuous(limits = c(xmin,1))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())
  lines=TRUE
  if(lines) p<-p+geom_line(size=0.2)
  if(xmin<0.5) p<-p+geom_abline(slope=0,intercept=0.5,color="black", linetype="dashed")
  p+facet_wrap(variable~type, scales = "free",ncol=4)+
    theme(strip.text  = element_text(size=12))
 
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

#####################

graph_precrec<-function(type,variable,x,path,col1,col2,col3){
  type2<-switch(type,"erdos"="Erdös","cluster"="Cluster","scale-free"="Scale-free","tree"="Tree")
  var2<-switch(variable,"d"="p","n"="n","prob"="Edge probability","dens"="density","r"="ratio")
  precrec <- as_tibble(readRDS(paste0(path,type, "/",variable,"/precrec/precrecSimple","_",x,".rds"))) %>% 
    filter(method%in%c("EMCond","gcodaResid"))
  precrecpool<- as_tibble(readRDS(paste0(path,type,"/",variable,"/precrec_pool/precrecPool_",x,".rds"))) %>% 
    filter(method%in%c("EMCond","gcodaResid")) %>%  mutate(method=paste0(method,"pool"))
 superimpose<-rbind(precrec[,c("prec","rec","method")], precrecpool[,c("prec","rec","method")])
  # triche : créer des filtres pour gérer les nuages de points séparément des courbes moyennes
  d_filtered <- superimpose %>%
    group_by(method) %>% 
    filter(method=="EMCondpool" || method=="gcodaResidpool" ) %>%
    ungroup()
  
  d_treeggm <- superimpose %>%
    group_by(method) %>% 
    filter(method=="EMCond") %>%
    ungroup()
  dgcodaResid<- superimpose %>%
    group_by(method) %>% 
    filter(method=="gcodaResid") %>%
    ungroup()
#  browser()
 p<- ggplot(superimpose) +
    theme_minimal()+
    geom_point(aes(rec,prec, group = method),data=dgcodaResid,colour = alpha(col1, 0.5),size=0.5) +
    geom_point(aes(rec,prec, group = method),data=d_treeggm,colour = alpha(col2, 0.2),size=0.5) +
    # colourise only the filtered data
    geom_line(aes(rec,prec, colour = method), data = d_filtered, size=1)+
    scale_color_manual(values=c(col3,col1),   breaks=c("EMCondpool","gcodaResidpool"),
                       labels=c("EMcv ","Gcoda"))+
    theme(legend.position="bottom", legend.title = element_blank(),
          legend.text = element_text(size=12))+
    labs(x="Recall",y="Precision", title=paste0(type2,": ",var2," = ",x))+
    theme(plot.title = element_text(hjust = 0.5))
 p
}
type<-c("erdos","cluster","scale-free","tree")
variable<-c("prob","r","d","n")
x<-c(0.25,21,10,30)
p1<-graph_precrec(type[1],variable[1],x[1],path,"steelblue4","goldenrod1","orange2")
p2<-graph_precrec(type[2],variable[2],x[2],path,"steelblue4","goldenrod1","orange2")
p3<-graph_precrec(type[3],variable[3],x[3],path,"steelblue4","goldenrod1","orange2")
p4<-graph_precrec(type[4],variable[4],x[4],path,"steelblue4","goldenrod1","orange2")

grid_arrange_shared_legend(p1,p2,p3,p4,nrow=2,ncol=2)

makegraph<-function(type, variable, x){
  pdf(paste0(path,"images/precrec/",type,"_",variable,x,".pdf"),
      width=6,
      height=4,onefile=TRUE)
  print(visu_precrec(type,variable,x))
  dev.off()
}

for(type in c("erdos","cluster","scale-free","tree")){
  cparam<-switch(type,"erdos"=c("n","prob"),"tree"=c("n"),
                 "cluster"=c("n","dens","r"),"scale-free"=c("n"))
  for(variable in cparam){
    min_max<-c(min(parameters[[variable]]),max(parameters[[variable]]))
    for(x in min_max){
      makegraph(type, variable, x)
    }
  }
}

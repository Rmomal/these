
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

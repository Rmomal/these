
#######
# amend precrec for each graph with spiecResid
#######
library(ROCR)
parameters<-list(c(seq(10,30,2)),c(seq(10,120,10)), c(seq(0.5,5,0.5)/20),
                 c(seq(1,30,5)),c(seq(0.1,0.4,0.05)))
names(parameters)<-c("d","n","prob","r","dens")
Bgraph<-100
#####
#precrec simple
#####
build_precrec<-function(){
  for(type in c("erdos","tree","scale-free","cluster")){
    cparam<-switch(type,"erdos"=c("n","prob"),"tree"=c("n"),
                   "cluster"=c("n","dens","r"),"scale-free"=c("n"))
    path<-paste0("/Users/raphaellemomal/simulations/Simu/PLNcov/")
    ### variable
    for(variable in cparam){
      ### valeur de la variable
      sapply(parameters[[variable]],function(x){
        print(paste0(type," // ", variable," // ",parameters[[variable]]))
        precrec<-data.frame(prec=double(),rec=double(),method=character(),
                            var=double(),B=integer(),param=integer(),
                            stringsAsFactors=FALSE)
        colnames(precrec)<-c("prec","rec","method","var","B","param")
        path2<-paste0(path,type,"/",variable)
        ### le numero du graph
        sapply(1:Bgraph, function(nbgraph){
          ### la methode utilisee
          sapply(c("_treeggm_","_spiecResid","_oracle"), function(method){
            vec<-build_vec(path2, nbgraph,x,variable,method)
            #####
            
            vec_pred<-vec[[1]]
            vec_obs<-vec[[2]]
            prediction<-prediction(vec_pred,vec_obs)
            ROC_precision<-performance(prediction,"prec")
            ROC_recal<-performance(prediction,"rec")
            #####
            tmp<-data.frame(ROC_precision@y.values,ROC_recal@y.values,method,x,B=1,nbgraph)
            colnames(tmp)<-c("prec","rec","method","var","B","param")
            precrec<<- rbind(precrec,tmp)
          })
        }) 
        
        saveRDS(precrec,paste0(path,type,"/",variable,"/precrec/precrec","_",x,".rds"))
      })
      
    }
  }
}
build_precrec()
#####
# courbes moyennes superosée au nuage de points
#####
path<-"/Users/raphaellemomal/simulations/Simu/PLNcov/"
visu_precrec<-function(type,variable,x,path="/Users/raphaellemomal/simulations/Simu/PLNcov/"){
  
  precrec <- readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLNcov/",type,
                            "/",variable,"/precrec/precrec","_",x,".rds"))
  indices<-which(precrec$method=="_treeggm_" | precrec$method=="_spiecResid")
  precrecpool<- readRDS(paste0(path,type,"/",variable,"/precrec_pool/precrec_",x,".rds"))
  
  indices2<-which(precrecpool$method=="_treeggm_" | precrecpool$method=="_spiecResid")
  precrecpool$method<-paste0(precrecpool$method,"pool")
  
  superimpose<-rbind(precrec[indices,c("prec","rec","method")],precrecpool[indices2,c("prec","rec","method")])
  
  # triche : créer des filtres pour gérer les nuages de points séparément des courbes moyennes
  d_filtered <- superimpose %>%
    group_by(method) %>% 
    filter(method=="_treeggm_pool" || method=="_spiecResidpool" || method=="_oraclepool") %>%
    ungroup()
  d_treeggm <- superimpose %>%
    group_by(method) %>% 
    filter(method=="_treeggm_") %>%
    ungroup()
  d_spiecResid<- superimpose %>%
    group_by(method) %>% 
    filter(method=="_spiecResid") %>%
    ungroup()
  
  ggplot(superimpose) +
    theme_bw()+
    geom_point(aes(rec,prec, group = method),data=d_spiecResid,colour = alpha("#ff7f0e", 0.3),size=0.5) +
    geom_point(aes(rec,prec, group = method),data=d_treeggm,colour = alpha("#d62728", 0.3),size=0.5) +
    # colourise only the filtered data
    geom_line(aes(rec,prec, colour = method), data = d_filtered, size=1)+
    scale_color_manual(values=c("#d62728","#ff7f0e","#fc5f94"), #jaune,bleu, rose
                       breaks=c("_treeggm_pool","_spiecResidpool","_oraclepool"),
                       labels=c("EM ","SpiecResid","Oracle" ))+
    theme(legend.position="bottom", legend.title = element_blank(),
          legend.text = element_text(size=12))+
    labs(x="Recall",y="Precision")
  
  
}


makegraph<-function(type, variable, x){
  pdf(paste0(path,"images/SansOracle",type,"_",variable,x,".pdf"),
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


#################################################################################
# POOL
##############
vec_obs_pred<-function(obs, pred){
  # browser()
  nvar<-ncol(obs)
  obs[which(abs(obs)<1e-16)]<-0
  indices_nuls<-which(obs==0)
  label<-matrix(1,nrow=nrow(obs),ncol=ncol(obs))
  label[indices_nuls]<-0
  
  vec_pred<-as.vector(pred[upper.tri(pred)])
  vec_obs<-as.vector(label[upper.tri(label)])
  
  return(list(vec_pred,vec_obs))
}
build_vec<-function(path2, nbgraph,x,variable,method){
  if(variable=="n"){
    obs<-readRDS(paste0(path2,"/Sets_param/Graph",nbgraph,".rds"))$omega
  }else{
    obs<-readRDS(paste0(path2,"/Sets_param/Graph",nbgraph,"_",x,".rds"))$omega  
  } 
  
  if(method=="_treeggm_"|| method=="_oracle"){
    pred<- readRDS(paste0(path2,"/Scores/Graph",nbgraph,method,x,".rds"))[["probaCond"]] 
  }else{
    pred<- readRDS(paste0(path2,"/Scores/Graph",nbgraph,method,x,".rds"))
  }
  return(vec_obs_pred(obs,pred))
}

build_precrecPool<-function(){
  for(type in c("erdos","tree","scale-free","cluster")){
    # path<-paste0(getwd(),"/Simu/PLNcov/")
    cparam<-switch(type,"erdos"=c("n","prob"),"tree"=c("n"),
                   "cluster"=c("n","dens","r"),"scale-free"=c("n"))
    path<-paste0("/Users/raphaellemomal/simulations/Simu/PLNcov/")
    
    
    for(variable in cparam){
      #fill df
      sapply(parameters[[variable]],function(x){
        print(paste0("type: ", type," // var: ", variable," // valeur :",x))
        fatListe_methodes<-lapply(c("_treeggm_","_spiecResid","_oracle"),function(method){#c("treeggm","one_step","glasso","oracle","spiecResid")
          grosseListe<-lapply(1:Bgraph,function(nbgraph){
            # browser()
            path2<-paste0(path,type,"/",variable)
            build_vec(path2,nbgraph,x,variable,method=method)
          })
          # browser()
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
        saveRDS(fatfatListe,paste0(path,type,"/",variable,"/precrec_pool/precrec_",x,".rds"))
      })
    }
  }
}
build_precrecPool()
##### LOOK PERFORMANCES
# type<-"cluster"
# variable<-"n"
# x<-100
# path<-"/Users/raphaellemomal/simulations/Simu/PLNcov/"
# precrecpool<- readRDS(paste0(path,type,"/",variable,"/precrec_",x,".rds"))
# indices2<-which(precrecpool$method=="treeggm" | precrecpool$method=="glasso")
# ggplot(precrecpool[indices2,],aes(rec,prec,colour=method,shape=method))+
#   geom_point()+
#   scale_shape_manual(values=c(16,15,9,8,17),
#                      breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
#                      labels=c("EM ","1 step","SpiecEasi","spiecResid", "oracle" ) )+
#   scale_color_manual(values=c("#E69F00","#076443", "#8037c9","#fc5f94" ,"#56B4E9"),
#                      breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
#                      labels=c("EM ","1 step","SpiecEasi","spiecResid", "oracle" ) 
#   )+
#   guides(shape = guide_legend(override.aes = list(size = 3)))+
#   labs(title="")+
#   scale_y_continuous(limits = c(0,1))+
#   theme_bw()
# 
# ggplot(precrec,aes(rec,prec,colour=method,linetype=method))+
#   geom_line(size=1)+
#   # scale_type_manual(values=c(16,15,9,8,17),
#   #                    breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
#   #                    labels=c("EM ","1 step","SpiecEasi","spiecResid", "oracle" ) )+
#   scale_linetype_manual(values=c("twodash", "solid", "dashed", "dotted", "dotdash" ),breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
#                         labels=c("EM ","1 step","SpiecEasi","SpiecResid", "Oracle" ) )+
#   scale_color_manual(values=c("#E69F00","#076443", "#8037c9","#fc5f94" ,"#56B4E9"),
#                      breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
#                      labels=c("EM ","1 step","SpiecEasi","SpiecResid", "Oracle" ) 
#   )+
#   guides(shape = guide_legend(override.aes = list(size = 2)))+
#   labs(title="")+
#   scale_y_continuous(limits = c(0,1))+
#   theme_bw()
# 
# 

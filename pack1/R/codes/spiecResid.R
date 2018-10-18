######
# Parameters
######
set.seed(3)
type<-"erdos"
variable<-"d"
nbgraph<-"10"
x<-"16"
path<-paste0("/Users/raphaellemomal/simulations/Simu/PLNcov/")
n<-100
covariables<-cbind(rep(c(0,1),each=n/2),rnorm(n,8,0.5),        round(runif(n)*10))

######
# spieceasi on residuals
######
inf_spiecresid<-function(path, type, variable, nbgraph, x,covariables){
  if(variable=="n"){
    param<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,".rds"))
  }else{
    param<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",x,".rds"))
  } 
  
  n<-ifelse(variable=="n",x,100)
  covariables<-cbind(rep(c(0,1),each=n/2),rnorm(n,8,0.5),        round(runif(n)*10))
  
  Y<-generator_PLN(param$sigma,covariables)[[1]]
  U<-clr(Y)
  model<-lm(U~covariables)
  summary(model)
  resid<-model$residuals
  inf<-inf_spieceasi(resid)
  return(inf)
}

######
# créer et enregistrer les scores
######

parameters<-list(c(seq(10,30,2)),c(seq(10,120,10)), c(seq(0.5,5,0.5)/20),
                 c(seq(1,30,5)),c(seq(0.1,0.4,0.05)))
names(parameters)<-c("d","n","prob","r","dens")
Bgraph<-100

T1<-Sys.time()
for(type in c("erdos","tree","scale-free","cluster")){
  # cparam<-switch(type,"erdos"=c("d","n","prob"),"tree"=c("d","n"),
  #                "cluster"=c("d","n","dens","r"),"scale-free"=c("d","n"))
  cparam<-"n"
  for(variable in cparam){
    sapply(parameters[[variable]],function(x){
      print(paste0("type: ", type," // var: ", variable," // valeur :",x))
      mclapply(1:Bgraph,function(nbgraph){
        inf<-inf_spiecresid(path, type, variable, nbgraph, x, covariables)
        save_file<-paste0(path,type,"/",variable,"/Scores/")
        saveRDS(inf,paste0(save_file,"Graph",nbgraph,"_spiecResid",x,".rds"))
      }, mc.cores=3)
    })
  }
}
T2<-Sys.time()
difftime(T2, T1)


######
# amend auc
######
T1<-Sys.time()
for(type in c("erdos","tree","scale-free","cluster")){
  #path<-paste0(getwd(),"/Simu/PLNcov/")
  # cparam<-switch(type,"erdos"=c("d","n","prob"),"tree"=c("d","n"),
  #                "cluster"=c("d","n","dens","r"),"scale-free"=c("d","n"))
  cparam<-"n"
  for(variable in cparam){
    df<-data.frame(var=rep(parameters[[variable]],Bgraph),param=rep(1:Bgraph,each=length(parameters[[variable]])),
                   spiecResid=rep(0,Bgraph*length(parameters[[variable]])))
    #fill df
    sapply(parameters[[variable]],function(x){
      print(paste0("type: ", type," // var: ", variable," // valeur :",x))
      mclapply(1:Bgraph,function(nbgraph){
        print(nbgraph)
        path2<-paste0(path,type,"/",variable)
        df$spiecResid[which(df$var==x & df$param==nbgraph)]<<-build_auc(path2,nbgraph,x,variable)
      }, mc.cores=1)

    })
    #save df in previous auc file
   
    auc <- readRDS(paste0(path,type,"/",variable,"/auc.rds"))
    auc<-auc[,-ncol(auc)]
    auc<-left_join(auc,df,by=c("var","param"))
    # colnames(auc)[4]<-"spiecEasi"
    saveRDS(auc,paste0(path,type,"/",variable,"/auc.rds"))
  }
}
T2<-Sys.time()
difftime(T2, T1)


build_auc<-function(path, nbgraph,x,variable){
  if(variable=="n"){
    obs<-readRDS(paste0(path,"/Sets_param/Graph",nbgraph,".rds"))$omega
  }else{
    obs<-readRDS(paste0(path,"/Sets_param/Graph",nbgraph,"_",x,".rds"))$omega  
  } 
  pred<- readRDS(paste0(path,"/Scores/Graph",nbgraph,"_spiecResid",x,".rds"))
  return(diagnost_auc(obs,pred))
}
diagnost_auc<-function(obs, pred){
  nvar<-ncol(obs)
  obs[which(abs(obs)<1e-16)]<-0
  indices_nuls<-which(obs==0)
  label<-matrix(1,nrow=nrow(obs),ncol=ncol(obs))
  label[indices_nuls]<-0
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),
                         as.vector(label[upper.tri(label)]))
  obs<-as.vector(label[upper.tri(label)])
  # Run the AUC calculations
  ROC_auc <- performance(prediction,"auc")
  res<-round(ROC_auc@y.values[[1]],digits=3)
  return(res)
}
### see performances
diagnostics<-function(file){
  file<-paste0(path,type,"/", variable,"/auc.rds" )
  tab<-data.frame(readRDS(file))
  lignes<-which(is.na(tab[,"EMCond"]))
  if (length(lignes)!=0) tab<-tab[-lignes,]
  tab<- gather(tab,key=method,value=value,EMCond,EM1Cond,spiecEasi, spiecResid)
  tab<-summarise(group_by(tab,var,method),mns=median(value),
                 inf=quantile(value,0.25),sup=quantile(value,0.75))
  elmts<-unlist(strsplit(file,"/"))
  variable<-elmts[length(elmts)]
  param<-elmts[length(elmts)-1]
  type<-elmts[length(elmts)-2]
  variable<-substr(variable,1,nchar(variable)-4)
  variable<-switch(variable,"auc"="AUC","sens"="Sensitivity","spec"="Specificity")

  tab$var<-as.numeric(as.character(tab$var))
  p<-ggplot(tab, aes(y=mns,x=as.numeric(as.character(var)),shape=method,color=method))+
    geom_errorbar(aes(ymin=inf, ymax=sup), width=0,position=position_dodge((max(tab$var)-min(tab$var))/100))+
    geom_point(size=2)+
    geom_line(size=0.2)+
   labs(y=variable,x=param)+
    scale_shape_manual(values=c(15,16,17,9),
                       breaks=c("EMCond","EM1Cond","spiecResid", "spiecEasi" ),
                       labels=c("EM ","1 step","SpiecResid" ,"SpiecEasi") )+
    scale_color_nejm(breaks=c("EMCond","EM1Cond","spiecResid", "spiecEasi" ),
                      labels=c("EM ","1 step","SpiecResid","SpiecEasi" ))+
    scale_y_continuous(limits = c(0.5,1))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())

  print(p)
}
type<-"tree"
variable<-"n"
propNA<-function(type, variable){
  file<-paste0(path,type,"/", variable,"/auc.rds" )
tab<-readRDS(file)
print(length(which(is.na(tab$EM1Cond)))/nrow(tab))
}
for(type in c("erdos","tree","cluster","scale-free")){
  cparam<-switch(type,"erdos"=c("prob","d", "n"),"tree"=c("n"),
                 "cluster"=c("n","dens","r"),"scale-free"=c("n"))
  
  for(variable in cparam){
    propNA(type,variable)
  }
}
diagnostics(file)
for(type in c("erdos","tree","cluster","scale-free")){
  cparam<-switch(type,"erdos"=c("prob", "n"),"tree"=c("n"),
                 "cluster"=c("n","dens","r"),"scale-free"=c("n"))
  for(variable in cparam){
    pdf(paste0(path,"images/SansOracle/",type,"_",variable,".pdf"),
        width=6,
        height=4,onefile=TRUE)
    file<-paste0(path,type,"/", variable,"/auc.rds" )
    diagnostics(file)
    dev.off()
  }
}
######
# precrec pool
######
#  pooler sur 100 graphes

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
  if(method=="oracle" || method=="spiecResid"){
    pred<- readRDS(paste0(path2,"/Scores/Graph",nbgraph,"_",method,x,".rds"))
  }else{
    pred<- readRDS(paste0(path2,"/Scores/Graph",nbgraph,"_",method,"_",x,".rds")) 
  }
  return(vec_obs_pred(obs,pred))
}

for(type in c("erdos","tree","scale-free","cluster")){
  # path<-paste0(getwd(),"/Simu/PLNcov/")
  # cparam<-switch(type,"erdos"=c("d","n","prob"),"tree"=c("d","n"),
  #                "cluster"=c("d","n","dens","r"),"scale-free"=c("d","n"))
  path<-paste0("/Users/raphaellemomal/simulations/Simu/PLNcov/")
  
  cparam<-"n"
  for(variable in cparam){
    #fill df
    sapply(parameters[[variable]],function(x){
      print(paste0("type: ", type," // var: ", variable," // valeur :",x))
      fatListe_methodes<-lapply(c("treeggm","one_step","glasso","oracle","spiecResid"),function(method){
        grosseListe<-lapply(1:Bgraph,function(nbgraph){
          # browser()
          path2<-paste0(path,type,"/",variable)
          build_vec(path2,nbgraph,x,variable,method=method)
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
      saveRDS(fatfatListe,paste0(path,type,"/",variable,"/precrec_pool/precrec_",x,".rds"))
    })
  }
}

##### LOOK PERFORMANCES
type<-"cluster"
variable<-"n"
x<-100
precrec<- readRDS(paste0(path,type,"/",variable,"/precrec/precrec_",x,".rds"))

ggplot(precrec,aes(rec,prec,colour=method,shape=method))+
  geom_point()+
  scale_shape_manual(values=c(16,15,9,8,17),
                     breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
                     labels=c("EM ","1 step","SpiecEasi","spiecResid", "oracle" ) )+
  scale_color_manual(values=c("#E69F00","#076443", "#8037c9","#fc5f94" ,"#56B4E9"),
                      breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
                      labels=c("EM ","1 step","SpiecEasi","spiecResid", "oracle" ) 
                     )+
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  labs(title="")+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()

ggplot(precrec,aes(rec,prec,colour=method,linetype=method))+
  geom_line(size=1)+
  # scale_type_manual(values=c(16,15,9,8,17),
  #                    breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
  #                    labels=c("EM ","1 step","SpiecEasi","spiecResid", "oracle" ) )+
  scale_linetype_manual(values=c("twodash", "solid", "dashed", "dotted", "dotdash" ),breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
                   labels=c("EM ","1 step","SpiecEasi","SpiecResid", "Oracle" ) )+
  scale_color_manual(values=c("#E69F00","#076443", "#8037c9","#fc5f94" ,"#56B4E9"),
                     breaks=c("treeggm","one_step", "glasso","spiecResid" ,"oracle"),
                     labels=c("EM ","1 step","SpiecEasi","SpiecResid", "Oracle" ) 
  )+
  guides(shape = guide_legend(override.aes = list(size = 2)))+
  labs(title="")+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()





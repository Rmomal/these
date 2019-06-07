
parameters<-list(c(seq(10,30,2)),c(seq(10,120,10)), c(seq(1,5,0.5)/20),
                 c(seq(1,30,5)),c(seq(0.1,0.4,0.05)))
names(parameters)<-c("d","n","prob","r","dens")
Bgraph<-100
path<-"/Users/raphaellemomal/simulations/Simu/PLNcov/"
library(parallel)
library(mvtnorm)
library(PLNmodels)
T1<-Sys.time()
for(type in c("erdos","tree","scale-free","cluster")){
  cparam<-switch(type,"erdos"=c("d","n","prob"),"tree"=c("d","n"),
                 "cluster"=c("d","n","dens","r"),"scale-free"=c("d","n"))
#   cparam<-"r"
  for(variable in cparam){
  
    sapply(parameters[[variable]],function(x){
      print(paste0("type: ", type," // var: ", variable," // valeur :",x))
      mclapply(1:Bgraph,function(nbgraph){
        print(nbgraph)
        save_file<-paste0(path,type,"/",variable,"/Scores/")
        inf<-from_stored_graphs(type,variable,nbgraph,x) 
        saveRDS(inf,paste0(save_file,"Graph",nbgraph,"_one_step_",x,".rds"))
      }, mc.cores=3)
    })
  }
}
T2<-Sys.time()
difftime(T2, T1) # 5h50 vs 1h20

type<-"cluster"
variable<-"r"
nbgraph<-"15"
x<-16
readRDS(paste0(paste0(path,type,"/",variable,"/Scores/"),"Graph",nbgraph,"_treeggm_",x,".rds"))

########
# oracle version
########
library(parallel)
T1<-Sys.time() #all but tree n
for(type in c("cluster")){
  # cparam<-switch(type,"erdos"=c("prob","d","n"),"tree"=c("n","d"),
  #                "cluster"=c("dens","d","n","r"),"scale-free"=c("d","n"))
 cparam<-"r"
  for(variable in cparam){
    sapply(parameters[[variable]],function(x){
      print(paste0("type: ", type," // var: ", variable," // valeur :",x))
      mclapply(1:Bgraph,function(nbgraph){
        inf<-from_stored_graphs(type,variable,nbgraph,x,oracle=TRUE) 
        save_file<-paste0(path,type,"/",variable,"/Scores/")
        saveRDS(inf,paste0(save_file,"Graph",nbgraph,"_oracle",x,".rds"))
      }, mc.cores=3)
      
    })
  }
}
T2<-Sys.time()
difftime(T2, T1) # 2h sans erdos


########
# refaire les  AUC et precrec
#######
T1<-Sys.time()
for(type in c("scale-free")){
  cparam<-switch(type,"erdos"=c("n","prob"),"tree"=c("n"),
                 "cluster"=c("n","dens","r"),"scale-free"=c("n"))
     cparam<-"d"
  for(variable in cparam){
    
    length<-length(parameters[[variable]])

    df<-data.frame(var=rep(parameters[[variable]],Bgraph),param=rep(1:Bgraph,each=length),
                   MInt=rep(0,Bgraph*length))
    #fill df
    sapply(parameters[[variable]],function(x){
      print(paste0("type: ", type," // var: ", variable," // valeur :",x))
      lapply(1:Bgraph,function(nbgraph){
        print(nbgraph)
       path2<-paste0(path,type,"/",variable)
        df$MInt[which(df$var==x & df$param==nbgraph)]<<-build_auc(path2,nbgraph,x,
                                                                    variable,method="_mint_", 
                                                                    colonne="pred")
       
      })
      
    })

    #save df in previous auc file

    auc <- readRDS(paste0(path,type,"/",variable,"/auc.rds"))[,1:13]
    auc<-left_join(auc,df,by=c("var","param"))
    saveRDS(auc,paste0(path,type,"/",variable,"/auc.rds"))
  }
}
T2<-Sys.time()
difftime(T2, T1) # 40s, bug pour d= 28, 26 et 26 (erdos, tree et cluster)

########
# récupérer les temps et itermax proprement
#######
for(type in c("erdos","tree","scale-free","cluster")){
  cparam<-switch(type,"erdos"=c("d","n","prob"),"tree"=c("d","n"),
                 "cluster"=c("d","n","dens","r"),"scale-free"=c("d","n"))
  # cparam<-"n"
  for(variable in cparam){
    length<-length(parameters[[variable]])
    df<-data.frame(var=rep(parameters[[variable]],Bgraph),param=rep(1:Bgraph,each=length),
                   itermax=rep(0,Bgraph*length), times=rep(0,Bgraph*length))
    
    sapply(parameters[[variable]],function(x){
      lapply(1:Bgraph,function(nbgraph){
        save_file<-paste0(path,type,"/",variable,"/Scores/")
        
        df$itermax[which(df$var==x & df$param==nbgraph)]<<-readRDS(paste0(save_file,"Graph",nbgraph,"_treeggm_",x,".rds"))[["itermax"]]
        df$times[which(df$var==x & df$param==nbgraph)]<<-readRDS(paste0(save_file,"Graph",nbgraph,"_treeggm_",x,".rds"))[["time"]]
      })
    })
    saveRDS(df,paste0(path,type,"/",variable,"/temps/timeMaxIter.rds"))
  }
}

visu_time_itermax<-function(type, variable){
  data<-readRDS(paste0(path,type,"/",variable,"/temps/timeMaxIter.rds"))
  p1<-ggplot(data,aes(var,times, group=var))+
    geom_boxplot()+
    labs(x=variable,y="times(s)")
  p2<-ggplot(data,aes(var,itermax, group=var))+
    geom_boxplot()+
    labs(x=variable)
  p3<-ggplot(data,aes(times))+
    geom_density()+
    labs(x="times(s)")
  p4<-ggplot(data,aes(itermax))+
    geom_density()
  grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2,top = textGrob(paste0(type," ",variable),gp=gpar(fontsize=15,font=6)))
}
type<-"cluster"
variable<-"d"
visu_time_itermax(type,variable)

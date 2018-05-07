library(tidyverse)

# récupérer la constante
root<-"/home/momal/Git/these/pack1/R/Simu/PLN/"
final<-matrix(NA,1,11)

for(type in c("tree","erdos","cluster","scale-free")){
   if (type=="tree" ) cparam="d" else cparam = c("d","prob")
  for(param in cparam){
    path<-paste0(root,type,"/",param,"/Sets_param/")
    frame<-matrix(NA,1,3)
    for(i in 1:length(dir(path))){
      file<-dir(path)[i]
      cste <- readRDS(paste0(path,file))$cste
      frame<-rbind(frame,c(cste,unlist(strsplit(substr(file,6,nchar(file)-4),split="_"))))
    }
    
    frame<-as_tibble(frame[-1,])
    colnames(frame)<-c("constant","Graph","val_param")
    frame[,1:3]<-apply(frame[,1:3],2,function(x) as.numeric(as.character(x)))
    attach(frame,warn.conflicts = FALSE)
    frame<-frame[order(Graph,val_param),]
  
  
  ##############
  # les AUC
    path<-paste0(root,type,"/",param,"/auc.rds")
    auc<-readRDS(path,file)
              
    auc<-auc%>%
      group_by(var,param)%>%
      summarise(treeggm=mean(treeggm),glasso=mean(glasso),ggm1step=mean(ggm1step))
    colnames(auc)[1:2]<-c("val_param","Graph")             
    attach(auc,warn.conflicts = FALSE)
    auc<-auc[order(Graph,val_param),]
    auc$Graph<-as.numeric(as.character(auc$Graph))
  ##############
  # les nbtriangles et nb_edges
  
    path<-paste0(root,type,"/",param,"/Graphs_characteristics/")
    triangles<-dir(path)[sapply(dir(path),function(x) grepl("triangle",x))]
    edges<-dir(path)[!sapply(dir(path),function(x) grepl("triangle",x))]
    
    frame2<-matrix(NA,1,3)
    colnames(frame2)<-c("nb_triangles","val_param","Graph")
    for(i in 1:length(triangles)){
      file<-triangles[i]
      tab <- cbind(readRDS(paste0(path,file)),substr(file,6,nchar(file)-16))
      colnames(tab)<-c("nb_triangles","val_param","Graph")
      frame2<-rbind(frame2,tab)
    }
    frame2<-as_tibble(frame2[-1,])
    
    frame3<-matrix(NA,1,4)
    colnames(frame3)<-c("nb_edges_pred","nb_edges_obs","val_param","Graph")
    for(i in 1:length(edges)){
      file<-edges[i]
      tab <-readRDS(paste0(path,file)) %>% 
        filter(B==1) %>% 
        dplyr::select(-B) %>% 
        cbind(substr(file,6,nchar(file)-18))
      colnames(tab)<-c("nb_edges_pred","nb_edges_obs","val_param","Graph")
      frame3<-rbind(frame3,tab)
    }
    frame3<-as_tibble(frame3[-1,])
    frame2$Graph<-as.numeric(as.character(frame2$Graph))
    frame3$Graph<-as.numeric(as.character(frame3$Graph))
  ##############
  # join final
    framejoin<-frame2 %>% 
      full_join(frame3,by=c("Graph","val_param")) %>% 
      full_join(auc,by=c("Graph","val_param")) %>% 
      full_join(frame,by=c("Graph","val_param")) %>% 
      as_tibble()
    framejoin$type<-type
    framejoin$param<-param
    colnames(final)<-colnames(framejoin)
    final<-as_tibble(rbind(final,framejoin))
  }
  
}
final<-final[-which(is.na(final$Graph)),]
saveRDS(final,paste0(root,"/extract_data.rds"))

data<-final %>% 
  gather(method,auc,c(treeggm,glasso,ggm1step))

ggplot(data,aes(x=nb_triangles,y=auc,colour = factor(type)))+
  geom_point()

ggplot(data,aes(x=nb_edges_obs,y=nb_edges_pred,colour = factor(type)))+
  geom_point()+
  geom_abline()

ggplot(data,aes(x=nb_triangles,y=constant,colour = factor(type)))+
  geom_point()






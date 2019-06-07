
type<-"cluster"
variable<-"dens"
nbgraph<-15
valeur<-0.35
path<-"/Users/raphaellemomal/simulations/Simu/PLNcov/"
####################
res<-from_stored_graphs(type, variable, nbgraph, valeur)

inf_spieresid<-readRDS(paste0(path,type,"/",variable,"/Scores/Graph",nbgraph,"_spiecResid",valeur,".rds"))
dataPLN<-res[[1]]
dataPLNcond<-res[[2]]
length<-length(which(upper.tri(dataPLN,diag=FALSE)))

dataggplot<-data.frame(SpiecResid=sort(inf_spieresid[upper.tri(inf_spieresid,diag=FALSE)]),
                       PLN=sort(dataPLN[upper.tri(dataPLN,diag=FALSE)]),
                       PLNcond=sort(dataPLNcond[upper.tri(dataPLNcond,diag=FALSE)]),
                       index=seq(1,length)
)
dataggplot<-gather(dataggplot,method,scores,-"index")

##############@ coude 
ggplot(dataggplot,aes(index,scores))+
  geom_point()+
  facet_grid(rows=vars(method))

# les 1000 plus forts sont-ils les mêmes ?

seuil<-100
compute_intersect<-function(seuil){
  seuil_spiec<-dataggplot[length-seuil,3]
  seuil_cond<-dataggplot[nrow(dataggplot)-seuil,3]
  seuil_PLN<-dataggplot[2*length-seuil,3]
  indices_cond<-which(as.matrix(dataPLNcond[upper.tri(dataPLNcond)])>seuil_cond)
  indices_PLN<-which(as.matrix(dataPLN[upper.tri(dataPLN)])>seuil_PLN)
  indices_spiec<-which(as.matrix(inf_spieresid[upper.tri(inf_spieresid)])>seuil_spiec)
  # print(list(indices_PLN,indices_cond))
  vec1<-length(intersect(indices_cond,indices_PLN))*100/length(union(indices_cond,indices_PLN))
  vec2<-length(intersect(indices_cond,indices_spiec))*100/length(union(indices_cond,indices_spiec))
  vec3<-length(intersect(indices_spiec,indices_PLN))*100/length(union(indices_spiec,indices_PLN))
  return(c(vec1,vec2,vec3))
}

res_seuil<-sapply(1:100,function(x) compute_intersect(x))
res_seuil<-data.frame(t(res_seuil))
colnames(res_seuil)<-c("PLN/cond","spiec/cond", "spiec/PLN" )
res_seuil$index<-1:nrow(res_seuil)
res_seuil<-gather(res_seuil,comparison,percent,-"index")

ggplot(res_seuil,aes(index,percent,color=comparison))+
  geom_point(size=0.6)+
  geom_line()+
  labs(x = "quantity of most probable edges",y="in common (%)",title=paste0(type," ",variable, " = ",valeur," du graph ", nbgraph))


path<- "~/these/pack1/R/Simu/PLN/"
cparam<-c("d","prob")
nbgraph<-100
n<-100
type<-"erdos"
parameters<-list(c(seq(10,30,2)),c(seq(10,30,10)),c(seq(0,1.5,0.2)),c(seq(0.5,1.5,0.5)/20),c(seq(1,50,10)))
names(parameters)<-c("d","n","u","prob","r")


ptitestimreg<-function(sigma,K){
  Y<-generator_PLN_nocov(param$sigma,n)[[1]]
  n<-nrow(Y)
  PLN = PLN(Y ~ 1)
  Sigma<-PLN$model_par$Sigma
        
  return(estim_reg(generator_PLN_nocov(Sigma,n)[[1]],K))
}
B<-2
for(variable in cparam){
  seq<-parameters[[variable]]
 for(nbgraph in 1:nbgraph){
    for(x in seq){
        print(paste0("variable ",variable,"// nbgraph ",nbgraph," // x ",x))
        param<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",x,".rds"))
        K<-param$omega
        obj<-lapply(1:B, function(x) ptitestimreg(sigma,K))
        estim_nb_edges<-do.call(rbind,lapply(obj, function(x){x[[1]]}))
        path2<-paste0(path,type,"/",variable,"/")
        record(estim_nb_edges,x,c("nb_pred","nb_obs"),paste0(path2,"Graphs_characteristics/Graph",nbgraph),B)
    }
  }
}

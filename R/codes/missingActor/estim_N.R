library(EMtree)
library(PLNmodels)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(mvtnorm)
library(useful)
library(MASS)
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(parallel)
library(sparsepca)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")


estim_N<-function(p,B=100,nSF=200,n,cores,r){
  cliques<-lapply(1:nSF, function(seed){
    cat(paste0("\n seed ",seed, " : "))
    set.seed(seed) ; type="scale-free" ; O=1:p ; plot=FALSE 
    # Data
    missing_data<-missing_from_scratch(n,p,r,type,plot)
    counts=missing_data$Y; ZH=missing_data$ZH ; sigmaO= missing_data$Sigma; 
    omega=missing_data$Omega; trueClique=missing_data$TC[[1]]; hidden=missing_data$H
    # Observed parameters
    PLNfit<-PLN(counts~1, control=list(trace=0))
    MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; theta=PLNfit$model_par$Theta ; 
    matcovar=matrix(1, n,1) ; sigma_obs=PLNfit$model_par$Sigma
    
    #1 missing actors
    t1<-Sys.time()
    cliques_spca <- boot_FitSparsePCA(scale(MO),B,r=r, cores=3, unique=FALSE)
    t2<-Sys.time()
    time_boots=difftime(t2, t1)
    cat(paste0("\nseed ", seed," in ",round(time_boots,3), attr(time_boots, "units"),"\n"))
    return(cliques_spca)
  })
  return(cliques)
}
List_cliques<-estim_N(p=14, n=200,r=1)
saveRDS(List_cliques, "/Users/raphaellemomal/these/R/codes/missingActor/SimResults/200_listcliques.rds")


List_nbocc<-lapply(List_cliques, function(graph){
  graph$nb_occ
})

data_nbocc=do.call(rbind, lapply(List_nbocc, function(graph_nbocc){
  f1<-sum(graph_nbocc == 1)
  f2<-sum(graph_nbocc==2)
  return(c(f1=f1, f2=f2, nbcliques=length(graph_nbocc)))
})) %>% as_tibble() %>% mutate(B=100, seed=1:200, estimN = B+f1^2 / (2*f2))

#--- traitement résultat

N19=data_nbocc$estimN[19]
q75 = quantile(data_nbocc$estimN,0.75)
q90 = quantile(data_nbocc$estimN,0.90)
data_nbocc %>% ggplot(aes(estimN))+geom_density(color="salmon2", fill="salmon2", alpha=0.2)+theme_light()+
  geom_segment(x=N19, y=0, xend=N19, yend=0.0012, color="skyblue4", lineend="round",size=1)+
  geom_segment(x=q75, y=0, xend=q75, yend=0.0052, color="salmon2", lineend="round", size=1)+
  geom_segment(x=q90, y=0, xend=q90, yend=0.0012, color="salmon2", lineend="round", size=1)


H=15
simus=list()
badseed=c()
simus=lapply(1:200, function(seed){
  file= paste0("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/15nodes_1r_400_filtreWg/SF_seed",seed,".rds")
  if(file.exists(file)){
    readRDS(file)
  }else{badseed<<-c(badseed,seed)}
})
misAct=build_misAct(simus, H=15)
misAct=misAct %>% mutate(influence=unlist(purrr::map(nH, function(x){
  if(x<5) res="minor"
  if(x>=5 & x<7) res="medium"
  if(x>=7) res="major"
  return(res)})))
Qual<-data.frame(auc=do.call(rbind, lapply(simus, function(seed){
  Pg=seed$VEM_1$Pg
  ome=seed$omega
  diag(ome)=0
  auc=round(auc(pred = Pg, label = ome),4)
  return(auc)
})),ppvh=do.call(rbind, lapply(simus, function(seed){
  Pg=seed$VEM_1$Pg ;ome=seed$omega ;diag(ome)=0 ;H=15
  ppvh=accppvtpr(Pg,ome,h=H,seuil=0.5)[5]
  return(ppvh)
})),tprh=do.call(rbind, lapply(simus, function(seed){
  Pg=seed$VEM_1$Pg ;ome=seed$omega ; diag(ome)=0 ; H=15
  ppvh=accppvtpr(Pg,ome,h=H,seuil=0.5)[8]
  return(ppvh)
})), seed = (1:200)) %>% as_tibble() %>%
  mutate( nH =misAct$nH, influence=misAct$influence, estimN=data_nbocc$estimN,
          nbcliques=data_nbocc$nbcliques)
#--estimN
Qual %>% filter(is.finite(estimN)) %>% 
  ggplot(aes(estimN, auc))+geom_point()+
  geom_hline(yintercept = 0.5, color="red")+geom_point()+ 
  facet_wrap(~influence, scales="free")+mytheme.dark("")
which(!is.finite(Qual$estimN))
data_nbocc[30,]
Qual %>% 
  ggplot(aes(x=estimN))+geom_density()+ facet_wrap(~influence)+mytheme.dark("")

#-- nbcliques

Qual %>% 
  ggplot(aes(nbcliques, auc))+geom_point()+geom_hline(yintercept = 0.5, color="red")+geom_point()+ facet_wrap(~influence)+mytheme.dark("")

Qual %>% 
  ggplot(aes(x=nbcliques))+geom_density()+ facet_wrap(~influence)+mytheme.dark("")


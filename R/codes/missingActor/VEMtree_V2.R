library(EMtree)
library(PLNmodels) 
library(ggraph)
library(tidygraph)
library(tidyverse)
library(mvtnorm)
library(useful) 
library(MASS)
library(parallel)
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(harrypotter)
library(sparsepca)
library(tictoc)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-V2.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/modif_pkg.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/VEM_tools.R")

#-- data simulation
set.seed(20)
n=200 ;p=14;r=1;type="scale-free";plot=TRUE
O=1:p
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y
sigmaO= missing_data$Sigma
upsilon=missing_data$Upsilon
trueClique=missing_data$TC
hidden=missing_data$H
PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M  
SO<-PLNfit$var_par$S  
sigma_obs=PLNfit$model_par$Sigma
plot(sigmaO ,sigma_obs)
abline(0,1)
G=(missing_data$G)
ggimage(G)
#-- normalize the PLN outputs
D=diag(sigma_obs)
matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
MO=MO*matsig
SO=SO*matsig^2
 
#-- VEM initialisation
cliques_spca<- boot_FitSparsePCA(scale(MO),B=100,r=1, cores=3)
init=initVEM(counts = counts,initviasigma=cliques_spca$cliqueList[[1]], cov2cor(sigma_obs),MO,r = 1) 
Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit

##-- single VEM
resVEM<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.1,
                 maxIter=30, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=FALSE)
plotVEM(resVEM$Pg,G,r=1,seuil=0.5)
resVEM$features
tail(resVEM$lowbound$J,1)
M=resVEM$M
plot(missing_data$UH, M[,15])
abline(0,1)
ggimage(resVEM$Pg)

##-- the big picture
get.ListVEM<-function(seed, eps=1e-3){
  set.seed(seed) ; p=14 ; r=1 ;n=200 # 2 faible influence
  type="scale-free" ; O=1:p ;H=(p+1):(p+r); plot=FALSE 
  # Data
  missing_data<-missing_from_scratch(n,p,r,type,plot)
  counts=missing_data$Y;  trueClique=missing_data$TC[[1]]; hidden=missing_data$H
  # Observed parameters
  PLNfit<-PLN(counts~1, control=list(trace=0))
  MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; sigma_obs=PLNfit$model_par$Sigma
  #-- normalize the PLN outputs
  D=diag(sigma_obs)
  matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
  MO=MO*matsig
  SO=SO*matsig^2
  R=cov2cor(sigma_obs)
  #------
  cliques_spca<- boot_FitSparsePCA(scale(MO),B=100,r=1, cores=3)
  ListVEM<-List.VEM(cliquesObj =cliques_spca, counts, R, MO,SO,r=1,alpha=0.1,
                    eps=eps,maxIter=100, cores=3 )
  res=ListVEM
  
  return(res)
}
get_data<-function(ListVEM,seed, r=1){
  set.seed(seed)
  p=14
  missing_data<-missing_from_scratch(n=n,p=p,r=r,type="scale-free",plot=FALSE)
  G=missing_data$G; q=ncol(G)
  hidden=missing_data$H
  badinit=which(lapply(ListVEM, function(vem) length(vem)) !=14)
  if(length(badinit)!=0) ListVEM=ListVEM[-badinit]
  J=do.call(rbind,lapply(ListVEM, function(vem){tail(vem$lowbound$J,1)}))
  Jcor_diff=do.call(rbind, lapply(ListVEM, function(vem){
    getJcor(vem,p)
  }))
  qual=do.call(rbind, lapply(ListVEM, function(vem){
    Pg=vem$Pg
    auc=round(auc(pred = Pg, label = G),4) 
    ppvh=accppvtpr(Pg,G,h=15,seuil=0.5)[5] 
    tprh=accppvtpr(Pg,G,h=15,seuil=0.5)[8]
    return(data.frame(auc=auc, ppvh=ppvh, tprh=tprh))
  }))
  pen=do.call(rbind, lapply(ListVEM, function(vem){
    Pg=vem$Pg; Wg=vem$Wg; S=vem$S ;H=(p+1):(p+r) 
    pen_T=-( 0.5*sum( Pg * log(Wg+(Wg==0)) ) - logSumTree(Wg)$det) 
    pen_UH=  0.5*sum(apply(S,2,function(x){ log(sum(x))}))  
    return(data.frame(penT=pen_T, penUH=pen_UH))
  }))
 
  data=data.frame(J,qual,pen,Jcor_diff) %>%  as_tibble() %>% 
    mutate(num=1:length(ListVEM))
  return(data)
}
tic()
ListVEM20<-get.ListVEM(20)
toc()
ListVEM19<-get.ListVEM(19)
toc()
data20=get_data(ListVEM20, 20)
data19=get_data(ListVEM19, 19)
data20=data20 %>% mutate(ICL=J-penT-penUH)
data20%>%  mutate(maxJcor=Jcor==max(Jcor, na.rm=TRUE)) %>% 
  gather(key, value, -diff,-Jcor,-J,-maxJcor,-num,-detEg,-ICL,-penUH,-penT,-delta) %>% 
  ggplot(aes(value,Jcor))+geom_vline(xintercept=0.5, color="gray", linetype="dashed")+
  geom_point(aes(size=ifelse(maxJcor,2.5,2), shape=ifelse(maxJcor,15,20)))+
  scale_shape_identity() +scale_size_identity() + facet_wrap(~key)+mytheme.dark("")



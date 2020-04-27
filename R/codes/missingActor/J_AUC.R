library(EMtree)
library(PLNmodels)
library(LITree)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(mvtnorm)
library(useful)
library(mclust)
library(MASS)
library(parallel)
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(harrypotter)
library(sparsepca)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")

#------ simu parameters
set.seed(3)
n=200 ;p=14;r=1;type="scale-free";plot=TRUE
O=1:p
 
#------ Data simulation
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y ;sigmaO= missing_data$Sigma;omega=missing_data$Omega
trueClique=missing_data$TC; hidden=missing_data$H

PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M;SO<-PLNfit$var_par$S  ;sigma_obs=PLNfit$model_par$Sigma
theta=PLNfit$model_par$Theta ;matcovar=matrix(1, n,1)
ome=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), hidden)]
diag(ome)=0

#------ cliques et VEMtee

cliques_spca <- boot_FitSparsePCA(scale(MO),100,r=1, cores=3)
 
ListVEM<-List.VEM(cliquesObj=cliques_spca, counts, sigma_obs,
                  MO,SO,r=1,eps=1e-3, maxIter=200, alpha = 0.1,cores=3,
                  nobeta = FALSE)
goodPrec=!do.call(rbind,lapply(ListVEM, function(x) x$max.prec))
J=do.call(rbind,lapply(ListVEM, function(vem){tail(vem$lowbound$J,1)}))
maxJ_good=which(J==max(J[J<min(J[!goodPrec])]))
AUC=do.call(rbind, lapply(ListVEM, function(vem){
  Pg=vem$Pg
  AUC=round(auc(pred = Pg, label = ome),4)
}))
ppvh=do.call(rbind, lapply(ListVEM, function(vem){
  Pg=vem$Pg
  ppvh=accppvtpr(Pg,ome,h=15,seuil=0.5)[5]
}))  
ICL<-do.call(rbind, lapply(ListVEM, function(vem){
  J=tail(vem$lowbound$J,1)
  Wg=vem$Wg
  Pg=Kirshner(Wg)
  pen_T=-( sum( Pg * log(Wg+(Wg==0)) ) - logSumTree(Wg)$det) 
  ICL = J-pen_T
  return(ICL)
}))
data= data.frame(goodPrec, J,ICL, AUC,ppvh) 
data%>% 
  group_by(goodPrec) %>%  summarise(maxJ=max(J), mean.auc=mean(auc),
                                    mean.ppvh=mean(ppvh),
                                    indexmaxJ = which.max(J),
                                    indexmaxICL = which.max(ICL))
VEM_1=ListVEM[[maxJ_good]]
maxJ_good=which(!goodPrec)[which.max(ICL[!goodPrec])]
data=data %>% mutate(penT = -ICL + J)

data %>% ggplot(aes(auc,J, color=goodPrec))+geom_point()+mytheme.dark("")
data %>% ggplot(aes(ppvh,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
data %>% ggplot(aes(ppvh,J, color=goodPrec))+geom_point()+mytheme.dark("")
data %>% ggplot(aes(ICL,J,color=auc>0.6, shape=goodPrec))+geom_point()+
  mytheme.dark("auc>0.6")+geom_abline()
data %>% ggplot(aes(ppvh,auc, color=goodPrec))+geom_point()+mytheme.dark("")
data %>% ggplot(aes(auc,(penT), color=goodPrec))+geom_point()+mytheme.dark("")

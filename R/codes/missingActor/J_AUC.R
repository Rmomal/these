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
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-exactDet.R")



J_AUC<-function(seed, p,r,alpha,eig.tol=1e-6,cliques_spca=NULL,B=100,cores=3,plot=FALSE,type="scale-free",n=200){
  #------ Data simulation
  set.seed(seed)
  O=1:p
  H=(p+1):(p+r)
  missing_data<-missing_from_scratch(n,p,r,type,plot)
  counts=missing_data$Y ; sigmaO= missing_data$Sigma ; omega=missing_data$Omega
  trueClique=missing_data$TC ; hidden=missing_data$H
  PLNfit<-PLN(counts~1, control=list(trace=0))
  MO<-PLNfit$var_par$M ; SO<-PLNfit$var_par$S ; sigma_obs=PLNfit$model_par$Sigma
  theta=PLNfit$model_par$Theta ; matcovar=matrix(1, n,1)
  ome=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), hidden)]
  diag(ome)=0
  
  #------ cliques et VEMtee
  if(is.null(cliques)){
     cliques_spca <- boot_FitSparsePCA(scale(MO),B=B,r=1, cores=3)
  }
 
  ListVEM_filtre<-List.VEM(cliquesObj=cliques_spca, counts, sigma_obs,
                    MO,SO,r=1,eps=1e-3, maxIter=200, alpha = alpha,cores=cores,
                    nobeta = FALSE, filterDiag = TRUE,filterWg=TRUE)
  ListVEM_nofiltre<-List.VEM(cliquesObj=cliques_spca, counts, sigma_obs,
                           MO,SO,r=1,eps=1e-3, maxIter=200, alpha =alpha,cores=cores,
                           nobeta = FALSE, filterDiag = FALSE,filterWg=FALSE)
  ListVEM=c(ListVEM_filtre,ListVEM_nofiltre)
  #------ shape des résultats
  filtre=rep(c(TRUE,FALSE), each=length(unique(cliques_spca$cliqueList)))
  goodPrec=!do.call(rbind,lapply(ListVEM, function(x) x$max.prec))
  J=do.call(rbind,lapply(ListVEM, function(vem){tail(vem$lowbound$J,1)}))
  nbocc=do.call(rbind, lapply(ListVEM, function(vem){ vem$nbocc}))
  AUC=do.call(rbind, lapply(ListVEM, function(vem){
    Pg=vem$Pg
    AUC=round(auc(pred = Pg, label = ome),4) }))
  ppvh=do.call(rbind, lapply(ListVEM, function(vem){
    Pg=vem$Pg
    ppvh=accppvtpr(Pg,ome,h=15,seuil=0.5)[5] }))  
  Icl<-do.call(rbind, lapply(ListVEM, function(vem){
    J=tail(vem$lowbound$J,1) ;Wg=vem$Wg ; Pg=vem$Pg
    pen_T=-( sum( Pg * log(Wg+(Wg==0)) ) - logSumTree(Wg)$det) 
    ICL = J-pen_T
    return(ICL) }))
  vBICs<-do.call(rbind, lapply(ListVEM, function(vem){
    J=tail(vem$lowbound$J,1)
   r=1;d=1 ; q=p+r
    nbparam<-p*(d) + (p*(p+1)/2 +r*p)+(q*(q-1)/2 - 1)
    nbparam_mixedmodel<-round(SumTree(vem$W),1)*(1+q^2)-1
    vbic0=J- nbparam*log(n)/2
    vbic1=J- nbparam*log(n)/2 -nbparam_mixedmodel
    vbic2=J- (nbparam+nbparam_mixedmodel)*log(n)/2
    return(c(vbic0=vbic0,vbic1=vbic1,vbic2=vbic2)) }))
  JPLN<-do.call(rbind, lapply(ListVEM, function(vem){
    EhZZ=t(vem$M[,O])%*%vem$M[,O] + diag(colSums(vem$S[,O]))
    sigTilde = (1/n)*EhZZ
    omega=vem$omega
    EsO=vem$Pg*vem$omega+diag(diag(vem$omega))
    EgOm = EsO[O,O] - matrix(EsO[O,H],p,r)%*%matrix(EsO[H,O],r,p)/EsO[H,H]
    EgOm = nearPD(EgOm, eig.tol=eig.tol)$mat
    JPLN_SigT = part_JPLN(sigTilde,EhZZ=EhZZ)
    JPLN_EgOm = part_JPLN(EgOm,EhZZ=EhZZ, var=FALSE)
    return(data.frame(JPLN_SigT=JPLN_SigT,JPLN_EgOm=JPLN_EgOm))}))
 
  data= data.frame(goodPrec, J,Icl, AUC,ppvh, vBICs,filtre=filtre, JPLN_SigT=JPLN$JPLN_SigT,JPLN_EgOm=JPLN$JPLN_EgOm,
                   index=rep(1:length(unique(cliques_spca$cliqueList)),2)) 
  return(data)
}

part_JPLN<-function(mat_var,EhZZ, var=TRUE){
  if(var){
    partJPLN=-n*0.5*(det.fractional(mat_var, log=TRUE)) - 0.5*sum(EhZZ*solve(mat_var))
  }else{# si on donne une matrice de précision
    partJPLN=n*0.5*(det.fractional(mat_var, log=TRUE)) - 0.5*sum(EhZZ*(mat_var))
  }
  return(partJPLN)
}
# maxJ_good=which(J==max(J[J<min(J[!goodPrec])]))
# maxJ=which.max(J)
# choice=which.max(ICL[goodPrec])
#ICL[choice]

#--- experiments
cliques_spca19=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/cliques_spca19_10000.rds")
cliques_spca1=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/cliques_spca1_10000.rds")
small_cliques1=list()
small_cliques1$cliqueList=cliques_spca1$cliqueList[1:10]
small_cliques1$nb_occ = cliques_spca1$nb_occ[1:10]
tic()
seed1<-J_AUC(seed = 1,p = 14,r=1, alpha=0.1,eig.tol=1e-6, cliques_spca=cliques_spca1) # 17 min
toc()
tic()
seed19<-J_AUC(seed = 19,p = 14,r=1, alpha=0.1,eig.tol=1e-6, cliques_spca=cliques_spca19) # 71 min
toc()


plot(seed19$JPLN_SigT,seed19$JPLN_EgOm) 


seed1$seed=1 ; seed19$seed=19
seed_filtre=rbind(seed1,seed19)
saveRDS(seed_filtre, "/Users/raphaellemomal/these/R/codes/missingActor/SimResults/seed_filtre_eigtol_1e-6.rds")
#--- plots
seed_filtre %>% ggplot(aes( AUC,J, color=(goodPrec)))+geom_point()+
  geom_vline(xintercept = 0.5, linetype="dashed", color="gray")+
  facet_grid(filtre~seed)+mytheme.dark("Sans précision \nmachine:")

seed_filtre %>% ggplot(aes( AUC,Delta, color=(goodPrec)))+geom_point()+
  geom_vline(xintercept = 0.5, linetype="dashed", color="gray")+
  facet_grid(filtre~seed)+mytheme.dark("Sans précision \nmachine:")

seed_filtre %>% ggplot(aes( ppvh,Delta, color=(goodPrec)))+geom_point()+
  geom_vline(xintercept = 0.5, linetype="dashed", color="gray")+
  facet_grid(filtre~seed)+mytheme.dark("Sans précision \nmachine:")


seed_filtre=seed_filtre %>% group_by(seed,filtre) %>% mutate(medDelta = quantile(Delta,0.1),
                                                 goodDelta = Delta<medDelta) %>% 
  ungroup()

seed_filtre %>% ggplot(aes( AUC,J, color=(goodDelta)))+geom_point()+
  geom_vline(xintercept = 0.5, linetype="dashed", color="gray")+
  facet_grid(filtre~seed) + mytheme.dark("Delta < q10(Delta)")
seed_filtre %>% dplyr::select(J,filtre,index,seed) %>% 
  spread(filtre,J) %>% 
  ggplot(aes(`FALSE`, `TRUE`))+geom_point()+facet_wrap(~seed)+geom_abline()+
  mytheme.dark("")

seed_filtre %>% dplyr::select(Delta,filtre,index,seed) %>% 
  spread(filtre,Delta) %>% 
  ggplot(aes(`FALSE`, `TRUE`))+geom_point()+facet_wrap(~seed)+geom_abline()+
  mytheme.dark("")+labs(title="Delta moins élevé avec filtres",x="sans filtre", y="avec filtre")

seed_filtre=seed_filtre %>% 
  mutate(tron_vbic1=ifelse(!is.finite(vbic1),.Machine$double.xmax,vbic1))
seed_filtre %>% ggplot(aes( AUC,tron_vbic1, color=(goodDelta)))+geom_point()+
  geom_vline(xintercept = 0.5, linetype="dashed", color="gray")+
  facet_grid(filtre~seed) + mytheme.dark("Delta < q10(Delta)")



data   %>% filter(goodPrec) %>%   filter(ICL==max(ICL))
data   %>%filter(J==max(J))
data%>% 
  group_by(goodPrec) %>%  summarise(maxJ=max(J), mean.auc=mean(AUC),
                                    mean.ppvh=mean(ppvh),
                                    indexmaxJ = which.max(J),
                                    indexmaxICL = which.max(ICL))
# VEM_1=ListVEM[[maxJ_good]]
# maxJ_good=which(!goodPrec)[which.max(ICL[!goodPrec])]
data %>%mutate(nbocc=nbocc) %>%  
  ggplot(aes(ppvh,J, color=as.factor(nbocc)))+geom_point()+mytheme.dark("")



data=data %>% mutate(penT = -ICL + J)
data05=data

data05 %>%ggplot(aes(AUC,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
data1 %>% ggplot(aes(AUC,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
data2 %>% ggplot(aes(AUC,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
data3 %>% ggplot(aes(AUC,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
data3 %>%mutate(nbocc=nbocc) %>%  
  ggplot(aes(ppvh,J, color=(nbocc > 1 & goodPrec)))+geom_point()+mytheme.dark("")

data %>% ggplot(aes(AUC,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
data %>% ggplot(aes(ppvh,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
data %>% ggplot(aes(ppvh,J, color=goodPrec))+geom_point()+mytheme.dark("")
data %>% ggplot(aes(ICL,J,color=AUC>0.6, shape=goodPrec))+geom_point()+
  mytheme.dark("auc>0.6")+geom_abline()
data %>% ggplot(aes(ppvh,AUC, color=goodPrec))+geom_point()+mytheme.dark("")
data %>% ggplot(aes(AUC,(penT), color=goodPrec))+geom_point()+mytheme.dark("")



#########@
# tests locaux

which(ICL>-2200)
trueClique
clique=cliques_spca$cliqueList[[1]]
clique=ListVEM[[33]]$clique
init=initVEM(counts = counts,initviasigma=clique, sigma_obs,r = 1)
Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit

VEM=VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, eps=1e-3, 
            alpha=0.3, 
            maxIter=200, plot=TRUE,print.hist=FALSE,filterWg=TRUE,
            verbatim = TRUE,nobeta = FALSE, filterDiag = TRUE)

vem=ListVEM[[27]]
plotVEM(vem$Pg,ome,r=1,seuil=0.5)
tail(vem$lowbound$J,1)
#-2739.929 avec alpha=0.1 et filtres, max.prec=FALSE
#-2746.588 avec alpha=0.1 sans filtres, max.prec=TRUE
#-3105.91 avec alpha=0.05, avec filtres, max.prec=FALSE (sans filtres passe pas)
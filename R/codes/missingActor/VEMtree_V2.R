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
set.seed(22) # n voisins de 35  73  89 120 162 206 342 367 395
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
init=initVEM(counts = counts,initviasigma=NULL, cov2cor(sigma_obs),MO,r = 0) 
Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit

##-- single VEM
#test blockmodels
resVEM0<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.15,
                 maxIter=30, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=FALSE)
ggimage(resVEM0$Pg)
draw_network(resVEM0$Pg, curv=0, layout="fr", nodes_label = 1:p, size=5)
sbm.vem <- BM_poisson("SBM_sym",resVEM0$Pg)
sbm.vem$estimate()
paramEstimSBMPoisson <- extractParamBM(sbm.vem,3)
c1=which(paramEstimSBMPoisson$Z==1)
c2=which(paramEstimSBMPoisson$Z==2)
c3=which(paramEstimSBMPoisson$Z==3)
toc()
EM=EMtree(cov2cor(sigma_obs))
Ghat=1*(EM$edges_prob>0.15)
sbm.seuil <- BM_poisson("SBM_sym",Ghat)
sbm.seuil$estimate()
sbm.prob <- BM_poisson("SBM_sym",EM$edges_prob)
sbm.prob$estimate()
paramEstimSBMPoisson_seuil <- extractParamBM(sbm.seuil,3)
paramEstimSBMPoisson_em <- extractParamBM(sbm.em,3)
# fin test blockmodels
# comparaison probabilitiés EMtree et VEMtree_r0
comp_em_vem=data.frame(EMtree=F_Sym2Vec(EM$edges_prob), VEMtree=F_Sym2Vec(resVEM0$Pg))
g1<-comp_em_vem %>% ggplot(aes(VEMtree, EMtree))+geom_point()+theme_light()
g2<-comp_em_vem %>% gather(key, value) %>% ggplot(aes(value, key, color=key, fill=key))+geom_density_ridges(alpha=0.5)+
  labs(x="edges probabilities",y="")+mytheme.dark("")+guides(color=FALSE, fill=FALSE)
grid.arrange(g1, g2, ncol=2)
# fin comparaison

# single VEM - bootstrap
cliques_spca<- boot_FitSparsePCA(scale(MO),B=100,r=1, cores=3)
# single VEM - plurisuers axes
cliques_spca<-FitSparsePCA(counts, r=2)$cliques
complement=lapply(cliques_spca, function(clique){setdiff(1:p,clique)})
clique=list()
clique$cliqueList=lapply(c(cliques_spca,complement), function(cl) list(cl))
ListVEM<-List.VEM(cliquesObj =clique, counts, cov2cor(sigma_obs), MO,SO,r=1,alpha=0.1,
                  eps=1e-3,maxIter=100, cores=3, trackJ=FALSE)
vecJ=do.call(rbind, lapply(ListVEM, function(vem){
  if(length(vem)==12){
    J=tail(vem$lowbound$J, 1)
  }else{
    J=NaN}
}))
VEM1=ListVEM[[which.max(vecJ)]]
0.5*sum(apply(VEM1[[1]]$S,2,function(x){ log(sum(x))}))+15*n*0.5*(1+log(2*pi))
0.5*sum(log(VEM1[[1]]$S))+15*n*0.5*(1+log(2*pi))
plotVEM(VEM1$Pg, missing_data$G,r=1, seuil=0.5)
plotVEM(ListVEM[[2]]$Pg, missing_data$G,r=1, seuil=0.5)
plotVEM(ListVEM[[3]]$Pg, missing_data$G,r=1, seuil=0.5)
plotVEM(ListVEM[[4]]$Pg, missing_data$G,r=1, seuil=0.5)

##-- single VEM trueClique
init=initVEM(counts = counts,initviasigma=trueClique, cov2cor(sigma_obs),MO,r = 1) 
Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
resVEM<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.1,
                 maxIter=30, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=TRUE)

 

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
  cliques_spca<- boot_FitSparsePCA((MO),B=100,r=1, cores=3)
  ListVEM<-List.VEM(cliquesObj =cliques_spca, counts, R, MO,SO,r=1,alpha=0.1,
                    eps=eps,maxIter=100, cores=3,trackJ=FALSE )
  res=ListVEM
  
  return(res)
}
get_data<-function(ListVEM,seed, r=1){
  set.seed(seed) ; n=200
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
tic()
ListVEM1<-get.ListVEM(10)
toc()
tic()
ListVEM7<-get.ListVEM(7)
toc()
data20=get_data(ListVEM20, 20)
data19=get_data(ListVEM19, 19)
data1=get_data(ListVEM1,10)
data7=get_data(ListVEM7,7)
data7=data7 %>% mutate(Jcor2=J-penT*n)
data7%>%  mutate(maxJ=J==max(J, na.rm=TRUE)) %>% 
  gather(key, value, -Jcor_diff,-J,-maxJ,-num,-penUH,-penT) %>% 
  ggplot(aes(value,J))+geom_vline(xintercept=0.5, color="gray", linetype="dashed")+
  geom_point(aes(size=ifelse(maxJ,2.5,2), shape=ifelse(maxJ,15,20)))+
  scale_shape_identity() +scale_size_identity() + facet_wrap(~key)+mytheme.dark("")

data1%>% 
  ggplot(aes(tprh,Jcor2))+geom_vline(xintercept=0.5, color="gray", linetype="dashed")+
  geom_point()+
  scale_shape_identity() +scale_size_identity() +mytheme.dark("")

#####
# example to compare with nestor outputs
set.seed(7) ; n= 200 ; p=15
data=missing_from_scratch(n=n,p=p,r=1,type="scale-free", plot=TRUE)
counts=data$Y
sigmaO= data$Sigma
upsilon=data$Upsilon
trueClique=data$TC
hidden=data$H
PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M  
SO<-PLNfit$var_par$S  
sigma_obs=PLNfit$model_par$Sigma
plot(sigmaO ,sigma_obs)
abline(0,1)
G=(data$G)
ggimage(G)
#-- normalize the PLN outputs
D=diag(sigma_obs)
matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
MO=MO*matsig
SO=SO*matsig^2

#-- find initial clique
findclique= boot_FitSparsePCA(scale(MO),B=100,r=1, cores=3)
initClique=findclique$cliqueList

fitList=List.VEM(cliquesObj=findclique, counts, sigma_obs, MO,SO, r=1,alpha=0.1, cores=1,maxIter=30,
                 eps=1e-4, nobeta=FALSE,  trackJ=FALSE,save=FALSE)
do.call(rbind,lapply(fitList, function(vem){
 c(auc=auc(vem$Pg, data$G), J=tail(vem$lowbound$J,1))
})) %>% as_tibble() %>% ggplot(aes(J,auc))+geom_point()+theme_light()

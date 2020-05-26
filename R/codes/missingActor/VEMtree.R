# etude de la quantité de bruit à privilégier:
# B=30 ; coeff=seq(0.1,4,0.5)
# test=lapply(coeff,function(c){
#   res=cbind( data.frame( t(sapply(1:B, function(x){
#     init.mclust(sigma_obs,title="Sigma",
#                 trueClique = trueClique,n.noise=round(p*c))$FPN
#   }))),c)
# })
# test=do.call(rbind,test)
# colnames(test)=c("FN","FP","coeff")
# test %>% mutate(nul = (FN==1 & FP==0)) %>% filter(!nul) %>% dplyr::select(-nul) %>% 
#   gather(FNP,value,-coeff) %>% 
#   ggplot(aes(as.factor(coeff),value,fill=FNP))+geom_boxplot()+mytheme.light

# test %>% mutate(nul = (FN==1 & FP==0)) %>% filter(!nul) %>% dplyr::select(-nul) %>%
#   group_by(coeff) %>% summarise(sumFP=sum(FP), sumFN=sum(FN), sum=sum(FP+FN))

# test epsilon SNR
# seqEpsi=seq(1.01, 5, 0.05)
# varyEpsi=sapply(seqEpsi, function(epsilon){
#   initial.param<-initEM(sigma_obs,n=n,cliqueList = list(initviasigma),cst=epsilon, pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
#   omega_init=initial.param$K0
#   compute_nSNR(K = omega_init,indexmissing = ncol(omega_init))
# })
# data.frame(SNR=varyEpsi,epsilon=seqEpsi) %>% ggplot(aes(epsilon, SNR))+geom_point()+mytheme

# VEMtree
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
#source("/home/mmip/Raphaelle/these/R/codes/missingActor/fonctions-missing.R")

# simu parameters
set.seed(20)
n=200 ;p=14;r=1;type="scale-free";plot=TRUE
O=1:p
################
#----- DATA
# simulate graph and omega, then sigma0 and finally counts
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y
sigmaO= missing_data$Sigma
omega=missing_data$Omega
upsilon=missing_data$Upsilon

trueClique=missing_data$TC
hidden=missing_data$H

PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M  
SO<-PLNfit$var_par$S  
sigma_obs=PLNfit$model_par$Sigma
theta=PLNfit$model_par$Theta 
matcovar=matrix(1, n,1)
order_sigma<-sigma_obs[c(trueClique[[1]],setdiff(1:p,trueClique[[1]])),
                       c(trueClique[[1]],setdiff(1:p,trueClique[[1]]))]
ggimage(order_sigma)
ggimage(solve(order_sigma))
plot(sigmaO ,sigma_obs)
ome_init=upsilon[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), hidden)]
ome=ome_init ; diag(ome)=0
ggimage(ome)

get.ListVEM<-function(seed, eps=1e-3){
  set.seed(seed) ; p=14 ; r=1 ;n=200 # 2 faible influence
  type="scale-free" ; O=1:p ;H=(p+1):(p+r); plot=FALSE 
  # Data
  missing_data<-missing_from_scratch(n,p,r,type,plot)
  counts=missing_data$Y; ZH=missing_data$ZH ; sigmaO= missing_data$Sigma; 
  omega=missing_data$Omega; trueClique=missing_data$TC[[1]]; hidden=missing_data$H
  # Observed parameters
  PLNfit<-PLN(counts~1, control=list(trace=0))
  MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; theta=PLNfit$model_par$Theta ; 
  matcovar=matrix(1, n,1) ; sigma_obs=PLNfit$model_par$Sigma
  #------
  cliques_spca<- boot_FitSparsePCA(scale(MO),B=100,r=1, cores=3)
  ListVEM<-List.VEM(cliquesObj =cliques_spca, counts, sigma_obs, MO,SO,r=1,alpha=0.1,
                           eps=eps,maxIter=100, nobeta=FALSE, cores=3,updateSH=TRUE,
                           filterDiag = FALSE, filterWg=FALSE,save=FALSE)
  res=ListVEM
  
 return(res)
}
ListVEM30_1<-get.ListVEM(30)
ListVEM30_2<-get.ListVEM(30, eps=1e-2)
ListVEM33_1<-get.ListVEM(33)
ListVEM33_2<-get.ListVEM(33, eps=1e-2)
ListVEM20_1<-get.ListVEM(20)
ListVEM20_2<-get.ListVEM(20, eps=1e-2)
ListVEM31_1<-get.ListVEM(31)
ListVEM31_2<-get.ListVEM(31, eps=1e-2)

ListVEM30<-get.ListVEM(30)
ListVEM33<-get.ListVEM(33)
ListVEM20<-get.ListVEM(20)
ListVEM31<-get.ListVEM(31)

ListVEM19<-get.ListVEM(19,eps = 1e-3)
do.call(rbind, lapply(ListVEM19, function(vem){length(vem)}))
get_data<-function(ListVEM,seed){
  set.seed(seed)
  missing_data<-missing_from_scratch(n=n,p=14,r=1,type="scale-free",plot=FALSE)
  omega=missing_data$Omega; 
  hidden=missing_data$H
  sorted_omega=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden),hidden)]
  diag(sorted_omega)=0
 badinit=which(lapply(ListVEM, function(vem) length(vem)) !=15)
 ListVEM=ListVEM[-badinit]
  J=do.call(rbind,lapply(ListVEM, function(vem){tail(vem$lowbound$J,1)}))
  auc=do.call(rbind, lapply(ListVEM, function(vem){
    Pg=vem$Pg
    auc=round(auc(pred = Pg, label = sorted_omega),4)
  }))
  ppvh=do.call(rbind, lapply(ListVEM, function(vem){
    Pg=vem$Pg
    ppvh=accppvtpr(Pg,sorted_omega,h=15,seuil=0.5)[5]
  }))
  tprh=do.call(rbind, lapply(ListVEM, function(vem){
    Pg=vem$Pg  
    tprh=accppvtpr(Pg,sorted_omega,h=15,seuil=0.5)[8]
  }))
  Jcor_diff=do.call(rbind, lapply(ListVEM, function(vem){
    getJcor(vem,14)
  }))
  sumP=do.call(rbind, lapply(ListVEM, function(vem){
   abs(sum(vem$Pg)-28) 
  }))
  projL=do.call(rbind,lapply(ListVEM, function(vem) vem$projL))
  nvois=do.call(rbind,lapply(ListVEM, function(vem){
    sum(vem$Pg[,15]>0.5)
  }))
  data=data.frame(J,auc,ppvh,tprh,Jcor_diff, sumP, projL, nvois) %>%  as_tibble() %>% 
    mutate(num=1:length(ListVEM))
  return(data)
}

data30_1<-get_data(ListVEM30_1,30)
data30_2<-get_data(ListVEM30_2,30)
data33_1<-get_data(ListVEM33_1,33)
data33_2<-get_data(ListVEM33_2,33)
data20_1<-get_data(ListVEM20_1,20)
data20_2<-get_data(ListVEM20_2,20)
data31_1<-get_data(ListVEM31_1,31)
data31_2<-get_data(ListVEM31_2,31)


data20<-get_data(ListVEM20[[1]],20)
data30<-get_data(ListVEM30[[1]],30)
data31<-get_data(ListVEM31[[1]],31)
data33<-get_data(ListVEM33[[1]],33)
data20_SH<-get_data(ListVEM20,20)
data30_SH<-get_data(ListVEM30,30)
data31_SH<-get_data(ListVEM31,31)
data33_SH<-get_data(ListVEM33,33)
data19=get_data(ListVEM19,19)

data20_SH=data20_SH %>% mutate(maxJcor=Jcor==max(Jcor[sumP<1e-1], na.rm=TRUE))
data20_SH %>% 
  gather(key, value, -sumP,-diff,-Jcor,-J,-maxJcor,-num,-detEg,-projL,-nvois) %>% 
  ggplot(aes(value,Jcor, color=(sumP>1e-1)))+geom_vline(xintercept=0.5, color="gray", linetype="dashed")+
  geom_point(aes(size=ifelse(maxJcor,2.5,2), shape=ifelse(maxJcor,15,20)))+
  scale_shape_identity() +scale_size_identity() + facet_wrap(~key)+mytheme.dark("")
data20_SH %>% filter(tprh>0.5)

ListVEM31_t=ListVEM31[-which(lapply(ListVEM31, function(vem) length(vem))!=15)]
plotVEM(ListVEM20_t[[31]]$Pg, ome, r=1, seuil=0.5)
ggimage(ListVEM20_t[[31]]$Wg)
hist(log(ListVEM20_t[[31]]$Wg[O,O]))
hist(log(ListVEM20_t[[31]]$Wg[O,H]))
data31_SH %>% mutate(maxJcor=ifelse(is.na(Jcor), FALSE,Jcor==max(Jcor, na.rm=TRUE))) %>% 
  gather(key, value, -sumP,-diff,-Jcor,-J,-maxJcor,-num,-detEg,-projL,-nvois) %>% 
  ggplot(aes(value,Jcor, color=sumP))+geom_vline(xintercept=0.5, color="gray", linetype="dashed")+
  geom_point(aes(size=ifelse(maxJcor,2.5,2), shape=ifelse(maxJcor,15,20)))+
  scale_shape_identity() +scale_size_identity() + facet_wrap(~key)+mytheme.dark("")
lapply(ListVEM20[[1]],function(vem) vem$time)

vem1 =ListVEM20_F[[24]]
vem2 =ListVEM20_F_SH[[12]]
getJcor(vem1,14, 1e-14, 1e-14)
getJcor(vem2,14, 1e-14, 1e-14)
tail(vem1$lowbound$J,1)


vem2=ListVEM_boot[[18]]
vem_s = ListVEM_boot[[25]]
tail(vem1$lowbound$J,1)
tail(vem2$lowbound$J,1)
plot(vem1$Pg, vem2$Pg)
table(F_Sym2Vec(vem1$Pg)>0.5,F_Sym2Vec(vem2$Pg)>0.5)
plot(vem1$omega, vem2$omega)
plot(vem1$M[,15], vem2$M[,15])
plot(log(vem1$W+(vem1$W==0)), log(vem2$W+(vem2$W==0)))
plot(log(vem1$Wg+(vem1$Wg==0)), log(vem2$Wg+(vem2$Wg==0)))
logSumTree(vem1$Wg)
logSumTree(vem2$Wg)
plot(vem1$Pg*log(vem1$Wg+(vem1$Wg==0)),vem2$Pg*log(vem2$Wg+(vem2$Wg==0)))
# entropie des arbres
-(0.5*sum(vem1$Pg*log(vem1$Wg+(vem1$Wg==0))) - logSumTree(vem1$Wg)$det)
-(0.5*sum(vem2$Pg*log(vem2$Wg+(vem2$Wg==0))) - logSumTree(vem2$Wg)$det)
# terme de trace dans Z|T
(- 0.5)* sum( ((vem$Pg+diag(15))*vem$omega)*( t(vem$M)%*%vem$M + diag(colSums(vem$S)) ) ) 
L1=LowerBound(vem1$Pg ,vem1$omega, vem1$M, vem1$S, vem1$W, vem1$Wg,vem1$p,
           logSumTree(vem1$W)$det, logSumTree(vem1$Wg)$det)
L2=LowerBound(vem2$Pg ,vem2$omega, vem2$M, vem2$S, vem2$W, vem2$Wg,vem2$p,
           logSumTree(vem2$W)$det, logSumTree(vem2$Wg)$det)
L1-L2


vem1$lowbound  %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>%
    ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(parameter)), size=3)+geom_line()+
    facet_wrap(~key, scales="free")+
    labs(x="iteration",y="", title="Lower bound and components")+mytheme+
    scale_color_discrete("")


vem2$features
 
vem1 =ListVEM20_F[[10]]
testclique=vem1$clique
cliques=boot_FitSparsePCA(counts,B=100,r = 1,cores=3)
init=initVEM(counts = counts,initviasigma=trueClique, sigma_obs,r = 1) #cliques_spca$cliqueList[[4]]
Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit

test_SH=VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, eps=1e-3, 
        alpha=0.1, 
        maxIter=40, plot=TRUE,print.hist=FALSE,filterWg=FALSE,updateSH = TRUE,
        verbatim = TRUE,nobeta = FALSE, filterDiag = FALSE)
plotVEM(test_SH$Pg,ome,1,0.5)
sum(test_SH$Pg)
test_SH$features
tail(test_SH$lowbound$J,1)
getJcor(test_SH,14,1e-8,1e-7)
colSums(test_SH$S)
diag(test_SH$omega)
testclique2=ListVEM_boot[[21]]$clique
init=initVEM(counts = counts,initviasigma=testclique, sigma_obs,r = 1) #cliques_spca$cliqueList[[4]]
Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit

bontest=VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, eps=1e-3, 
             alpha=0.1, 
             maxIter=50, plot=TRUE,print.hist=FALSE,filterWg=FALSE,
             verbatim = TRUE,nobeta = FALSE, filterDiag = FALSE)
plotVEM(bontest$Pg,ome,1,0.5)
bontest$features
tail(bontest$lowbound$J,1)
getJcor(bontest,14)
# data%>% 
#   group_by(goodPrec) %>%  summarise(maxJ=max(J), mean.auc=mean(auc),
#                                     mean.ppvh=mean(ppvh),
#                                     indexmaxJ = which.max(J),
#                                     indexmaxICL = which.max(ICL))
# VEM_1=ListVEM[[maxJ_good]]
# maxJ_good=which(!goodPrec)[which.max(ICL[!goodPrec])]
# data=data %>% mutate(penT = -ICL + J)
# # data1=data
#  data1=data
# data %>% ggplot(aes(auc,J, color=goodPrec))+geom_point()+mytheme.dark("")
# data %>% ggplot(aes(ppvh,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
# data %>% ggplot(aes(ppvh,J, color=goodPrec))+geom_point()+mytheme.dark("")
# data %>% ggplot(aes(ICL,J,color=auc>0.6, shape=goodPrec))+geom_point(size=3)+
#   mytheme.dark("auc>0.6")+geom_abline()
# data19 %>% ggplot(aes(ppvh,auc, color=goodPrec))+geom_point()+mytheme.dark("")
# data19 %>% ggplot(aes(auc,(penT), color=goodPrec))+geom_point()+mytheme.dark("")
# 
# 
# data1 %>% ggplot(aes(auc,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
# data1 %>% ggplot(aes(ppvh,ICL, color=goodPrec))+geom_point()+mytheme.dark("")
# data1 %>% ggplot(aes(auc,J, color=goodPrec))+geom_point()+mytheme.dark("")
# data1 %>% ggplot(aes(ICL,J,color=auc>0.6, shape=goodPrec))+geom_point(size=3)+
#   mytheme.dark("auc>0.6")+geom_abline()
# data1 %>% ggplot(aes(ppvh,auc, color=goodPrec))+geom_point()+mytheme.dark("")
#   data1 %>% ggplot(aes(auc,(penT), color=goodPrec))+geom_point()+mytheme.dark("")
#   
  # sélectionner le premier modèle sans problème nimériques dont la J
  # est inférieure au plus petit J des modèles qui ont eu un pb numérique ?
  
# vBICs<-(criteria(ListVEM,counts,theta, matcovar,r))
# vBICs$J
# hist(vBICs$J, breaks=10)
# best=which.max(vBICs$J)
# VEM_1<-ListVEM[[best]]
# computeFPN(VEM_1$clique,trueClique = trueClique[[1]],p=p) 
# VEM_1$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>%
#   ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(parameter)), size=3)+geom_line()+
#   facet_wrap(~key, scales="free")+
#   labs(x="iteration",y="", title="Lower bound and components")+mytheme+
#   scale_color_discrete("")
#  plotVEM(VEM_1$Pg,ome,r=1, 0.5)
# VEM_1$clique
# trueClique
# TJ=True_lowBound(Y=counts, M=VEM_1$M,VEM_1$S,theta=theta,X=matcovar,
#               W=VEM_1$W,Wg=VEM_1$Wg,Pg=VEM_1$Pg, omega=VEM_1$omega)
# ICL_T(TJ, Pg=VEM_1$Pg, Wg=VEM_1$Wg)
# crit=criteria(ListVEM,counts=counts,theta = theta,matcovar = matcovar, r=r)
# crit%>%
#   gather(key, value) %>%
#   ggplot(aes(key, value, color=key, fill=key))+geom_boxplot(alpha=0.3)+
#   mytheme.dark("")


####################
#-----  VEM
# find alpha on the observed part of the initial "non beta" quantities needed to compute the beta tilde
# alpha<-computeAlpha(omegainit[O,O], MO, SO)
# alpha=1
cliques_spca<- boot_FitSparsePCA(scale(MO),B=100,r=1, cores=3)
init=initVEM(counts = counts,initviasigma=cliques_spca$cliqueList[[5]], sigma_obs,r = 1) #cliques_spca$cliqueList[[4]]
Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
test=omegainit ;diag(test)=0

# #-- alpha
# D=.Machine$double.xmax
# f<-function(u,q){u*exp(u)-D^(1/(q-1))/sqrt(q*(q-1))} # f est monotone croissante
# fprim=function(u,q){
#   exp(u)*(u-1)-u*D^(1/(q-1))/sqrt(q*(q-1))
# }
# minimum=optimize(fprim, c(0,100),maximum = FALSE, q=16) # max 3 missing actors
# alpha_sup=minimum$minimum/n

#-- VEM
resVEM<- VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, eps=1e-3, 
                alpha=0.1, 
                maxIter=50, plot=TRUE,print.hist=FALSE,filterWg=FALSE,updateSH = TRUE,
                verbatim = TRUE,nobeta = FALSE, filterDiag = FALSE)#},
                #error=function(e){e},finally={})
sum(resVEM$Pg)
-(0.5*sum(resVEM$Pg*log(resVEM$Wg+(resVEM$Wg==0)))-logSumTree(resVEM$Wg)$det)
tail(resVEM$lowbound$J,1)
resVEM$features
plotVEM(resVEM$Pg,ome,r=1,seuil=0.5)
 plot(missing_data$ZH,resVEM$M[,15],pch=20)
rvalue=summary(lm(missing_data$ZH~resVEM$M[,15]))$r.squared

values=courbes_seuil(probs = resVEM$Pg,omega = ome,
                        h = 15,seq_seuil = seq(0,1,0.02))
values %>% mutate(crit = PPV+TPR) %>% filter(crit==max(crit, na.rm=TRUE)) %>%
  summarise(mins=min(seuil), maxs=max(seuil), PPV=max(PPV), PPVO=max(PPVO),PPVH=max(PPVH), 
            TPR=max(TPR),TPRO=max(TPRO),TPRH=max(TPRH))
plotVerdict(values, seuil)+guides(color=FALSE)
 
#ggsave("precrec_missing.png", plot=p, width=8, height=4,path= "/Users/raphaellemomal/these/R/images")

# lower bound check
ListVEM[[1]]$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-V6) %>% 
  ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(V6)), size=3)+geom_line()+
  facet_wrap(~key, scales="free")+
  labs(x="iteration",y="", title="Lower bound and components")+mytheme+
  scale_color_discrete("")


#====
#ajustement omega

# vrai vs initialisation vs estimation
# ome_init = vrai omega reordonne pour comparaison (acteur manquant en fin de matrice)
# omega_init = matrice d'initialisation obtenue par initEM et mclust
# omega_res = omega estime par VEMtree
omegas = data.frame(vrai = F_Sym2Vec(ome_init) , init = F_Sym2Vec(omega_init),
                    estimation = F_Sym2Vec(resVEM$omega))
omegas %>% gather(key, value, -vrai) %>% 
  ggplot(aes(as.factor(vrai), value,  fill=key, clor=key))+geom_boxplot()+
  mytheme#+coord_cartesian(ylim = c(-1,2))

omegas %>% group_by(vrai) %>% summarise(q75=quantile(estimation, 0.75),
                                        q25=quantile(estimation, 0.25))
quantile(omegas$estimation[omegas$vrai==1], 0.25) - quantile(omegas$estimation[omegas$vrai==0], 0.75)

Diagomegas = data.frame(vrai = diag(ome_init) , init = diag(omega_init),
                        estimation = diag(resVEM$omega))
Diagomegas[-nrow(Diagomegas),] %>% gather(key, value, -vrai) %>% 
  ggplot(aes((vrai), value,  color=key))+geom_point()+ theme_light()+
  geom_abline()

rvalue=summary(lm(Diagomegas$estimation~Diagomegas$vrai))$r.squared

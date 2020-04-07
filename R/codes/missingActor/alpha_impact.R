library(EMtree)
library(PLNmodels)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(mvtnorm)
library(useful)
library(mclust)
library(MASS)
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(kableExtra)
library(parallel)
library(sparsepca)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")

##### function



set.seed(seed)
type="scale-free" ; O=1:p ; plot=FALSE 
# Data
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y; ZH=missing_data$ZH ; sigmaO= missing_data$Sigma; 
omega=missing_data$Omega; trueClique=missing_data$TC[[1]]; hidden=missing_data$H
# Observed parameters
PLNfit<-PLN(counts~1, control=list(trace=0))
MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; theta=PLNfit$model_par$Theta ; 
matcovar=matrix(1, n,1) ; sigma_obs=PLNfit$model_par$Sigma
if(r!=0){
  sorted_omega=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), 
                                                           hidden)]
}else{
  sorted_omega=omega
}
#cliques for initialisations r=1 and r=2
t1<-Sys.time()
cliques_spca1 <- boot_FitSparsePCA(scale(counts),B,r=1)
cliques_spca2 <- boot_FitSparsePCA(scale(counts),B,r=2)
t2<-Sys.time()
time_boots=difftime(t2, t1)

eps=1e-3
init0=initVEM(counts = counts, initviasigma = NULL,  sigma_obs,r = 0)
Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 

# alpha=1
alpha=1
VEM_01<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
                Wg_init =Wginit,plot = plot, maxIter = maxIter,print.hist = FALSE,
                vraiOm = NULL, alpha=alpha, verbatim=FALSE, filterPg = TRUE, filterWg = TRUE )
ListVEM11<-List.VEM(cliqueList=cliques_spca1, counts, sigma_obs, MO,SO,alpha,r=1, cores=3,
                    eps=eps,maxIter)
ListVEM21<-List.VEM(cliqueList=cliques_spca2, counts, sigma_obs, MO,SO,alpha,r=2, cores=3,
                    maxIter, eps=eps)
# alpha=1/n
alpha=1/n
VEM_02<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
                Wg_init =Wginit,plot = plot, maxIter = maxIter,print.hist = FALSE,
                vraiOm = NULL, alpha=alpha, verbatim=FALSE, filterPg = TRUE, filterWg = TRUE )
ListVEM12<-List.VEM(cliqueList=cliques_spca1, counts, sigma_obs, MO,SO,alpha,r=1, cores=3,
                    eps=eps,maxIter)
ListVEM22<-List.VEM(cliqueList=cliques_spca2, counts, sigma_obs, MO,SO,alpha,r=2, cores=3,
                    maxIter, eps=eps)
# alpha tuned
alpha<-tryCatch(expr={computeAlpha(omegainit,default =1/n, MO, SO, plot=plot)},
                error=function(e){message("sythetic alpha")
                  return(1/n)}) #0.46
VEM_03<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
                Wg_init =Wginit,plot = plot, maxIter = maxIter,print.hist = FALSE,
                vraiOm = NULL, alpha=alpha, verbatim=FALSE, filterPg = TRUE, filterWg = TRUE )
ListVEM13<-List.VEM(cliqueList=cliques_spca1, counts, sigma_obs, MO,SO,alpha,r=1, cores=3,
                    eps=eps,maxIter)
ListVEM23<-List.VEM(cliqueList=cliques_spca2, counts, sigma_obs, MO,SO,alpha,r=2, cores=3,
                    maxIter, eps=eps)
#--
crit0<-criteria(list(VEM_01),counts,theta, matcovar,r=0)
crit1<-criteria(ListVEM11,counts,theta, matcovar,r=1)
crit2<-criteria(ListVEM21,counts,theta, matcovar,r=2)
crit_alpha1=rbind(crit0, crit1,crit2)
crit_alpha1$alpha="1"
crit0<-criteria(list(VEM_02),counts,theta, matcovar,r=0)
crit1<-criteria(ListVEM12,counts,theta, matcovar,r=1)
crit2<-criteria(ListVEM22,counts,theta, matcovar,r=2)
crit_alpha2=rbind(crit0, crit1,crit2)
crit_alpha2$alpha="1/n"
crit0<-criteria(list(VEM_03),counts,theta, matcovar,r=0)
crit1<-criteria(ListVEM13,counts,theta, matcovar,r=1)
crit2<-criteria(ListVEM23,counts,theta, matcovar,r=2)
crit_alpha3=rbind(crit0, crit1,crit2)
crit_alpha3$alpha="tuned"
############
crit=rbind(crit_alpha1, crit_alpha2, crit_alpha3)
plot=crit %>% gather(key, value, -r, -alpha) %>%
  ggplot(aes(as.factor(r), value, color=as.factor(alpha), fill=as.factor(alpha)))+geom_boxplot(alpha=0.3)+
  facet_wrap(~as.factor(key))+mytheme.dark("alpha=")+
  labs(x="number of missing actors",y="values",title="r=0 : impact of alpha on the lower bound (and ICL, vBIC)")
ggsave(filename = "alpha_impacts_J.png", plot = plot,
       path ="/Users/raphaellemomal/these/R/images/", width = 9, height = 4)
############

crit_alpha1 %>% as_tibble() %>% 
  mutate(penICL = J-ICL) %>% group_by(r) %>% 
  summarise(mvBIC=median(vBIC),mICL=median(ICL),mJ=median(J), mpenICL=median(penICL))
crit_alpha1 %>% gather(key, value, -r, -alpha) %>%
  ggplot(aes(as.factor(r), value, color=as.factor(key), fill=as.factor(key)))+geom_boxplot(alpha=0.3)+
  facet_wrap(~as.factor(key))+mytheme.dark("")+
  labs(x="number of missing actors",y="values",title=)

 
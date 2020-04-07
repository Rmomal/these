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


seed=1
set.seed(seed)
p=14 ; B=100;
type="scale-free" ; O=1:p ; plot=FALSE 
# Data
r=0 ; n=200 ; maxIter=100
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
 ome=sorted_omega ; diag(ome)=0
eps=1e-3
init0=initVEM(counts = counts, initviasigma = NULL,  sigma_obs,r = 0)
Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 

# alpha=1
alpha=1
VEM_0<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
                Wg_init =Wginit,plot = TRUE, maxIter = maxIter,print.hist = FALSE,
                vraiOm = NULL, alpha=alpha, verbatim=FALSE, filterPg = TRUE, filterWg = TRUE )
ggimage(VEM_0$Pg) ; ggimage(ome)
VEMr<-lapply(1:6, function(r){
  cat(paste0(r," missing actors: "))
  t1<-Sys.time()
  cliques_spca <- boot_FitSparsePCA(scale(counts),B,r=r)
  ListVEM<-List.VEM(cliqueList=cliques_spca, counts, sigma_obs, MO,SO,alpha,r=r, cores=1,
                    eps=eps,maxIter)
  t2<-Sys.time()
  runtime=difftime(t2,t1)
  cat(paste0(round(runtime,3), attr(runtime, "units"),"\n"))
  return(ListVEM)
})
saveRDS(VEMr,
        file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/seed1_r0_1-6.rds")        

cliques_spca6 <- boot_FitSparsePCA(scale(counts),100,r=6)
init=initVEM(counts = counts, initviasigma=cliques_spca6[[1]], sigma_obs,r = 6)
Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
VEM_6<-VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, 
             eps=eps, alpha=1,verbatim = FALSE,
             maxIter=maxIter, plot=TRUE,vraiOm=NULL, print.hist=FALSE, filterPg=FALSE,
             filterWg = FALSE)
 


lapply(cliques_spca6, function(clique){
  length(unique(clique))==length(clique)
})
length((cliques_spca6))
#--
crit0<-criteria(list(VEM_0),counts,theta, matcovar,r=0)
critr<-do.call(rbind, lapply(seq_along(VEMr), function(r){
  criteria(VEMr[[r]],counts,theta, matcovar,r=r)
}))
crit=rbind(crit0, critr)

VEMr[[6]][[1]]$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>%
  ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(parameter)), size=3)+geom_line()+
  facet_wrap(~key, scales="free")+
  labs(x="iteration",y="", title="Lower bound and components")+mytheme+
  scale_color_discrete("")

############
 
stat=crit %>% as_tibble() %>% 
  mutate(penICL = J-ICL) %>% group_by(r) %>% 
  summarise(mvBIC=median(vBIC),mICL=median(ICL),mJ=median(J), mpenICL=median(penICL), nr=n())
model=lm(stat$mJ~stat$r)
quant=summary(model)
quant$coefficients[2]/(p*n)
crit %>% gather(key, value, -r ) %>% filter(key%in%c("ICL","J")) %>% 
  ggplot(aes(as.factor(r), value, color=as.factor(key), fill=as.factor(key)))+
  geom_boxplot(alpha=0.3)+mytheme.dark("")+
  labs(x="number of missing actors",y="values",title=)


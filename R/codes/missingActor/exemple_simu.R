library(EMtree)
library(PLNmodels)
library(tidyverse)
library(mvtnorm)
library(useful)
library(MASS)
library(ROCR)
library(reshape2)
library(parallel)
library(sparsepca)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")

simu_varyr<-function(p,n,trueR,B,maxH,maxIter,seed,type,eps=1e-3,alpha=1, cores=3, plot=FALSE){
  # -----------------------------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   Applique le VEMtree initialisées avec & à maxH acteurs manquants à des données
  #   simulées avec trueR
  # INPUT
  #   p,n     :  nombre de variables et d'individus
  #   trueR : vrai nombre d'acteurs manquants
  #   B     : nombre d'échantillons bootstrap pour sPCA
  #   maxH   : nombre maximal d'acteurs manquants testé
  #   eps : précision pour l'arrêt de VEMtree
  #   maxIter : nombre maximal d'itérations de VEMtree
  #   seed : graine
  #   type : type de graph pour la structure de dépendance ("scale-free","erdos" ou "cluster")
  #   cores : nombre de coeurs pour mclapply
  # OUTPUT
  #   VEM_r0    : résultats de VEMtree sans acteurs manquants
  #   VEMr        : liste de taille maxH dont chaque élément est la liste des résultats de VEMtree
  #                 pour chaque initialisation de sPCA et pour le nombre correspondant d'acteurs 
  #                 manquants
  #   counts    : comptages originaux simulés
  #   TrueOmega : matrice omega d'origine simulée
  #   theta : matrice des coefficients de régression estimée par PLN sur les comptages originaux
  # -----------------------------------------------------------------------------------------------------------------------------
  
  set.seed(seed)
  # sim data
  missing_data<-missing_from_scratch(n,p,r=trueR,type,plot)
  counts=missing_data$Y; omega=missing_data$Omega ;
  if(r!=0){
    omega=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden),hidden)]
  } 
  # Observed parameters
  PLNfit<-PLN(counts~1, control=list(trace=0))
  MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; 
  theta=PLNfit$model_par$Theta; sigma_obs=PLNfit$model_par$Sigma
  #initialize 
  init0=initVEM(counts , initviasigma = NULL,  sigma_obs,r = 0)
  Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 
  # vem r 0
  VEM_r0<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
                  Wg_init =Wginit,plot = plot, maxIter = maxIter,print.hist = FALSE,
                  vraiOm = NULL, alpha=alpha, verbatim=FALSE, filterPg = TRUE, 
                  filterWg = TRUE )
  #  vary r
  VEMr<-lapply(1:maxH, function(r){
    cat(paste0(r," missing actors: "))
    t1<-Sys.time()
    cliques_spca <- boot_FitSparsePCA(scale(counts),B,r=r)
    ListVEM<-List.VEM(cliqueList=cliques_spca, counts, sigma_obs, MO,SO,alpha,r=r, cores=cores,
                      eps=eps,maxIter)
    t2<-Sys.time()
    runtime=difftime(t2,t1)
    cat(paste0(round(runtime,3), attr(runtime, "units"),"\n"))
    return(ListVEM)
  })
  return(list(VEM_r0=VEM_r0, VEMr=VEMr, counts=counts, TrueOmega=omega, theta=theta))
}
#----- exemple
#-- simu
p=14 ; n=200 ; r=1 ; B=20 ; maxH=2 ; seed=1 ; type="scale-free"
sim_r1<-simu_varyr(p,n,r,B,maxH,maxIter,seed,type,eps=1e-3,alpha=1, cores=3) # en 20s

#- valeurs
#sim_r1$VEM_r0 contient les objets M, S, Pg, Wg, W, omega, lowbound, features et finalIter
ggimage(sim_r1$VEM_r0$Pg)
#sim_r1$VEMr[[1]] contient tous les ajustements de VEMtree pour 1 acteur manquant
#sim_r1$VEMr[[2]] pour deux etc.
length(sim_r1$VEMr[[1]]) # nombre d'initialisations différentes essayées pour 1 acteur manquant
#sim_r1$VEMr[[r]][[x]] contient les objets M, S, Pg, Wg, W, omega, lowbound, features et finalIter

#-- criteres
crit0<- criteria(list(sim_r1$VEM_r0),sim_r1$counts,sim_r1$theta,
                 matcovar=matrix(1, nrow(sim_r1$counts),1),r=0)
critr<-do.call(rbind, lapply(seq_along(sim_r1$VEMr), function(r){
  criteria(sim_r1$VEMr[[r]],sim_r1$counts,sim_r1$theta, 
           matcovar=matrix(1, nrow(sim_r1$counts),1),r=r)
}))
critr0=rbind(crit0, critr)

#-- plot
critr0  %>%  gather(key, value, -r ) %>% 
  filter(key%in%c("ICL","J")) %>% 
  ggplot(aes(as.factor(r), value, color=as.factor(r==0), fill=as.factor(r==0)))+
  geom_boxplot(alpha=0.3)+mytheme.dark("")+
  facet_wrap(~key, scale="free")+guides(color=FALSE,fill=FALSE)+
  labs(x="number of missing actors",y="",title="seed=1, r=0")

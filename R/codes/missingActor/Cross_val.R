# tirer un arbre
# heuristique
library(ape)
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
source("/Users/raphaellemomal/these/Pgm_SR/TreeSampling.R")
# rtree<-function(P){
#   # créer une matrice de poids en tirant uniformément entre 0 et chaque proba
#   p=ncol(P)
#   weights=apply(P, 2, function(x) Vectorize(runif)(1,0,x))
#   #trouver max spann tree (package ape)
#   max_tree=matrix(mst(-weights), p,p)
#   return(max_tree)
# }

split_fit<-function(Y,v=0.8,r=1){
  p=ncol(Y) ;n=nrow(Y); O=1:p ; H=(p+1):(p+r)
  # data
  ntrain=round(n*v) ; ntest=n-ntrain ; samp = sample(1:n, n*v)
  Ytrain =counts[samp,] ; Ytest=counts[-samp,]
  # Ytrain=Ytest=Y
  # ntrain=ntest=n
  #-- normalized PLN outputs
  PLNfit<-PLN(Ytrain~1, control=list(trace=0))
  MO<-PLNfit$var_par$M  ;SO<-PLNfit$var_par$S  ;sigma_obs=PLNfit$model_par$Sigma
  D=diag(sigma_obs)
  matsig=(matrix(rep(1/sqrt(D),ntrain),ntrain,p, byrow = TRUE))
  MO=MO*matsig
  SO=SO*matsig^2
  
  if(r!=0){
    # init
    cliques_spca<-FitSparsePCA(Ytrain, r=2)$cliques
    complement=lapply(cliques_spca, function(clique){setdiff(1:p,clique)})
    clique=list()
    clique$cliqueList=lapply(c(cliques_spca,complement), function(cl) list(cl))
    #VEM
    ListVEM<-List.VEM(cliquesObj =clique, Ytrain, cov2cor(sigma_obs), MO,SO,r=r,alpha=0.1,
                      eps=1e-3,maxIter=100, cores=3, trackJ=FALSE)
    vecJ=do.call(rbind, lapply(ListVEM, function(vem){
      if(length(vem)==12){ J=tail(vem$lowbound$J, 1)
      }else{ J=NaN} }))
    VEM=ListVEM[[which.max(vecJ)]]
  }else{
    init0=initVEM(Ytrain , initviasigma = NULL,  cov2cor(sigma_obs),r = 0)
    Wginit= init0$Wginit; Winit= init0$Winit; upsinit=init0$upsinit 
    VEM<-tryCatch({VEMtree(Ytrain, MO, SO, MH=NULL,upsinit,W_init =Winit,eps=1e-3,
                           Wg_init =Wginit,plot = FALSE, maxIter = 100,print.hist = FALSE,
                           alpha=0.1, verbatim=FALSE, trackJ=FALSE )},
                  error=function(e){e}, finally={})
  }
  return(list(Ytest=Ytest,Pg=VEM$Pg,beta=VEM$Wg,Omega_hat=VEM$Upsilon,D=D))
}

pY_condT<-function(Ytest,beta, prob, Omega_hat,D, r=1, plot=FALSE){
  p=ncol(Ytest) ;n=nrow(Ytest); O=1:p ; H=(p+1):(p+r)
  prob[prob>1]=1
  beta <- beta / SumTree(beta)^(1/(ncol(beta)-1))
 
  #--- calcul critere : tirer T, calculer omega et sigma puis Uo et p(Y)
  # obj.tree = rSpanTreeV1(beta=beta,prob=prob)
  # tree=obj.tree$tree
  tree=rSpanTreeV2(beta)
  OmegaT=(tree+diag(ncol(prob)))*Omega_hat
  #  OmegaT = (Pg+diag(ncol(Pg)))*VEM$Upsilon
  if(r!=0){
    SigmaTm= solve(OmegaT[O,O]- OmegaT[O,H]%*%solve(OmegaT[H,H])%*%OmegaT[H,O])
  }else{
    SigmaTm=solve(OmegaT[O,O])
  }
  UO<-rmvnorm(n=n,sigma=SigmaTm)
  
  PLNfit_test<-PLN(Ytest~1, control=list(trace=0)) # calcul des theta test
  X=matrix(1, n, 1) ; theta= PLNfit_test$model_par$Theta
  lambda = (X%*%t(theta)+UO*(rep(1, nrow(UO))%o%D))
  logpY= sum(dpois(Ytest, exp(lambda), log=TRUE))
  if(plot) plot((lambda), log(Ytest+1), pch=20);abline(0,1)
  
  return(list(logpY=logpY, tree=tree))
}

#-- data simulation
set.seed(1) 
n=200 ;p=14;r=1;type="scale-free";plot=TRUE
O=1:p
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y
tic()
VEMfit=split_fit(counts, v=0.8,r=0)
toc()#6sec
saveRDS(VEMfit, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/VEMfit.R")
tic()
res=lapply(1:50, function(x){#1min
  pY_condT(Ytest = VEMfit$Ytest,beta = VEMfit$beta,prob = VEMfit$Pg,
           Omega_hat = VEMfit$Omega_hat,D=VEMfit$D, r = 0,plot = TRUE)
})
toc()

#toujours le même arbre
colSums(do.call(rbind,lapply(res, function(x){
  F_Sym2Vec(x$tree)
})))

# des proba numériquement nulles
vec_logpY = (do.call(rbind, lapply(res, function(x) x$logpY)))
summary(vec_logpY)

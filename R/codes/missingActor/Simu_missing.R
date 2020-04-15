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

Simu_missing<-function(p,B,N,n,cores,r, maxIter, eps){
  lapply(1:N, function(seed){
    cat(paste0("\n seed ",seed, " : "))
    T1<-Sys.time()
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
    
    #1 missing actors
    t1<-Sys.time()
    cliques_spca <- boot_FitSparsePCA(scale(MO),B,r=1)
    t2<-Sys.time()
    time_boots=difftime(t2, t1)
    # best VEM with 1 missing actor
  
    # alpha<-tryCatch(expr={computeAlpha(omegainit,default =1/n, MO, SO, plot=plot)},
    # error=function(e){message("sythetic alpha")
    #   return(1/n)}) #0.46
    
    ListVEM<-List.VEM(cliquesObj =cliques_spca, counts, sigma_obs, MO,SO,r=1, cores=cores,
                      eps=eps,maxIter=maxIter)
    crit1<-criteria(ListVEM,counts,theta, matcovar,r=1)
    
    ############
    nbconv=length(ListVEM)
    best=which.max(crit1$J)
    VEM_1<-ListVEM[[best]]
    
    #---- end
    
    T2<-Sys.time()
    runtime=difftime(T2, T1)
    cat(paste0(round(runtime,3), attr(runtime, "units"),"\n"))
    return(list(omega=sorted_omega,ZH=ZH,VEM_1=VEM_1,time_boots=time_boots, nbconv=nbconv ))
  })
}

######### run
 
Sim15<-Simu_missing(p = 14, n = 200, B = 100,N = 200, cores=1,eps=1e-3,r=1,maxIter=100)
saveRDS(Sim15, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15_r1_200SF.rds")
# Sim30<-Simu_missing(p = 29, n = 200, B = 40,N = 200)
# saveRDS(Sim15, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim30.rds")
Sim15_r0<-Simu_missing(p = 14, n = 200, B = 40,N = 30, cores=3,r=0,maxIter=100)
saveRDS(Sim15_r0, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15_r0.rds")

#pb seed 33 avec 30 noeuds
seed=33

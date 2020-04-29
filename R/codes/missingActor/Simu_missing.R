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
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-exactDet.R")

##### function

Simu_missing<-function(p,B,N,n,cores,r, maxIter, eps){
  lapply(363:N, function(seed){
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
    cliques_spca <- boot_FitSparsePCA(scale(MO),B,r=1, cores=3)
    t2<-Sys.time()
    time_boots=difftime(t2, t1)
    # best VEM with 1 missing actor
    ListVEM<-List.VEM(cliquesObj =cliques_spca, counts, sigma_obs, MO,SO,r=1,alpha=0.1,
                      eps=eps,maxIter=maxIter, nobeta=FALSE, cores=cores)
    
    goodPrec=!do.call(rbind,lapply(ListVEM, function(x) x$max.prec))
    J=do.call(rbind,lapply(ListVEM, function(vem){tail(vem$lowbound$J,1)}))
    
    # if(sum(goodPrec)!=0){ # tri spécifique si filtreWg est FALSE
    #   if(sum(J<min(J[!goodPrec]))!=0){
    #     maxJ_good=which(J==max(J[J<min(J[!goodPrec])])) 
    #   }else{
    #     maxJ_good = which(J==max(J[goodPrec]))
    #   }
    # }else{
    #   maxJ_good = which.max(J)
    # } 
    maxJ_good=which.max(J)
    VEM_1=ListVEM[[maxJ_good]]
    
    ############
    nbinit=length(ListVEM)
    
    #---- end
    
    T2<-Sys.time()
    runtime=difftime(T2, T1)
    cat(paste0("\nseed ", seed," in ",round(runtime,3), attr(runtime, "units"),"\n"))
    Sim=list(omega=sorted_omega,ZH=ZH,VEM_1=VEM_1,time_boots=time_boots, nbinit=nbinit )
    saveRDS(Sim, file=paste0("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/15nodes_1r_400_filtreWg/SF_seed",
                             seed,".rds"))
    
    return(Sim)
  })
}

######### run
t1<-Sys.time()
Sim15<-Simu_missing(p = 14, n = 200, B = 100,N =400, cores=3,eps=1e-3,r=1,maxIter=200)
t2<-Sys.time()
difftime(t2,t1)
saveRDS(Sim15, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15_r1_200SF.rds")
# Sim30<-Simu_missing(p = 29, n = 200, B = 40,N = 200)
# saveRDS(Sim15, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim30.rds")
Sim15_r0<-Simu_missing(p = 14, n = 200, B = 40,N = 30, cores=3,r=0,maxIter=100)
saveRDS(Sim15_r0, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15_r0.rds")

#pb seed 33 avec 30 noeuds
seed=33

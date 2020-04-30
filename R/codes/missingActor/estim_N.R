library(EMtree)
library(PLNmodels)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(mvtnorm)
library(useful)
library(MASS)
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(parallel)
library(sparsepca)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")


estim_N<-function(p,B=100,nSF=200,n,cores,r, maxIter, eps){
  cliques<-lapply(1:nSF, function(seed){
    cat(paste0("\n seed ",seed, " : "))
    set.seed(seed) ; type="scale-free" ; O=1:p ; plot=FALSE 
    # Data
    missing_data<-missing_from_scratch(n,p,r,type,plot)
    counts=missing_data$Y; ZH=missing_data$ZH ; sigmaO= missing_data$Sigma; 
    omega=missing_data$Omega; trueClique=missing_data$TC[[1]]; hidden=missing_data$H
    # Observed parameters
    PLNfit<-PLN(counts~1, control=list(trace=0))
    MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; theta=PLNfit$model_par$Theta ; 
    matcovar=matrix(1, n,1) ; sigma_obs=PLNfit$model_par$Sigma
    
    #1 missing actors
    t1<-Sys.time()
    cliques_spca <- boot_FitSparsePCA(scale(MO),B,r=1, cores=3, unique=FALSE)
    t2<-Sys.time()
    time_boots=difftime(t2, t1)
    cat(paste0("\nseed ", seed," in ",round(time_boots,3), attr(time_boots, "units"),"\n"))
    return(cliques_spca)
  })
  return(cliques)
}


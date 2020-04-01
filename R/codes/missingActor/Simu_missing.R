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

Simu_missing<-function(p,B,N,n){
  sapply(1:N, function(seed){
    cat(paste0("\n seed ",seed, " : "))
    T1<-Sys.time()
    set.seed(seed)
    type="scale-free" ; O=1:p ; plot=FALSE ; r=1
    # Data
    missing_data<-missing_from_scratch(n,p,r,type,plot)
    counts=missing_data$Y; sigmaO= missing_data$Sigma; omega=missing_data$Omega;
    trueClique=missing_data$TC[[1]]; hidden=missing_data$H
    # Observed parameters
    PLNfit<-PLN(counts~1, control=list(trace=0))
    MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; theta=PLNfit$model_par$Theta ; 
    matcovar=matrix(1, n,1) ; sigma_obs=PLNfit$model_par$Sigma
    # initalisation with 1 missing actors
    init=initVEM(counts = counts, initviasigma = trueClique,  sigma_obs,r = r)
    Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
    ome_init=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), hidden)]
    
    t1<-Sys.time()
    cliques_spca <- boot_FitSparsePCA(scale(counts),B,r=r)
    t2<-Sys.time()
    time_boots=difftime(t2, t1)
    # best VEM with 1 missing actor
    ListVEM<-List.VEM(cliqueList=cliques_spca, counts, sigma_obs, MO,SO,r=r, cores=3)
    time_1<-mean(do.call(rbind, lapply(ListVEM, function(vem){vem$time})))
    vBICs<-as.numeric(vec.vBIC(ListVEM,counts,theta, matcovar,r))
    best=which.max(vBICs)
    vBIC_1<-vBICs[best]
    VEM_1<-ListVEM[[best]]
    VEM_1$time=time_1
    #VEM with no missing actor
    init0=initVEM(counts = counts, initviasigma = NULL,  sigma_obs,r = 0)
    Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 
    alpha<-tryCatch(expr={computeAlpha(omegainit,default =0.3, MO, SO, plot=plot)},
                    error=function(e){message("sythetic alpha")
                      return(0.3)})
    VEM_0<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,
                   Wg_init =Wginit,plot = plot, maxIter = 30,print.hist = FALSE,
                   vraiOm = NULL, alpha=alpha, verbatim=FALSE )
    VEM_0$alpha=alpha
    J0<-True_lowBound(counts,VEM_0$M,VEM_0$S, theta, matcovar,VEM_0$W, VEM_0$Wg,
                      VEM_0$Pg, VEM_0$omega )
    vBIC_0<-VBIC(J0, ncol(counts),r=0, d=1, n=nrow(counts))
    #---- end
    T2<-Sys.time()
    runtime=difftime(T2, T1)
    cat(paste0(round(runtime,3), attr(runtime, "units"),"\n"))
    return(list(omega=ome_init,VEM_0=VEM_0,VEM_1=VEM_1,vBIC_0=vBIC_0,vBIC_1=vBIC_1,time_boots=time_boots ))
  })
}

######### run
Sim15<-Simu_missing(p = 14, n = 200, B = 40,N = 200)
saveRDS(Sim15, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15.rds")




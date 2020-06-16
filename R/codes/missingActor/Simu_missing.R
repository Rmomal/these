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
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-V2.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/modif_pkg.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/VEM_tools.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-exactDet.R")

##### function

Simu_missing<-function(p,B,N,n,cores,r, maxIter, eps){
  lapply(1:N, function(seed){
    cat(paste0("\n seed ",seed, " : "))
    T1<-Sys.time()
    set.seed(seed)
    type="scale-free" ; O=1:p ;H=(p+1):(p+r); plot=FALSE 
    # Data
    missing_data<-missing_from_scratch(n,p,r,type,plot)
    counts=missing_data$Y; UH=missing_data$UH  
    G=missing_data$G; trueClique=missing_data$TC[[1]]; hidden=missing_data$H
    # Observed parameters
    PLNfit<-PLN(counts~1, control=list(trace=0))
    MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; sigma_obs=PLNfit$model_par$Sigma
  
    #-- normalize the PLN outputs
    D=diag(sigma_obs)
    matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
    MO=MO*matsig
    SO=SO*matsig^2
    #1 missing actors
    t1<-Sys.time()
  #  cliques_spca <- boot_FitSparsePCA(scale(MO),B,r=1, cores=3)
   
    # best VEM with 1 missing actor
    # cliques_spca<-FitSparsePCA(counts, r=2)$cliques
    # complement=lapply(cliques_spca, function(clique){setdiff(1:p,clique)})
    # clique=list()
    # t2<-Sys.time()
    # time_spca=difftime(t2, t1)
    # clique$cliqueList=lapply(c(cliques_spca,complement), function(cl) list(cl))
    #--- init oracle
    clique=list()
    clique$cliqueList=list(list(trueClique))

    ListVEM<-List.VEM(cliquesObj =clique, counts, cov2cor(sigma_obs), MO,SO,r=1,alpha=0.1,
                      eps=eps,maxIter=maxIter, cores=3, trackJ = FALSE)
    
    # Jcor<-do.call(rbind, lapply(ListVEM, function(vem){
    #   if(length(vem)==15){
    #     res=getJcor(vem,14)[1]
    #   }else{
    #     res=NA
    #   }
    #   return(res)
    # }))
    # if(sum(!is.na(Jcor))!=0){
    #   sumP=do.call(rbind, lapply(ListVEM, function(vem){(sum(vem$Pg)-28)<1e-10}))
    #   data=data.frame(Jcor,sumP, num=1:length(ListVEM))
    #   maxJ_good=unlist(data %>% filter(sumP) %>% filter(Jcor==max(Jcor, na.rm=TRUE)) %>%
    #                 dplyr::select(num))
    # }else{
    #   maxJ_good = which.max(do.call(rbind, lapply(ListVEM, function(vem){
    #    tail(vem$lowbound$J,1)
    #   })))
    # }
    # 
    # VEM_1=ListVEM[[maxJ_good]]
    # 
    ############
   # nbinit=length(ListVEM)
    
    #---- end
    
    T2<-Sys.time()
    runtime=difftime(T2, T1)
    cat(paste0("\nseed ", seed," in ",round(runtime,3), attr(runtime, "units"),"\n"))
    Sim=list(G=G,UH=UH,
             ListVEM=ListVEM)#,#VEM_1=VEM_1,
             #time_spca=time_spca)#, nbinit=nbinit )
    saveRDS(Sim, file=paste0("/Users/raphaellemomal/simulations/15nodes_V4_oracle/SF_seed",
                             seed,".rds"))
    
    return(Sim)
  })
}

######### run
tic()
Sim15<-Simu_missing(p = 14, n = 200, B = 100,N = 403,eps = 1e-3, cores=3,r=1,maxIter=200)
toc()
saveRDS(Sim15, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15_r1_200SF.rds")
# Sim30<-Simu_missing(p = 29, n = 200, B = 40,N = 200)
# saveRDS(Sim15, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim30.rds")
Sim15_r0<-Simu_missing(p = 14, n = 200, B = 40,N = 30, cores=3,r=0,maxIter=100)
saveRDS(Sim15_r0, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15_r0.rds")

#pb seed 33 avec 30 noeuds
seed=33

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

Simu_missing<-function(p,B,N,n,cores,r, maxIter){
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
    #VEM with no missing actor
    init0=initVEM(counts = counts, initviasigma = NULL,  sigma_obs,r = 0)
    Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 
    
   #alpha=1/n
    # VEM_0$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>%
    #   ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(parameter)), size=3)+geom_line()+
    #   facet_wrap(~key, scales="free")+
    #   labs(x="iteration",y="", title="Lower bound and components")+mytheme+
    #   scale_color_discrete("") 
    # VEM_0$features  %>%  rowid_to_column() %>%
    #   pivot_longer(-rowid,names_to="key",values_to = "values") %>%
    #   ggplot(aes(rowid,values, color=key))+
    #   geom_point()+geom_line() + facet_wrap(~key, scales="free")+
    #   labs(x="",y="", title="Parameters")+ mytheme.dark+guides(color=FALSE)
    # ome=sorted_omega ; diag(ome)=0
    # g1<-ggimage(VEM_0$Pg) ; g2<-ggimage(ome) ; grid.arrange(g1, g2, ncol=2)
    # 
    #r missing actors
    t1<-Sys.time()
    cliques_spca1 <- boot_FitSparsePCA(scale(counts),B,r=1)
    cliques_spca2 <- boot_FitSparsePCA(scale(counts),B,r=2)
    t2<-Sys.time()
    time_boots=difftime(t2, t1)
    # best VEM with 1 missing actor
    eps=1e-3
    # alpha=1
    alpha=1
    VEM_01<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
                   Wg_init =Wginit,plot = plot, maxIter = maxIter,print.hist = FALSE,
                   vraiOm = NULL, alpha=alpha, verbatim=FALSE, filterPg = TRUE, filterWg = TRUE )
    ListVEM11<-List.VEM(cliqueList=cliques_spca1, counts, sigma_obs, MO,SO,alpha,r=1, cores=3,
                       eps=eps,maxIter)
    ListVEM21<-List.VEM(cliqueList=cliques_spca2, counts, sigma_obs, MO,SO,alpha,r=2, cores=3,
                      maxIter, eps=eps)
    crit0<-criteria(list(VEM_01),counts,theta, matcovar,r=0)
    crit1<-criteria(ListVEM11,counts,theta, matcovar,r=1)
    crit2<-criteria(ListVEM21,counts,theta, matcovar,r=2)
    crit_alpha1=rbind(crit0, crit1,crit2)
    crit_alpha1$alpha="1"
    # alpha=1/n
    alpha=1/n
    VEM_02<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
                    Wg_init =Wginit,plot = plot, maxIter = maxIter,print.hist = FALSE,
                    vraiOm = NULL, alpha=alpha, verbatim=FALSE, filterPg = TRUE, filterWg = TRUE )
    ListVEM12<-List.VEM(cliqueList=cliques_spca1, counts, sigma_obs, MO,SO,alpha,r=1, cores=3,
                        eps=eps,maxIter)
    ListVEM22<-List.VEM(cliqueList=cliques_spca2, counts, sigma_obs, MO,SO,alpha,r=2, cores=3,
                        maxIter, eps=eps)
    crit0<-criteria(list(VEM_02),counts,theta, matcovar,r=0)
    crit1<-criteria(ListVEM12,counts,theta, matcovar,r=1)
    crit2<-criteria(ListVEM22,counts,theta, matcovar,r=2)
    crit_alpha2=rbind(crit0, crit1,crit2)
    crit_alpha2$alpha="1/n"
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
    crit0<-criteria(list(VEM_03),counts,theta, matcovar,r=0)
    crit1<-criteria(ListVEM13,counts,theta, matcovar,r=1)
    crit2<-criteria(ListVEM23,counts,theta, matcovar,r=2)
    crit_alpha3=rbind(crit0, crit1,crit2)
    crit_alpha3$alpha="tuned"
    ############
    crit=rbind(crit_alpha1, crit_alpha2, crit_alpha3)
    crit %>% gather(key, value, -r, -alpha) %>%
      ggplot(aes(as.factor(r), value, color=as.factor(alpha), fill=as.factor(alpha)))+geom_boxplot(alpha=0.3)+
      facet_wrap(~as.factor(key))+mytheme.dark+
      labs(x="number of missing actors",y="values",title="Crit from boot.sPCA (B=100)")
    
    time_1<-mean(do.call(rbind, lapply(ListVEM, function(vem){vem$time})))
    best=which.max(crit1$vBIC)
    vBIC_1<-crit1$vBIC[best]
    VEM_1<-ListVEM[[best]]
    VEM_1$time=time_1
      #---- end
    T2<-Sys.time()
    runtime=difftime(T2, T1)
    cat(paste0(round(runtime,3), attr(runtime, "units"),"\n"))
    return(list(omega=sorted_omega,ZH=ZH,alpha, VEM_0=VEM_0,VEM_1=VEM_1,
                vBIC_0=vBIC_0,vBIC_1=vBIC_1,time_boots=time_boots ))
  })
}

######### run
Sim15<-Simu_missing(p = 14, n = 200, B = 40,N = 200, cores=3,r=1,maxIter=100)
saveRDS(Sim15, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15conv.rds")
# Sim30<-Simu_missing(p = 29, n = 200, B = 40,N = 200)
# saveRDS(Sim15, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim30.rds")
Sim15_r0<-Simu_missing(p = 14, n = 200, B = 40,N = 30, cores=3,r=0,maxIter=100)
saveRDS(Sim15_r0, file="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15_r0.rds")

#pb seed 33 avec 30 noeuds
seed=33

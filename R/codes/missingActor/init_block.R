library(mvtnorm)
library(PLNmodels)
library(tictoc)
library(EMtree)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(ROCR)
library(blockmodels)
# bad inits for those with NA results and those that did not converge
simus=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/simus_V2_400_4init.rds")
badseeds=which(do.call(rbind, lapply(simus, function(seed){
  length(seed)
}))!=4)
finalIter=do.call(rbind, lapply(simus, function(graine){
  finalIter=graine$VEM_1$finalIter
}))
noconv=which(finalIter==100)
badinits=c(badseeds,noconv ) # 11 : 35  73  89 120 162 206 342 367 395  45 273
nH=do.call(rbind,lapply(badinits, function(seed){
  set.seed(seed); p=14;n=200;r=1
  type="scale-free" ; O=1:p ;H=(p+1):(p+r); plot=FALSE 
  # Data
  missing_data<-missing_from_scratch(n,p,r,type,plot)
  counts=missing_data$Y; UH=missing_data$UH  
  G=missing_data$G
  return(sum(G[,15]))
}))
badinits[7]
#------ function

Simu_block<-function(seed, k=3){
  cat(paste0("\n seed ",seed, " : "))
  #----- generate data
  T1<-Sys.time()
  set.seed(seed); p=14;n=200;r=1
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
  
  #------- fit VEMTree
  T3<-Sys.time()
  cliques<-init_blockmodels(k, counts, sigma_obs, MO, SO)
  T4<-Sys.time()
  time_block=difftime(T4, T3)
  ListVEM<-List.VEM(cliquesObj =cliques, counts, cov2cor(sigma_obs), MO,SO,r=1,alpha=0.1,
                    eps=1e-3,maxIter=100, cores=3, trackJ = FALSE)
  vecJ=do.call(rbind, lapply(ListVEM, function(vem){
    if(length(vem)==12){
      J=tail(vem$lowbound$J, 1)
    }else{
      J=NaN}
  }))
  cat(paste0("\n", vecJ))
  T2<-Sys.time()
  runtime=difftime(T2, T1)
  cat(paste0("\nseed ", seed," in ",round(runtime,3), attr(runtime, "units"),"\n"))
  Sim=list(G=G,UH=UH,
           ListVEM=ListVEM,#VEM_1=VEM_1,
           time_block=time_block)#, nbinit=nbinit )
  saveRDS(Sim, file=paste0("/Users/raphaellemomal/simulations/15nodes_V2_blockinit/SF_seed",
                           seed,".rds"))
  
  return(Sim)
}

#---- run
lapply(badinits[-3], function(seed){
  Simu_block(seed)
})



VEM1=ListVEM[[which.max(vecJ)]]
plotVEM(ListVEM[[1]]$Pg, G, 1, 0.5)
VEM1$lowbound %>% rowid_to_column() %>%  
  ggplot(aes(rowid,J ))+geom_line()+geom_point(aes(color=as.factor(parameter)),size=2, alpha=0.8)+
  labs(x="iteration",y="", title="Lower bound")+theme_light()+ scale_color_manual("",values="#2976d6")+
  guides(color=FALSE)
VEM1$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>% 
  ggplot(aes(rowid,value, group=key))+geom_line()+geom_point(aes(color=as.factor(parameter)), size=2, alpha=0.8)+
  facet_wrap(~key, scales="free")+  guides(color=FALSE)+ labs(x="sub-iteration",y="", title="Lower bound and components")+mytheme.dark("")

#--------- test single
init=initVEM(counts = counts,initviasigma=list(trueClique), cov2cor(sigma_obs),MO,r = 1) 
Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
resVEM<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.1,
                 maxIter=30, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=TRUE)
plotVEM(resVEM$Pg, G, 1, 0.5)




#---------- test single r=0
set.seed(342); p=14;n=200;r=1
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
init=initVEM(counts = counts,initviasigma=NULL, cov2cor(sigma_obs),MO,r = 0) 
Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
resVEM0<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.1,
                  maxIter=100, plot=FALSE,print.hist=FALSE, verbatim = TRUE,trackJ=FALSE)
resVEM0$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>% 
  ggplot(aes(rowid,value, group=key))+geom_line()+geom_point(aes(color=as.factor(parameter)), size=2, alpha=0.8)+
  facet_wrap(~key, scales="free")+  guides(color=FALSE)+ labs(x="sub-iteration",y="", title="Lower bound and components")+mytheme.dark("")
resVEM0$features$diffUpsi

#------ etude 342 nH=4
SF_seed162 <- readRDS("/Users/raphaellemomal/simulations/15nodes_V2_blockinit/SF_seed162.rds")
G=SF_seed162$G
vecJ=do.call(rbind, lapply(SF_seed162$ListVEM, function(vem){
  if(length(vem)==12){
    J=tail(vem$lowbound$J, 1)
  }else{
    J=NaN}
}))
VEM1=SF_seed162$ListVEM[[which.max(vecJ)]]
VEM1$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>% 
  ggplot(aes(rowid,value, group=key))+geom_line()+geom_point(aes(color=as.factor(parameter)), size=2, alpha=0.8)+
  facet_wrap(~key, scales="free")+  guides(color=FALSE)+ labs(x="sub-iteration",y="", title="Lower bound and components")+mytheme.dark("")
plotVEM(VEM1$Pg, G, 1, 0.5)

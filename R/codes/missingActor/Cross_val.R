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
library(sna)
library(poilog)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-V2.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/modif_pkg.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/VEM_tools.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-exactDet.R")


# rtree<-function(P){
#   # créer une matrice de poids en tirant uniformément entre 0 et chaque proba
#   p=ncol(P)
#   weights=apply(P, 2, function(x) Vectorize(runif)(1,0,x))
#   #trouver max spann tree (package ape)
#   max_tree=matrix(mst(-weights), p,p)
#   return(max_tree)
# }

split_fit<-function(Y,v=0.8,r=1,alpha=0.1, bloc=1,eps=1e-3){#shuffle data prior to this
  p=ncol(Y) ;n=nrow(Y); O=1:p ; H=(p+1):(p+r)
  # data
  T1<-Sys.time()
  ntrain=round(n*v) ; ntest=n-ntrain ; nblocs=round(n/ntest)
  vec_blocs=rep(0,n)
  vec_blocs[1:((nblocs-1)*ntest)]=rep(1:(nblocs-1), each=ntest)
  vec_blocs[vec_blocs==0]=nblocs
  samp=which(vec_blocs==bloc)
  Ytrain =Y[-samp,] ; Ytest=Y[samp,]
  ntrain=nrow(Ytrain)

  #-- normalized PLN outputs
  PLNfit<-PLN(Ytrain~1, control=list(trace=0))
  MO<-PLNfit$var_par$M  ;SO<-PLNfit$var_par$S  ;sigma_obs=PLNfit$model_par$Sigma
  D=diag(sigma_obs)
  matsig=(matrix(rep(1/sqrt(D),ntrain),ntrain,p, byrow = TRUE))
  MO=MO*matsig
  SO=SO*matsig^2
  
  if(r!=0){
    if(r==1){minV=2
    }else{minV=0}
      clique1=boot_FitSparsePCA(Y, B=200,r=r,minV=minV)
      clique2=boot_FitSparsePCA(Ytrain, B=200,r=r,minV=minV)
      clique=list()
      clique$cliqueList=unique(c(clique1$cliqueList,clique2$cliqueList))
 
      # if(length(clique$cliqueList)>50){
      #   clique_spca<-FitSparsePCA(counts, r=2)$cliques
      #   complement=lapply(clique_spca, function(clique){setdiff(1:p,clique)})
      #   clique=list()
      #   clique$cliqueList=lapply(c(clique_spca,complement), function(cl) list(cl))
      # }
    #VEM
    cat(paste0("\n   Fitting ",length(clique$cliqueList)," VEM..."))
    ListVEM<-List.VEM(cliquesObj =clique, Ytrain, cov2cor(sigma_obs), MO,SO,r=r,alpha=alpha,
                      eps=eps,maxIter=100, cores=3, trackJ=FALSE)

    if(r>1){ # select best J among models which gives r Mh with reasonabe sd 
      vec_sdMH_J=do.call(rbind, lapply(ListVEM, function(vem){
        if(length(vem)>5){
          J=tail(vem$lowbound$J, 1)
          sdM=apply(vem$M,2, function(x) sd(x))
          sdMH=tail(sdM,r)
          res=c(sdMH,J=J)
        }else{
          J=NaN ;sdMH =rep(NaN,r)
          res=c(sdMH,J=J)}
        return(res)
      }))
      data_sd=data.frame(vec_sdMH_J) %>% dplyr::select(-J)
      AllActors=apply(data_sd, 1, function(x){
        sum(log(x)>-20)==r
      })
      if(sum(AllActors, na.rm=TRUE)!=0){
        index=unlist(data.frame(vec_sdMH_J,AllActors=AllActors,num=1:nrow(data_sd)) %>% as_tibble() %>% 
                     filter(AllActors) %>% filter(J==max(J,na.rm=TRUE)) %>% dplyr::select(num))
      }else{ index=which.max(vec_sdMH_J$J)}
      
    }else{
      vec_J=do.call(rbind, lapply(ListVEM, function(vem){
        if(length(vem)>5){
          J=tail(vem$lowbound$J, 1)
        }else{
          J=NaN }
      }))
      index=which.max(vec_J)
    }
    
    VEM=ListVEM[[index]]
  }else{
    cat(paste0("\n   Fitting 1 VEM with r=0..."))
    init0=initVEM(Ytrain , initviasigma = NULL,  cov2cor(sigma_obs),r = 0)
    Wginit= init0$Wginit; Winit= init0$Winit; upsinit=init0$upsinit 
    VEM<-tryCatch({VEMtree(Ytrain, MO, SO, MH=NULL,upsinit,W_init =Winit,eps=eps,
                           Wg_init =Wginit,plot = FALSE, maxIter = 100,print.hist = FALSE,
                           alpha=alpha, verbatim=FALSE, trackJ=FALSE )},
                  error=function(e){e}, finally={})
  }
  T2<-Sys.time()
  runtime=difftime(T2, T1)
  cat(paste0(round(runtime,3)," ", attr(runtime, "units")))
  return(list(Ytest=Ytest,Pg=VEM$Pg,beta=VEM$Wg,M=VEM$M,S=VEM$S,
              Omega_hat=VEM$Upsilon,D=D, theta=PLNfit$model_par$Theta ))
}

pY_condT<-function(Ytest,beta, prob,M,S, Omega_hat,D,theta, r=1, plot=FALSE){
  p=ncol(Ytest) ;n=nrow(Ytest); O=1:p
  prob[prob>1]=1 ; q=ncol(beta)
  if(is.finite(SumTree(beta))){
    beta <- beta / SumTree(beta)^(1/(q-1))
  }else{
    beta=beta/100
    beta <- beta / SumTree(beta)^(1/(q-1))
  }
#browser()

  #--- calcul critere : tirer T, calculer omega et sigma et p(Y)
  obj.tree = rSpanTreeV1(beta=beta,prob=prob)
  tree=obj.tree$tree
  #tree=rSpanTreeV1(beta)
  OmegaT=tree*Omega_hat
  ntrain = nrow(M)
  Rho = cov2cor((t(M)%*%M+diag(colSums(S)))/ntrain)
  quantity<-tree*Rho^2/(1-Rho^2)
  diag(quantity)=NA
  vectorSuml=colSums(quantity, na.rm=TRUE)
  diag(OmegaT) = 1+vectorSuml
  attr(OmegaT,"class")="matrix"
  
  if(r!=0){
    H=(p+1):(p+r)
    SigmaTm=cov2cor(solve(nearPD(OmegaT[O,O]- OmegaT[O,H]%*%solve(OmegaT[H,H])%*%OmegaT[H,O], 
                                 eig.tol = 1e-14, posd.tol = 1e-14)$mat))
  }else{
    SigmaTm=cov2cor(solve(nearPD(OmegaT, eig.tol = 1e-14, posd.tol = 1e-14)$mat))
  }
  #pour chaque site test i et chaque paire d’especes (j,k),
  #on calcule la densit́e Poissonlog-normale
  mat_dens=matrix(0,p,p)
  sapply(1:(p-1), function(i){
    sapply((i+1):p, function(j){
      sig=sqrt(D[i]) ; sig2=sqrt(D[j]) ; rho=SigmaTm[i,j]
      if(rho<0) rho = max(SigmaTm[i,j],-0.9999)
      if(rho>0) rho = min(SigmaTm[i,j],0.9999)
      mu1 = theta[i] ; mu2 = theta[j]
      mat_dens[i,j]<<-sum(pmax(log(dbipoilog(n1=Ytest[,i],n2=Ytest[,j],mu1=mu1,mu2=mu2,
                                             sig=sig,sig2=sig2,rho=rho)),-709))
      mat_dens[j,i]<<-mat_dens[i,j]
    })
  })
  VCp<-sum(mat_dens)/2
  return(list(VCp=VCp, tree=tree))
}

#cross_val0 pour r=0 et trueR=1
cross_val<-function(counts,r, nblocs=10,B=50,alpha=0.1,eps=1e-3, cores=3,dir){
  T0<-Sys.time()
  res=lapply(1:nblocs, function(i){
    cat(paste0("Bloc ",i,":"))
    file=paste0("/Users/raphaellemomal/simulations/",dir,"/VEMfit_r",r,"_bloc",i,".rds")
    if(!file.exists(file)){
      VEMfit=split_fit(counts, v=1-(1/nblocs),r=r,alpha=alpha, bloc=i,eps=eps)
      saveRDS(VEMfit, file=file)
    }  
    cat(paste0("\n   Estimating pairwise composite likelihood..."))
    T1<-Sys.time()
    res=mclapply(1:B, function(x){#1min
    
      VEMfit=readRDS(paste0("/Users/raphaellemomal/simulations/",dir,"/VEMfit_r",r,"_bloc",i,".rds"))
      pY_condT(Ytest = VEMfit$Ytest,beta = VEMfit$beta,prob = VEMfit$Pg,M = VEMfit$M,S = VEMfit$S,
              Omega_hat = VEMfit$Omega_hat,D=VEMfit$D,theta = VEMfit$theta,
               r = r,plot = FALSE)
    }, mc.cores=cores)
    edge_freq=colSums(do.call(rbind,lapply(res, function(x){
      F_Sym2Vec(x$tree)
    })))
    vec_logpY = do.call(rbind, lapply(res, function(x) x$VCp))
    T2<-Sys.time()
    runtime=difftime(T2, T1)
    cat(paste0(round(runtime,3), attr(runtime, "units"),"\n"))
    return(list(vec_logpY=vec_logpY,edge_freq=edge_freq))
  })
  Tn<-Sys.time()
  runtime=difftime(Tn, T0)
  cat(paste0(" in ",round(runtime,3), attr(runtime, "units"),"\n"))
  return(res)
}
#-- data simulation
set.seed(1) ;n=200 ;p=14;type="scale-free";plot=TRUE;O=1:p
r=0
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y
reorder=sample(1:n, n,replace=FALSE)
counts0=counts[reorder,]
cross_val00<-cross_val(counts0,r=0,dir="CV_example/trueR0")#5min
cross_val01<-cross_val(counts0,r=1,dir="CV_example/trueR0")#6min
cross_val02<-cross_val(counts0,r=2,dir="CV_example/trueR0")#10.5min

r=1;set.seed(1)
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y
reorder=sample(1:n, n,replace=FALSE)
counts1=counts[reorder,]

cross_val10<-cross_val(counts1,r=0,dir="CV_example/trueR1")#5.6min
cross_val11<-cross_val(counts1,r=1,dir="CV_example/trueR1")#6.8min
cross_val12<-cross_val(counts1,r=2,dir="CV_example/trueR1")#12.8min


get_plot<-function(listmodels){
  res=data.frame(VCp=numeric(), r=numeric())
  rMax=length(listmodels)-1
  for(r in 0:rMax){
    data=listmodels[[r+1]]
    VCp=do.call(rbind, lapply(data , function(x){
      m=mean(x$vec_logpY)
      if(!is.finite(m)) m=NA
      return(m)}))
    res=rbind(res,data.frame(VCp, r=r))
  }
  data = res %>% as_tibble() %>% filter(VCp>-200000) %>% group_by(r) %>% summarize(mean.PCL=mean(VCp), sd.PCL= sd(VCp),
                                                           inf=mean.PCL - 1.96*sd.PCL/sqrt(10),
                                                           sup=mean.PCL + 1.96*sd.PCL/sqrt(10))
  g1<-data %>% ggplot( aes(y=mean.PCL,x=r))+
    geom_errorbar(aes(ymin=inf, ymax=sup), width=0,position=position_dodge(0))+
    geom_point(size=2.5)+theme_light()
  g2<-data %>% ggplot( aes(y=mean.PCL,x=r))+ geom_line()+
    geom_point(size=2.5)+theme_light()
  return(list(g1, g2, data))
}

g1<-get_plot(list(cross_val00,cross_val01,cross_val02))[[2]]+labs(title="True r=0")
g2<-get_plot(list(cross_val10,cross_val11,cross_val12))[[2]]+labs(title="True r=1")
g=grid.arrange(g1, g2, top="Selection model : simulation seed=1")
ggsave(plot=g,filename = "selecR_simu.png", path =  "/Users/raphaellemomal/these/R/images", width=5, height=5 )

#########
# Fatala
library(ade4)
data(baran95)
counts=as.matrix(baran95$fau)
p=ncol(counts); n=nrow(counts)
set.seed(1)
reorder=sample(1:n, n,replace=FALSE)
counts=counts[reorder,]
CV_fatala_0<-cross_val(counts,r=0,B=100,nblocs=10,alpha=0.05,cores=3,dir="Fatala_missing")
CV_fatala_1<-cross_val(counts,r=1,B=100,nblocs=10,alpha=0.05,cores=3,dir="Fatala_missing")
CV_fatala_2<-cross_val(counts,r=2,B=100,nblocs=10,alpha=0.05,cores=3,dir="Fatala_missing")#5h
CV_fatala_3<-cross_val(counts,r=3,B=100,nblocs=10,alpha=0.05,cores=3,dir="Fatala_missing")

#saveRDS(CV_fatala_0, file="/Users/raphaellemomal/simulations/Fatala_missing/CV_fatala_0.rds")
# saveRDS(CV_fatala_1, file="/Users/raphaellemomal/simulations/Fatala_missing/CV_fatala_1.rds")
#saveRDS(CV_fatala_2, file="/Users/raphaellemomal/simulations/Fatala_missing/CV_fatala_2.rds")
saveRDS(CV_fatala_3, file="/Users/raphaellemomal/simulations/Fatala_missing/CV_fatala_3.rds")

#-- figure
g=get_plot(list(CV_fatala_0,CV_fatala_1,CV_fatala_2,CV_fatala_3))[[2]]+
  labs(title="Fatala selection model with pairwise composite likelihood")
ggsave(plot=g,filename = "selecR_Fatala.png", path =  "/Users/raphaellemomal/these/R/images",
       width=5, height=3 )

#############
# Barents
# load("/Users/raphaellemomal/these/Data_SR/BarentsFish.Rdata")
# counts=Data$count
# p=ncol(counts); n=nrow(counts)
# reorder=sample(1:n, n,replace=FALSE)
# counts=counts[reorder,]
#countsBar=counts # trace du réarangement
CV_barents_0<-cross_val(countsBar,r=0,B=100,nblocs=10,alpha=0.05,eps=1e-2,cores=3,dir="Barents_missing_2")
CV_barents_1<-cross_val(countsBar,r=1,B=100,nblocs=10,alpha=0.05,eps=1e-2,cores=3,dir="Barents_missing_3")
CV_barents_2<-cross_val(countsBar,r=2,B=100,nblocs=10,alpha=0.05,eps=1e-2,cores=3,dir="Barents_missing_2")#5.5h
CV_barents_3<-cross_val(countsBar,r=3,B=100,nblocs=10,alpha=0.05,eps=1e-2,cores=3,dir="Barents_missing_2")
#saveRDS(CV_barents_0, file="/Users/raphaellemomal/simulations/Barents_missing/CV_barents_0.rds")
saveRDS(CV_barents_1, file="/Users/raphaellemomal/simulations/Barents_missing_3/CV_barents_1.rds")
#saveRDS(CV_barents_2, file="/Users/raphaellemomal/simulations/Barents_missing/CV_barents_2.rds")
saveRDS(CV_barents_3, file="/Users/raphaellemomal/simulations/Barents_missing/CV_barents_3.rds")
#saveRDS(countsBar, file="/Users/raphaellemomal/simulations/Barents_missing/shuffle_counts.rds")

# see what's inside
p=ncol(countsBar)
(do.call(rbind,lapply(CV_barents_0, function(x) mean(x$vec_logpY))))
(do.call(rbind,lapply(1:10, function(x){
  VEMfit=readRDS(paste0("/Users/raphaellemomal/simulations/Barents_missing/VEMfit_r1_bloc",x,".rds"))
  sd(VEMfit$M[,31])
  })))

#-- figure
g=get_plot(list(CV_barents_0,CV_barents_1,CV_barents_2,CV_barents_3))[[2]]+
  labs(title="Barents selection model with pairwise composite likelihood")
ggsave(plot=g,filename = "selecR_Barents.png", path =  "/Users/raphaellemomal/these/R/images",
       width=5, height=3 )

dataFat=get_plot(list(CV_fatala_0,CV_fatala_1,CV_fatala_2,CV_fatala_3))[[3]] %>% mutate(dat="Fatala")
#######
# pmax model r=1 barents
vec1_allY=(do.call(rbind,lapply(CV_barents_1, function(x) mean(x$vec_logpY))))
CV_barents_1Ytrain=readRDS("/Users/raphaellemomal/simulations/Barents_missing/CV_barents_1.rds")
vec1_Ytrain=(do.call(rbind,lapply(CV_barents_1Ytrain, function(x) mean(x$vec_logpY))))


dataBar=get_plot(list(CV_barents_0,CV_barents_1,CV_barents_2,CV_barents_3))[[3]] %>% mutate(dat="Barents")
#dataBar$mean.PCL[2] = mean(pmax(vec1_allY,vec1_Ytrain))
g=rbind(dataFat, dataBar) %>% 
  ggplot(aes(r, mean.PCL))+geom_point()+geom_line()+facet_wrap(~dat, scales="free")+
  labs(y="PCL")+mytheme.dark("")
ggsave(plot=g,filename = "selec_model_applis.png", path =  "/Users/raphaellemomal/these/R/images",
       width=6, height=2 )

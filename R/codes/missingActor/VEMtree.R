# VEMtree
library(EMtree)
library(PLNmodels)
library(LITree)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(mvtnorm)
library(useful)
library(mclust)
library(MASS)
library(reshape2)#for ggimage
library(gridExtra)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")

# simu parameters
set.seed(7)
n=200
p=14
r=1
type="scale-free"
plot=TRUE


# Code for one hidden covariate


################
#----- DATA
# simulate graph and omega, then sigma0 and finally counts
data=data_from_scratch(type = type,p = p+r,n = n,signed = FALSE,prob = 5/p,v = 0.001)
omega=data$omega
hidden=which(diag(omega)%in%sort(diag(omega), decreasing = TRUE)[1:r])[1:r] # on cache les r plus gros
trueClique=which(omega[hidden,-hidden]!=0)
if(plot){
  G=draw_network(1*(omega==1),groupes=1*(diag(omega)==diag(omega)[hidden][1]), 
                 layout="nicely",curv=0,nb=2,pal="black",nodes_label = 1:(p+r))$G
  print(G)
}
Kh  <- omega[hidden,hidden]
Ko  <- omega[-hidden,-hidden]
Koh <- omega[-hidden,hidden]
Km  <- Ko - Koh %*%solve(Kh)%*% t(Koh)
sigmaO=solve(Km)
counts=generator_PLN(sigmaO,covariates = NULL,n=n, seuil=15)

ome_init=omega[c(2:15,1),c(2:15,1)] #ome is for testing
diag(ome)=0
####################
#----- PLN on counts
# optimization of theta and h(Z_O)

PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M # MiO = ith row of MO
SO<-PLNfit$var_par$S # SiO = diag(ith row of SO)
sigma_obs=PLNfit$model_par$Sigma

####################
#-----  VE step
#--  Initialize

# whole Z
# initviasigma=init.mclust(sigma_obs,title="Sigma",trueClique = trueClique,n.noise=p*3+5)
# initial.param<-initEM(sigma_obs,n=n,cliquelist = list(initviasigma),pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
# omega_init=initial.param$K0
# sigma_init=initial.param$Sigma0

# Tree
O=1:p
Wg_init <- matrix(1, p+r, p+r); diag(Wg_init) = 0; Wg_init =Wg_init / sum(Wg_init)
W_init <- matrix(1, p+r, p+r); diag(W_init) = 0;# W_init =W_init / sum(W_init)
W_init[O,O] <- EMtree_corZ(cov2cor(sigma_obs),n = n,maxIter = 20)$edges_weight



VE<-function(MO,SO,sigma_obs,omega,W,Wg,MH = matrix(1,n,r),maxIter,eps, alpha,theoretical=TRUE,beta.min=1e-10, plot=FALSE){
  t1=Sys.time()
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO); H=(p+1):ncol(omega);r=length(H); iter=0 ; 
  omegaH=omega[H,H]; 
  diffKL=-1 ; diff=c(1000); Wdiff=1; diffW=c(0.1); diffMH=1; diffM=c(0.1)
  
  M<-cbind(MO, MH)
  phi=CorOmegaMatrix(omega)
  SH <- 1/omegaH
  S<-cbind(SO, rep(SH,n)) # all SHi have same solution, depending only on Mstep
  # logWtree<-computeWtree(omega, W, Wg, MH, MO, SO,phi,alpha=alpha, trim=FALSE)  
  Wg= computeWg(phi,omega,W,MH,MO,alpha,theoretical=theoretical)
  print(theoretical)
  diag(Wg)=1
  KL<-numeric(maxIter)
  
  while( ((diffKL < 0) && (iter < maxIter)) || (iter < 3)){
    iter=iter+1
    #   cat(iter)
    #-- Probabilities estimates
    Pg = EdgeProba(Wg)
    if(sum(is.nan(Pg))!=0){
      browser()
      cat( ": adjust ", summary(c(logWtree)), "\n")
      logWtree=log(adjustW(exp(logWtree)))
      cat("to: ",summary(c(logWtree)), "\n")
      Pg = EdgeProba(exp(logWtree))
    } 
    #-- Updates
    #- MH et SH
    MH.new<-  (-MO) %*% (Pg[O,H] * omega[O,H]) / omegaH
    diffMH<-max(abs(MH-MH.new))
    
    # M<-cbind(MO, MH.new)
    # Mei=Meila(Wg) #justifier que Mei soit calculée avec Wg et pas Wgtree
    # lambda=SetLambda(Pg,Mei)
    #   browser()
    Wg.new= computeWg(phi,omega,W,MH,MO,alpha) #Pg/(Mei+lambda)
    diag(Wg.new)=1
    Wg.new[which(Wg.new< beta.min)] = beta.min
    
    #  Wg.new=trimW( Wg.new) # triggers nan in Pg
    Wdiff=max(abs(F_Sym2Vec(Wg.new)-F_Sym2Vec(Wg)))
    if(is.nan(Wdiff)) browser()
    
    KL[iter]<-argminKL(F_Sym2Vec(log(Wg.new)), Pg, M,S,omega,phi,W,p)
    
    if(iter>1) diffKL =(KL[iter] - KL[iter-1])
    
    MH=MH.new
    Wg=Wg.new
    
    #-- end
    
    if(iter>1){
      diff = c(diff,diffKL)
      diffW=c(diffW,Wdiff)
      diffM=c(diffM,diffMH)
    } 
  }
  KL=KL[1:iter]
  t2=Sys.time(); time=t2-t1
  cat(paste0("\nVE step converged in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
             "\nFinal W difference: ",round(diffW[iter],5),
             "\nFinal KL difference: ",round(diff[iter],4)))
  if(plot){
    g=data.frame(Diff.W=diffW,  diff.KL=diff,diff.MH=diffM, PartofKL=KL) %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+scale_color_brewer(palette="Dark2")+
      facet_wrap(~key, scales="free")+theme_light()+labs(x="iter",y="")+
      theme(strip.background=element_rect(fill="gray50",colour ="gray50"))
    print(g)
  }
  Pg = EdgeProba(Wg)
  
  res=list(Gprobs=Pg,Gweights=Wg,Hmeans=M,Hvar=S, KL=KL, diff=diff, diffW=diffW )
  return(res)
}
resVe=VE(MO,SO,sigma_obs,ome_init,W_init,Wg_init,maxIter=150,eps=1e-3, plot=TRUE, 
         form="id1",alpha=0.6)

h=15
ome=omega[c(2:h,1),c(2:h,1)]
diag(ome)=0
seuil=0.5
performance=accppvtpr(resVe$Gprobs,ome,h,seuil)
Acc=performance[1] ;AccH=performance[2] ;AccO=performance[3] 
PPV=performance[4] ;PPVH=performance[5] ; PPVO=performance[6]
TPR=performance[7] ;TPRH=performance[8] ;TPRO=performance[9] 
p1<-ggimage(resVe$Gprobs)+labs(title=paste0("G hat (thresh=",seuil,")"))
p2<-ggimage(ome)+labs(title="G")
grid.arrange(p1,p2,ncol=2, top=paste0("Tpr=",TPR," (TprO=",TPRO," , TprH=",TPRH,
                                      ")\n Ppv=",PPV," (PpvO=",PPVO," , PpvH=",PPVH,")"))
# auc 0.92 en partant du vrai omega, 
# et en annulant le terme B dans le calcul de logWtree


####################
#-----  M steps

Mstep<-function(M,S,Pg, omega,W,maxIter, beta.min, trim=TRUE,plot=FALSE,eps, verbatim=TRUE){
  t1=Sys.time()
  diffJ=(1); diff.J=c(1);  diff.W=c(0.05)
  maxJ=c()
  n=nrow(S) 
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  iter=0
  print("Iter n°:")
  while((diffJ>eps && iter < maxIter) || iter < 3){
    iter=iter+1
    print(iter)
    #-- Updates
    # beta
    
    Mei=Meila(W)  
    lambda=SetLambda(Pg,Mei)
    W.new= Pg/(Mei+lambda)
    W.new[which(W.new< beta.min)] = beta.min
    diffW=max(abs(F_Sym2Vec(W.new)-F_Sym2Vec(W)))
    
    # W=trimW(W)
    # omega
    maxi = 1e-3
    mini = 1e+3
    
    omegaDiag <- sapply(1:(p+r), function(i){
      dichotomie(mini, maxi, function(omega_ii)
        optimDiag( omega_ii,i, omega, SigmaTilde, Pg), 1e-5)
    })
    
    omega.new=computeOffDiag(omegaDiag,SigmaTilde)
    phi=CorOmegaMatrix(omega)
    maxJ[iter]<-argmaxJ(F_Sym2Vec(log(W.new)),Pg,omega.new,SigmaTilde,phi,lambda,n)
    if(is.nan(maxJ[iter])) browser()
    if(iter>1){
      diffJ = (maxJ[iter] - maxJ[iter-1])
      diff.J = c(diff.J,diffJ)
      diff.W = c(diff.W,diffW)
    } 
    omega=omega.new
    W=W.new
  }
  t2=Sys.time(); time=t2-t1
  cat(paste0("M step converged in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
             "\nFinal maxJ difference: ",round(diff.J[iter],4)))
  # browser()
  if(plot){
    g=data.frame(Diff.W=diff.W,Diff.J=diff.J, Jbound=maxJ) %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+scale_color_brewer(palette="Dark2")+
      facet_wrap(~key, scales="free")+theme_light()+labs(x="iter",y="")+
      theme(strip.background=element_rect(fill="gray50",colour ="gray50"))
    print(g)
  }
  res=list(W=W, omega=omega, diff=diff, diffW=diffW)
  return(res)
}
resM=Mstep(M,S,Pg, ome_init,W_init,maxIter=50, beta.min=1e-6, eps=1e-2 ,plot=TRUE)


M=resVe$Hmeans
S=resVe$Hvar
Pg=resVe$Gprobs

VEMtree<-function(MO,SO,sigma_obs,ome_init,W_init,Wg_init){
  MH = matrix(10,n,r)
  resVe<-VE(MO,SO,sigma_obs,ome_init,W_init,Wg_init,MH=MH,maxIter=5,eps=1e-3, plot=TRUE, alpha=0)
  M=resVe$Hmeans
  S=resVe$Hvar
  Pg=resVe$Gprobs
  resM<-Mstep(M,S,Pg, ome_init,W_init,maxIter=5, beta.min=1e-6, eps=1e-2 ,plot=FALSE)
}


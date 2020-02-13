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
initviasigma=init.mclust(cov2cor(sigma_obs),title="Sigma",trueClique = trueClique,n.noise=p*3+5)
#findCliques(cov2cor(sigma_obs),k = 2)
initial.param<-initEM(sigma_obs,n=n,cliquelist = list(initviasigma),pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
omega_init=initial.param$K0
sigma_init=initial.param$Sigma0

# Tree
O=1:p
Wg_init <- matrix(1, p+r, p+r); diag(Wg_init) = 0; Wg_init =Wg_init / sum(Wg_init)
W_init <- matrix(1, p+r, p+r); diag(W_init) = 0;# W_init =W_init / sum(W_init)
W_init[O,O] <- EMtree_corZ(cov2cor(sigma_obs),n = n,maxIter = 20)$edges_weight



resVe.Th=VE(MO ,SO,MH = matrix(100,n,r),sigma_obs,ome_init,W_init,Wg_init,maxIter=150,minIter=3,eps=1e-3, plot=TRUE, 
            form="theory",alpha=0.9)


h=15
plotVE(resVe.Th,ome,h,seuil=0.5) 
####################
#-----  M steps

M=resVe.Th$Hmeans
S=resVe.Th$Hvar
Pg=resVe.Th$Gprobs
resM=Mstep(M,S,Pg, omega_init,W_init,maxIter=2, beta.min=1e-6, eps=1e-2 ,plot=TRUE)


VEMtree<-function(MO,SO,sigma_obs,ome_init,W_init,Wg_init, verbatim=TRUE, plot=TRUE){
  MH = matrix(100,n,r); omega=omega_init;  W=W_init;  Wg=Wg_init
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO); H=(p+1):ncol(omega);r=length(H); iter=0 ; 
  KL=c() ; J=c()
  t1=Sys.time()
  while(iter<2){
    iter=iter+1
    #VE
    resVe<-VE(MO,SO,sigma_obs,omega,W,Wg,MH=MH,maxIter=5,minIter=3,eps=1e-3, plot=FALSE, 
              form="theory",alpha=0.9, verbatim=FALSE)
    KL[iter]=resVe$KL
    M=resVe$Hmeans ; MH=matrix(M[,H],n,r)
    S=resVe$Hvar
    Pg=resVe$Gprobs
    Wg=resVe$Gweights
    #M
    resM<-Mstep(M,S,Pg, ome_init,W_init,maxIter=5, beta.min=1e-6, eps=1e-2 ,plot=FALSE, verbatim=FALSE)
    W=resM$W
    omega=resM$omega
    J[iter]=resM$finalJ
  }
  t2=Sys.time()
  t2=Sys.time(); time=t2-t1
  if(verbatim) cat(paste0("\nVEMtree ran in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal Jbound difference: ",round(J[iter]-J[iter-1],5)))
  if(plot){
    g=data.frame(Jbound=J, KL=KL) %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+scale_color_brewer(palette="Dark2")+
      facet_wrap(~key, scales="free")+theme_light()+labs(x="iter",y="")+
      theme(strip.background=element_rect(fill="gray50",colour ="gray50"))
    print(g)
  }
  
  return(list(M=M,S=S,Pg=Pg,Wg=Wg,W=W,omega=omega))
}
resVEM<-VEMtree(MO,SO,sigma_obs,omega_init,W_init,Wg_init)
plotVEM(resVEM$Pg,ome,r=1,seuil=0.5)
values=courbes_seuil(probs = resVEM$Pg,omega = ome,h = 15,seq_seuil = seq(0,1,0.05))
values %>%as_tibble() %>%    gather(key,value, -seuil) %>% 
  ggplot(aes(seuil,value,color=key))+
  geom_point()+  geom_line()+
  facet_wrap(~key, scales="free")+theme_light()+scale_color_brewer(palette="Paired")+
  theme(strip.background=element_rect(fill="gray50",colour ="gray50"))

#TODO
# choice of alpha
# stop criterion
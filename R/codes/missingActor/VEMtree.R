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
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(harrypotter)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")

mytheme.light <- list(theme_light(), scale_color_brewer("",palette="Set3"),guides(color=FALSE),
                      theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                            plot.title = element_text(hjust = 0.5)))

mytheme.dark <- list(theme_light(), scale_color_brewer("",palette="Dark2"),guides(color=FALSE),
                     theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                           plot.title = element_text(hjust = 0.5)))
# mypal<-c(brewer.pal(3, "Blues"),brewer.pal(3, "Reds"),brewer.pal(3, "Greens"))
# mypal<-c(brewer.pal(8, "Dark2"),"blue","red")
mytheme <- list(theme_light(), scale_fill_brewer("",palette="Dark2"),#scale_colour_hp_d(option = "RavenClaw", name = ""), #scale_color_manual("",values=mypal),
                guides(color=FALSE),
                theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                      plot.title = element_text(hjust = 0.5)))

# simu parameters
set.seed(7)
n=200
p=14
r=1
type="scale-free"
plot=TRUE

################
#----- DATA
# simulate graph and omega, then sigma0 and finally counts
data=data_from_scratch(type = type,p = p+r,n = n,signed = FALSE,prob = 5/p,v = 0)
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
counts=generator_PLN(sigmaO,covariates = NULL,n=n)

ome_init=omega[c(2:15,1),c(2:15,1)] #ome is for testing
ome=ome_init
diag(ome)=0
####################
#----- PLN on counts
# optimization of theta and h(Z_O)

PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M # MiO = ith row of MO
SO<-PLNfit$var_par$S # SiO = diag(ith row of SO)
sigma_obs=PLNfit$model_par$Sigma

plot(sigmaO ,sigma_obs)
####################
#-----  VE step
#--  Initialize

# whole Z
initviasigma=init.mclust((cov2cor(sigma_obs)),title="Sigma",
                         trueClique = trueClique,n.noise=1)$init

# etude de la quantité de bruit à privilégier:
# B=30 ; coeff=seq(0.1,4,0.5)
# test=lapply(coeff,function(c){
#   res=cbind( data.frame( t(sapply(1:B, function(x){
#     init.mclust(sigma_obs,title="Sigma",
#                 trueClique = trueClique,n.noise=round(p*c))$FPN
#   }))),c)
# })
# test=do.call(rbind,test)
# colnames(test)=c("FN","FP","coeff")
# test %>% mutate(nul = (FN==1 & FP==0)) %>% filter(!nul) %>% dplyr::select(-nul) %>% 
#   gather(FNP,value,-coeff) %>% 
#   ggplot(aes(as.factor(coeff),value,fill=FNP))+geom_boxplot()+mytheme.light

# test %>% mutate(nul = (FN==1 & FP==0)) %>% filter(!nul) %>% dplyr::select(-nul) %>%
#   group_by(coeff) %>% summarise(sumFP=sum(FP), sumFN=sum(FN), sum=sum(FP+FN))

#findCliques(cov2cor(sigma_obs),k = 2)
initial.param<-initEM(sigma_obs,n=n,cliquelist = list(initviasigma),pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
omega_init=initial.param$K0
sigma_init=initial.param$Sigma0

# Tree
O=1:p
Wg_init <- matrix(1, p+r, p+r); diag(Wg_init) = 0; Wg_init =Wg_init / sum(Wg_init)
W_init <- matrix(1, p+r, p+r); diag(W_init) = 0;# W_init =W_init / sum(W_init)
W_init[O,O] <- EMtree_corZ(cov2cor(sigma_obs),n = n,maxIter = 20)$edges_weight


pr=prcomp(t(counts),scale. = FALSE)
MH = matrix(pr$rotation[,1]*pr$sdev[1],nrow=n,ncol=r)
resVe.Th=VE(MO ,SO,MH = matrix(100,n,r),sigma_obs,ome_init,W_init,Wg_init,maxIter=150,
            minIter=8,eps=1e-3, plot=TRUE,
            form="theory",alpha=0.9)
h=15
plotVE(resVe.Th,ome,h,seuil=0.5)
####################
#-----  M steps

# M=resVe.Th$Hmeans
# S=resVe.Th$Hvar
# Pg=resVe.Th$Gprobs
# resM=Mstep(M,S,Pg, omega_init,W_init,maxIter=2, beta.min=1e-6, eps=1e-2 ,plot=TRUE)


VEMtree<-function(counts,MO,SO,sigma_obs,ome_init,W_init,Wg_init, verbatim=TRUE,maxIter=20, 
                  plot=TRUE, eps=1e-2, alpha, vraiOm){
  
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO); H=(p+1):ncol(omega);r=length(H);
  pr=prcomp(t(counts),scale. = FALSE)
  MH = matrix(pr$rotation[,1]*pr$sdev[1],nrow=n,ncol=r)
  omega=omega_init;  W=W_init;  Wg=Wg_init
  iter=0 ; lowbound=list()
  KL=c() ; J=c();diffW=c();diffOm=c();diffWg=c();diffPg=c();diffMH=c();diffOmDiag=c();rvalue=c();diffquantile=c()
  t1=Sys.time()
  Pg=matrix(0.5, ncol(W),ncol(W))
  SH=matrix(1,nrow=n,ncol=r)
  while( iter<maxIter){ #(diffW[iter] > eps && diffOm[iter] > eps && iter < maxIter) || iter<1
    
    iter=iter+1 
    cat(paste0("\n Iter n°", iter))
    #VE
    resVe<-VE(MO,SO,SH, sigma_obs,omega,W,Wg,MH=MH,Pg=Pg,maxIter=1,minIter=1,eps=1e-3, plot=FALSE, 
              form="theory",alpha=alpha, verbatim=FALSE)
    KL[iter]=resVe$KL
    M=resVe$Hmeans ; 
    S=resVe$Hvar ; SH=matrix(S[,H],n,r)
    Pg.new=resVe$Gprobs
    Wg.new=resVe$Gweights
    MH.new<-matrix(M[,H],n,r)

    diffMH[iter]<-abs(max(MH.new-MH))
    diffWg[iter]<-abs(max(Wg.new-Wg))
    diffPg[iter]<-abs(max(Pg.new-Pg))
    Wg=Wg.new
    Pg=Pg.new
    MH=MH.new
    #M
    resM<-Mstep(M,S,Pg, omega,W,maxIter=5, beta.min=1e-6, eps=1e-3 ,plot=FALSE, verbatim=FALSE,
                Wg=Wg, p=p)
    W.new=resM$W
    diffW[iter]=abs(max(W.new-W))
    W=W.new
    omega.new=resM$omega
    diffOm[iter]=abs(max(F_Sym2Vec(omega.new)-F_Sym2Vec(omega)))
    diffOmDiag[iter]=abs(max(diag(omega.new)-diag(omega)))
    omega=omega.new
    rvalue[iter]=summary(lm(diag(omega)~diag(vraiOm)))$r.squared
    diffquantile[iter]=quantile(F_Sym2Vec(omega)[F_Sym2Vec(vraiOm)==1], 0.25) - quantile(F_Sym2Vec(omega)[F_Sym2Vec(vraiOm)==0], 0.75)

    J[iter]=resM$finalJ
    
    
    lowbound[[iter]] = LowerBound(Pg ,omega, M, S, W, Wg,p)
  }
  lowbound=do.call(rbind,lowbound)
  features<-data.frame(diffMH=diffMH, diffWg=diffWg, diffPg=diffPg, diffW=diffW, diffOm=diffOm, diffOmDiag=diffOmDiag,
                       adjustDiag=rvalue, diffquantile=diffquantile)
  t2=Sys.time()
  t2=Sys.time(); time=t2-t1
  if(verbatim) cat(paste0("\nVEMtree ran in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal Jbound difference: ",round(J[iter]-J[iter-1],5)))
  if(plot){
    g1<-features %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line() + facet_wrap(~key, scales="free")+
      labs(x="",y="", title="Stoping criteria")+ mytheme.dark
    
    
    g2<- lowbound %>% rowid_to_column() %>% gather(key,value,-rowid) %>% 
      ggplot(aes(rowid,value, color=key))+geom_point()+geom_line()+
      facet_wrap(~key, scales="free")+
      labs(x="iteration",y="", title="Lower bound and components")+mytheme
    
    grid.arrange(g1, g2, ncol=1)
  }
  
  
  return(list(M=M,S=S,Pg=Pg,Wg=Wg,W=W,omega=omega, lowbound=lowbound, features=features))
}

resVEM<-VEMtree(counts,MO,SO,sigma_obs,omega_init,W_init,Wg_init, eps=1e-3, alpha=1,
                maxIter=10, plot=TRUE,vraiOm=ome_init)
plotVEM(resVEM$Pg,ome,r=1,seuil=0.5)
values=courbes_seuil(probs = resVEM$Pg,omega = ome,h = 15,seq_seuil = seq(0,1,0.05))
plotVerdict(values, seuil)

# lower bound check
resVEM$lowbound %>% rowid_to_column() %>% gather(key,value,-rowid) %>% 
  ggplot(aes(rowid,value, color=key))+geom_point()+geom_line()+
  facet_wrap(~key, scales="free")+mytheme



#====
#ajustement omega

# vrai vs initialisation vs estimation
# ome_init = vrai omega reordonne pour comparaison (acteur manquant en fin de matrice)
# omega_init = matrice d'initialisation obtenue par initEM et mclust
# omega_res = omega estime par VEMtree
omegas = data.frame(vrai = F_Sym2Vec(ome_init) , init = F_Sym2Vec(omega_init),
                    estimation = F_Sym2Vec(resVEM$omega))
omegas %>% gather(key, value, -vrai) %>% 
  ggplot(aes(as.factor(vrai), value,  fill=key, clor=key))+geom_boxplot()+
  mytheme#+coord_cartesian(ylim = c(-1,2))

omegas %>% group_by(vrai) %>% summarise(q75=quantile(estimation, 0.75),
                                        q25=quantile(estimation, 0.25))
quantile(omegas$estimation[omegas$vrai==1], 0.25) - quantile(omegas$estimation[omegas$vrai==0], 0.75)

Diagomegas = data.frame(vrai = diag(ome_init) , init = diag(omega_init),
                        estimation = diag(resVEM$omega))
Diagomegas %>% gather(key, value, -vrai) %>% 
  ggplot(aes((vrai), value,  color=key))+geom_point()+ theme_light()+
  geom_abline()

rvalue=summary(lm(Diagomegas$estimation~Diagomegas$vrai))$r.squared

#TODO
# choice of alpha
# VEM stop criterion
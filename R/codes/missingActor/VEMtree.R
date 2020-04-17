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

# test epsilon SNR
# seqEpsi=seq(1.01, 5, 0.05)
# varyEpsi=sapply(seqEpsi, function(epsilon){
#   initial.param<-initEM(sigma_obs,n=n,cliqueList = list(initviasigma),cst=epsilon, pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
#   omega_init=initial.param$K0
#   compute_nSNR(K = omega_init,indexmissing = ncol(omega_init))
# })
# data.frame(SNR=varyEpsi,epsilon=seqEpsi) %>% ggplot(aes(epsilon, SNR))+geom_point()+mytheme

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
library(parallel)
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(harrypotter)
library(sparsepca)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")
#source("/home/mmip/Raphaelle/these/R/codes/missingActor/fonctions-missing.R")

# simu parameters
set.seed(7)
n=200 ;p=14;r=1;type="scale-free";plot=TRUE
O=1:p
################
#----- DATA
# simulate graph and omega, then sigma0 and finally counts
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y
sigmaO= missing_data$Sigma
omega=missing_data$Omega
trueClique=missing_data$TC
hidden=missing_data$H

PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M  
SO<-PLNfit$var_par$S  
sigma_obs=PLNfit$model_par$Sigma
theta=PLNfit$model_par$Theta 
matcovar=matrix(1, n,1)
order_sigma<-sigma_obs[c(trueClique[[1]],setdiff(1:p,trueClique[[1]])),
                       c(trueClique[[1]],setdiff(1:p,trueClique[[1]]))]
ggimage(order_sigma)
ggimage(solve(order_sigma))
plot(sigmaO ,sigma_obs)
ome_init=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), hidden)]
ome=ome_init ; diag(ome)=0


# alpha<-tryCatch(expr={computeAlpha(solve(sigma_obs),default =0.3, MO, SO, plot=plot)},
#                 error=function(e){message("sythetic alpha")
#                   return(0.3)})
# clique_mclust=init.mclust((cov2cor(sigma_obs)),title="Sigma", nb.missing = r,
#                           trueClique = trueClique,n.noise=1)
# clique_mclust$init
# plotInitMclust(res=clique_mclust,title = "")


cliques_spca <- boot_FitSparsePCA(scale(MO),100,r=2, cores=1)

# best VEM with 1 missing actor
ListVEM<-List.VEM(cliqueList=cliques_spca, counts, sigma_obs, MO,SO,r=1,eps=1e-3, cores=3,maxIter=100)

vBICs<-(criteria(ListVEM,counts,theta, matcovar,r))
vBICs$J
best=which.max(vBICs$J)
VEM_1<-ListVEM[[best]]
computeFPN(VEM_1$clique,trueClique = trueClique[[1]],p=p) 
VEM_1$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-V6) %>% 
  ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(V6)), size=3)+geom_line()+
  facet_wrap(~key, scales="free")+
  labs(x="iteration",y="", title="Lower bound and components")+mytheme+
  scale_color_discrete("")
plotVEM(VEM_1$Pg,ome,r=1, 0.5)
 
# TJ=True_lowBound(Y=counts, M=VEM_1$M,VEM_1$S,theta=theta,X=matcovar,
#               W=VEM_1$W,Wg=VEM_1$Wg,Pg=VEM_1$Pg, omega=VEM_1$omega)
# ICL_T(TJ, Pg=VEM_1$Pg, Wg=VEM_1$Wg)
# crit=criteria(ListVEM,counts=counts,theta = theta,matcovar = matcovar, r=r)
# crit%>%
#   gather(key, value) %>%
#   ggplot(aes(key, value, color=key, fill=key))+geom_boxplot(alpha=0.3)+
#   mytheme.dark("")


####################
#-----  VEM
# find alpha on the observed part of the initial "non beta" quantities needed to compute the beta tilde
# alpha<-computeAlpha(omegainit[O,O], MO, SO)
# alpha=1
init=initVEM(counts = counts,initviasigma=cliques_spca$cliqueList[[1]], sigma_obs,r = 2) #cliques_spca$cliqueList[[4]]
Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
test=omegainit ;diag(test)=0
ggimage(test)
r=2
q=p+r
D=.Machine$double.xmax
 
alpha = (1/n)*((1/(q-1))*log(D) - log(q))
#alpha2=(1/(n*q))*log(D/(q^(q/2)))
# curve((1/n)*((1/(x-1))*log(D) - log(x)),from=15, to=30)
# curve((1/(n*x))*log(D/(x^(x/2))),from=15, to=30, add=T, col="red")

resVEM<-VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, eps=1e-5, alpha=alpha, 
                maxIter=100, plot=TRUE,print.hist=FALSE, verbatim = TRUE,nobeta = TRUE)

plotVEM(resVEM$Pg,ome,r=1,seuil=0.5)
ggimage(resVEM$Pg)
resVEM$features
g1<-ggimage(resVEM$Pg)+labs(title="après la bosse")
g2<-ggimage(resVEM$Pg)+labs(title="avant la bosse")
ggimage(resVEM$omega)
ggimage(EdgeProba(resVEM$Wg))
ggimage(EdgeProba(resVEM$W))

g3<-ggimage(ome)+labs(title="vérité")

grid.arrange(g2, g1, g3, nrow=1)
 

valuesRaw=courbes_seuil(probs = resVEMraw$Pg,omega = ome,h = 15,seq_seuil = seq(0,1,0.02))
valuesRaw %>% mutate(crit = PPV+TPR) %>% filter(crit==max(crit, na.rm=TRUE)) %>%
  summarise(mins=min(seuil), maxs=max(seuil), PPV=max(PPV), PPVO=max(PPVO),PPVH=max(PPVH), 
            TPR=max(TPR),TPRO=max(TPRO),TPRH=max(TPRH))
plotVerdict(valuesFilter, seuil)+guides(color=FALSE)
 
#ggsave("precrec_missing.png", plot=p, width=8, height=4,path= "/Users/raphaellemomal/these/R/images")

# lower bound check
ListVEM[[1]]$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-V6) %>% 
  ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(V6)), size=3)+geom_line()+
  facet_wrap(~key, scales="free")+
  labs(x="iteration",y="", title="Lower bound and components")+mytheme+
  scale_color_discrete("")

# model selection
PLNfit = PLN(counts~1)
theta=PLNfit$model_par$Theta
int=matrix(1, nrow(counts),1) # intercept
matcovar=int
omegainit0=solve(PLNfit$model_par$Sigma)
MO=PLNfit$var_par$M
SO=PLNfit$var_par$S
VEM0<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit0,W_init = W_init[O,O],Wg_init = Wg_init[O,O],
              plot = FALSE, maxIter = 10,print.hist = FALSE, vraiOm = NULL, alpha=1, verbatim=FALSE )

VEM1<-VEMtree(counts,MO = MO,SO = SO,MH=MHinit,ome_init = omega_init,W_init = W_init,Wg_init = Wg_init,
              eps=1e-2, alpha=0.9, maxIter=10, plot=FALSE,vraiOm = NULL, print.hist = FALSE)

J0<-True_lowBound(counts,VEM0$M,VEM0$S, theta, matcovar,VEM0$W, VEM0$Wg, VEM0$Pg, VEM0$omega )
vBIC0<-VBIC(J0, ncol(counts),r=0, d=1, n=nrow(counts))

J1<-True_lowBound(counts,VEM1$M,VEM1$S, theta, matcovar,VEM1$W, VEM1$Wg, VEM1$Pg, VEM1$omega )
vBIC1<-VBIC(J1, ncol(counts),r=1, d=1, n=nrow(counts))


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
Diagomegas[-nrow(Diagomegas),] %>% gather(key, value, -vrai) %>% 
  ggplot(aes((vrai), value,  color=key))+geom_point()+ theme_light()+
  geom_abline()

rvalue=summary(lm(Diagomegas$estimation~Diagomegas$vrai))$r.squared

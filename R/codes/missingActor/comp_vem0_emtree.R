library(EMtree)
library(PLNmodels)
library(tidyverse)
library(mvtnorm)
library(gridExtra)
library(reshape2)
library(ROCR)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")

# comparaison VEMtree0 et EMtree
seed=1 ; r=0 ; p=14
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
ome=omega ; diag(ome)=0
#VEM with no missing actor
Wginit <- matrix(1, p+r, p+r); Wginit =Wginit / sum(Wginit)
Winit <- matrix(1, p+r, p+r); Winit =Winit / sum(Winit)
Winit[1:p,1:p] <- EMtree_corZ(cov2cor(sigma_obs),n = n,maxIter = 20,verbatim = FALSE)$edges_weight
diag(Wg_init) = 1;diag(W_init) = 1
omegainit=solve(sigma_obs)
alpha<-tryCatch(expr={computeAlpha(omegainit,default =0.3, MO, SO, plot=TRUE)},
                error=function(e){message("sythetic alpha")
                  return(0.3)})
VEM_0<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,
               Wg_init =Wginit,plot = TRUE,eps=1e-4, maxIter = 10,print.hist = FALSE,
               vraiOm = NULL, alpha=1/n, verbatim=TRUE, filterPg = TRUE, filterWg=TRUE )
plotVEM(VEM_0$Pg,ome,r=0,seuil=0.5)

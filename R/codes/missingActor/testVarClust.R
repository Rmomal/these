
rm(list=ls()); par(mfrow=c(1, 1))
library(mvtnorm); library(sna)
setwd("/Users/raphaellemomal/these/Pgm_SR")
source('Functions/FunctionVarClustPCA.R')
library(PLNmodels)
n = 200; p = 14; seed = 10; set.seed(seed) 
r=1;type="scale-free" 
# exemple fétiche

missing_data<-missing_from_scratch(n,p,r,type,plot=FALSE)
counts=missing_data$Y
sigmaO= missing_data$Sigma
omega=missing_data$Omega
trueClique=missing_data$TC
hidden=missing_data$H

PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M  
SO<-PLNfit$var_par$S  
sigma_obs=PLNfit$model_par$Sigma

vcPCA = F_VarClustPCA(sigma_obs)
filter_content<-list()
c=1
lapply(seq_along(vcPCA$clustContent), function(x){
  if(length(vcPCA$clustContent[[x]])>2){
    filter_content[[c]]<<-vcPCA$clustContent[[x]]
    c<<-c+1
  } 
})

computeFPN(filter_content,trueClique[[1]],p)

clique_mclust=init.mclust((cov2cor(sigma_obs)),title="Sigma", nb.missing = r,
                          trueClique = trueClique,n.noise=3*p)
clique_mclust$init
plotInitMclust(res=clique_mclust,title = "")
init=initVEM(counts = counts, trueClique = trueClique,initviasigma=clique_mclust$init, sigma_obs,r = r)
Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit

plot(sigmaO ,sigma_obs)
ome_init=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), hidden)]
ome=ome_init ; diag(ome)=0

computeFPN(res = clique_mclust,trueClique = trueClique)
####################


S = cov2cor(Snone)

library(EMtree)
library(PLNmodels)
library(mvtnorm)
library(tidyverse)
nbgraph=10
valeur=20
type="cluster"
param="d"
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"

#needs data_from_stored_graphs in Kit
#needs generator_PLN in Kit


test_EMtree<-function(type, param, nbgraph, valeur){
  data=data_from_stored_graphs(type,param,nbgraph,valeur,path=path, covar=NULL, fixe=TRUE)[[1]]
  covar<-readRDS(paste0(path, "fixeCovar.rds"))
  m<- model.matrix(~X1+X2+X3,covar)[,-1]
  model<-PLN(data ~ m,control = list("trace"=0))
  inf<-EMtree(model,  maxIter=30, cond.tol=1e-15, verbatim=FALSE, plot=FALSE)
  return(inf)
}


test_EMtree(type, param, nbgraph, valeur)

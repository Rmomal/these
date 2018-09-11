#############
# Required, Sources & Functions
#############
library(PLNmodels)
library(mvtnorm)

source('/home/momal/Git/these/pack1/R/codes/FunctionsMatVec.R')
source('/home/momal/Git/these/pack1/R/codes/FunctionsTree.R')
source('/home/momal/Git/these/pack1/R/codes/FunctionsInference.R')
source('/home/momal/Git/these/pack1/R/codes/TreeMixture-RML.R')

generator_PLN<-function(Sigma,covariates){
  # vraies abondances, log normales
  n<-nrow(covariates)
  c<-ncol(covariates)# nb covariables
  p<-ncol(Sigma) # nb esp??ces
  beta<-matrix(runif(c*p),c,p)
  # browser()
  Z<- rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
  Y = matrix(rpois(n*p, exp(Z+ covariates %*% beta)), n, p)
  return(list(Y,cor(Z)))
}

##############
##############

from_sigma_x<-function(sigma,covariates){
  Y<-generator_PLN(sigma,covariates)[[1]]
  PLN = PLN(Y ~ -1+covariates)                  # run PLN
  Sigma<-PLN$model_par$Sigma
  corVEM<-cov2cor(Sigma)
  inf_treeggm<-TreeGGM(corVEM,"FALSE",FALSE)$P  # run EM
  return(inf_treeggm)
}

from_stored_graphs<-function(type, variable, nbgraph, valeur){
  path<-"~/these/pack1/R/Simu/PLNcov/"
  if(variable!="n") {
    param<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds"))
  }else{
    param<-readRDS(paste0(path,type,"/n/Sets_param/Graph",nbgraph,".rds"))
  }
  covariates<-cbind(rep(c(0,1),each=n/2),rnorm(n,8,0.5),
                     round(runif(n)*10))
  Y<-generator_PLN(param$sigma,covariates)[[1]]
  PLN = PLN(Y ~ -1+covariates)
  Sigma<-PLN$model_par$Sigma
  corVEM<-cov2cor(Sigma)
  inf_treeggm<-TreeGGM(corVEM,"FALSE",FALSE)$P
  return(inf_treeggm)
}




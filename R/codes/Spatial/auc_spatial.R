
#########
# Simu spatialized PLN data with latent structure
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
library(PLNmodels)
library(mvtnorm)
library(EMtree)
library(ROCR)
source('/Users/raphaellemomal/simulations/codes/FunctionsMatVec.R')
source('/Users/raphaellemomal/simulations/codes/FunctionsTree.R')
source('/Users/raphaellemomal/simulations/codes/FunctionsInference.R')
source('/Users/raphaellemomal/simulations/codes/fonctions.R')
### covariance functions
lapply(c("sp", "mvtnorm", "tidyverse", "fields", "raster"), library, character.only = TRUE)

cov_exp <- function(distance, range = 1, sill = 1) {
  sill * sill * exp(- distance / range)
}
cov_matern <- function(distance, range = 1, sill = 1) {
  sill * sill * (1 + distance * sqrt(3) / range) * exp(-distance * sqrt(3) / range)
}

gener_Delta<-function(n){
  sqrt_n <- sqrt(n)
  grid <- data.frame(lon = rep(seq(0, 1, length.out = sqrt_n), times = sqrt_n), 
                     lat = rep(seq(0, 1, length.out = sqrt_n), each = sqrt_n)
  )
  Delta <- as.matrix(cov_matern(distance = dist(grid, method = "euclidean")))
  diag(Delta) <- 1
  return(Delta)
}



generator_spatial_PLN<-function(Sigma,covariates){
  # vraies abondances, log normales
  n<-nrow(covariates)
  p<-ncol(Sigma) # nb esp??ces
  m<- model.matrix(~X1+X2+X3,covariates)[,-1]
  mc<-ncol(m)
  beta<-matrix(runif(p*mc),mc,p)
  
  Delta<-gener_Delta(n)
  Z<- rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
  Epsi<- t(rmvnorm(p, rep(0, nrow(Delta)), Delta))
  
  Y = matrix(rpois(n*p, exp( m %*% beta + Z+ Epsi)), n, p)
  return(list(Y,cor(Z)))
}

spatial_data_from_stored_graphs<-function(type, variable, nbgraph, valeur, covar,path, fixe=TRUE){
  
  if(variable!="n") {
    param<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds"))
    n=nrow(covar)
    if (fixe){
      covar<-readRDS(paste0(path, "fixeCovar.rds"))
    }else{
      covar=covar
    }
  }else{
    covar<-readRDS(paste0(path, "fixeCovar_n",valeur,".rds"))
    param<-readRDS(paste0(path,type,"/n/Sets_param/Graph",nbgraph,".rds"))
    n=valeur
  }
  
  gener<-generator_spatial_PLN(param$sigma,covar) #gener = list(Y,corZ)
  return(gener)
}

###################
nbgraph=10
valeur=20
type="cluster"
variable="d"

get_auc<-function(type, param, nbgraph, valeur, sp_noise=FALSE){
  if(sp_noise){
    data=spatial_data_from_stored_graphs(type,param,nbgraph,valeur,path=path, covar=NULL, fixe=TRUE)[[1]]
  }else{
    data=data_from_stored_graphs(type,param,nbgraph,valeur,path=path, covar=NULL, fixe=TRUE)[[1]]
  }
  covar<-readRDS(paste0(path, "fixeCovar.rds"))
  m<- model.matrix(~X1+X2+X3,covar)[,-1]
  model<-PLN(data ~ m,control = list("trace"=0))
  inf<-EMtree(model,  maxIter=30, cond.tol=1e-15, verbatim=FALSE, plot=FALSE)
  # obs<-readRDS(paste0(path,type,"/",param,"/Sets_param/Graph",nbgraph,"_",valeur,".rds"))$omega
  # pred<-inf$ProbaCond
  # auc=diagnost_auc(obs,pred)
  return(inf)
}
B=200
#needs data_from_stored_graphs in Kit
#needs generator_PLN in Kit
#needs diagnost_auc in function
#needs vec_obs_pred in function
#PLN and ROCR packages
res=replicate(B, get_auc(type, param, nbgraph, valeur, sp_noise=FALSE))
res_noise=replicate(B, get_auc(type, param, nbgraph, valeur, sp_noise=TRUE))
tibble(brut=res, noise=res_noise) %>% 
  gather(key, value) %>% 
  ggplot(aes(key, value))+geom_beeswarm(aes(color=key))

################
data=spatial_data_from_stored_graphs(type,param,nbgraph,valeur,path=path, covar=NULL, fixe=TRUE)[[1]]
covar<-readRDS(paste0(path, "fixeCovar.rds"))
m<- model.matrix(~X1+X2+X3,covar)[,-1]
model<-PLN(data ~ m,control = list("trace"=0))

sigma=model$model_par$Sigma
M=model$var_par$M
S=model$var_par$S

sigHat=(t(M)%*%M)/100+diag(colMeans(S))
plot(sigma, sigHat)




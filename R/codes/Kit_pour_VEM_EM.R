#############
# Required, Sources & Functions
#############
library(PLNmodels)
library(mvtnorm)
source('/Users/raphaellemomal/simulations/codes/FunctionsMatVec.R')
source('/Users/raphaellemomal/simulations/codes/FunctionsTree.R')
source('/Users/raphaellemomal/simulations/codes/FunctionsInference.R')
source('/Users/raphaellemomal/simulations/codes/fonctions.R')
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
# n=100
# fixeCovar<-data.frame(X1=round(runif(n)*5),X2=rnorm(n,0,2),
#                   X3=as.character(round(runif(n)*5)))
# saveRDS(fixeCovar,paste0(path, "fixeCovar.rds"))
# for(n in 70){
#   fixeCovar<-data.frame(X1=round(runif(n)*5),X2=rnorm(n,0,2),
#                         X3=as.character(round(runif(n)*5)))
#   saveRDS(fixeCovar,paste0(path, "fixeCovar_n",n,".rds"))
# }
fixeCovar<-readRDS(paste0(path, "fixeCovar.rds"))
##############
# DATA
##############
generator_PLN<-function(Sigma,covariates){
  # vraies abondances, log normales
  n<-nrow(covariates)
  c<-ncol(covariates)# nb covariables
  p<-ncol(Sigma) # nb esp??ces
  m<- model.matrix(~X1+X2+X3,covariates)[,-1]
  mc<-ncol(m)
  beta<-matrix(runif(p*mc),mc,p)

  Z<- rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
  Y = matrix(rpois(n*p, exp(Z+ m %*% beta)), n, p)
  return(list(Y,cor(Z)))
}

data_from_stored_graphs<-function(type, variable, nbgraph, valeur, covar,path, fixe=TRUE){
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

  gener<-generator_PLN(param$sigma,covar) #gener = list(Y,corZ)
  return(gener)
}

data_from_scratch<-function(type, p=20, r=5, covar,prob=log(p)/p,dens=log(p)/p){
  #print(p)
  # make graph
  #browser()
  graph<- generator_graph(graph=type,d=p,prob=prob,dens=dens,r=r)

  param<-generator_param(as.matrix(graph))
  data<-generator_PLN(param$sigma,covar)[[1]]
  # as_tbl_graph(as.matrix(graph)) %>%
  #   ggraph(layout="kk")+
  #   geom_edge_link()+
  #   geom_node_point(size=3, color="blue")
  #
  return(list(data=data,omega= param$omega))
}
data_for_MInt<-function(Y,covar,path){ # à voir si on construit Y avant ou appel à data_from_stored_graphs
  Y <-cbind(1:nrow(Y),Y)
  Y<-rbind(c("Observations",1:(ncol(Y)-1)),Y)

  covariates <-cbind(1:nrow(covar),covar)
  #  browser()
  # covariates<-rbind(c("Observations","feature1","feature2","feature3"),covariates)
  covariates<-rbind(c("Observations","feature1","feature2"),covariates)

  pathY<-paste0(path,"mint_data/y.txt")
  pathX<-paste0(path,"mint_data/x.txt")
  write.table(Y, file = pathY, sep = " ", col.names = FALSE, row.names = FALSE)
  write.table(covariates, file = pathX, sep = " ", col.names = FALSE, row.names = FALSE)

  invisible(list(y=pathY,x=pathX))
}

##############
# INFERENCES
##############
from_sigma_x<-function(sigma,covariates){
  Y<-generator_PLN(sigma,covariates)[[1]]
  PLN = PLN(Y ~ -1+covariates)                  # run PLN
  Sigma<-PLN$model_par$Sigma
  corVEM<-cov2cor(Sigma)
  inf_treeggm<-TreeGGM(corVEM,"FALSE",FALSE)$P  # run EM
  return(inf_treeggm)
}

from_stored_graphs<-function(type, variable, nbgraph, valeur,covar=fixeCovar,step="FALSE", cond.tol=1e-10,path, maxIter=300){
  maxIter<-ifelse(step=="TRUE",1,maxIter)
  #Y<-data_from_stored_graphs(type, variable, nbgraph, valeur,covar,path)[[1]]
  Y<-readRDS(paste0(path,type,"/",variable,"/YData/Y_",nbgraph,"_",valeur,".rds"))[[1]]
  T1<-Sys.time()
  m<- model.matrix(~X1+X2+X3,covar)[,-1]
  model<-PLN(Y ~ m,control = list("trace"=0))
  #inf<-TreeGGM(cov2cor(Sigma),step=step,maxIter = 150, n=nrow(Y), cond.tol= cond.tol)
  inf<-EMtree(model,  maxIter=maxIter, cond.tol=cond.tol, verbatim=FALSE, plot=FALSE)
  T2<-Sys.time()

  time<-difftime(T2,T1)

  return(list(inf,timeFit=time))
}


#########################################################
#########################################################

roc_curve<-function(pred, obs){
  nvar<-ncol(obs)
  obs[which(abs(obs)<1e-16)]<-0
  indices_nuls<-which(obs==0)
  label<-matrix(1,nrow=nrow(obs),ncol=ncol(obs))
  label[indices_nuls]<-0
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),
                         as.vector(label[upper.tri(label)]))

  perf <- performance( prediction, "tpr", "fpr" )
  ROC_auc <- performance(prediction,"auc")
  title=round(ROC_auc@y.values[[1]],digits=3)
  # plot( perf, main=title )
  return(title)
}
# result<-sapply(1:100, function(x){
#   res<-from_stored_graphs("erdos","d",x,20)
#   c(roc_curve(res[[1]], res[[3]]),
#     roc_curve(res[[2]], res[[3]]))
# })
# median(result[1,]-result[2,])


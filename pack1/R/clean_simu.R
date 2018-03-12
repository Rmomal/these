
#############################
# The Code Idea :
# Here we initiate a graph with edge probability prob and number of nodes d, and compute the 
# parameters  sigma and K (with noise u) corresponding to a gaussian graphical model version 
# of the network and save their values.
# Then we use the value of sigma to infer the network according to three methods : glasso,
# treeGGM and treeGGM1step. Using the truth given by K, we compare the methods performances
# according to AUC, sensibility and sensitivity criteria. Each comparison is done B times
# and a graphical summary is given in the "images" folder, giving the mean and the IC95 of the latter
# for each method, each criteria and at each combination of tested parameters.

# In order to compare methods and the impact of the parameters u (noise on  omega=K),
# n (number of observations),prob (edge probability) and d (number of vertices),
# ** we have to notice that prob and d affect the shape of the graph, so we have one graph
# per value of these variables, and then one file of saved parameters per graph, and datasets
# are generated according to these files.
# ** As for the variable u, its role is to modify the parameters so we have only one graph, and 
# several parameter files, one per value of u. Then, datasets are generated with these files.
# ** Finally, as the variable n has only an impact on the size of the generated datasets, we have one 
# graph, one file of parameters, and then one type of dataset per value of n.

# By default : u=0.2, prob=0.4, d=20, n=100
############################

############
# packages #
############
rm(list=ls())
set.seed(2)
# install.packages('/home/momal/lib/igraph_1.1.2.tar.gz',repos=NULL,lib="/home/momal/R/lib/")
# install.packages('/home/momal/lib/igraph_1.0.1.tar.gz',lib='/home/momal/R/lib/',
#                  repos = c(CRAN = "https://cran.rstudio.com"))
library(igraph, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library(huge, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library(glasso)
library(ROCR)
library(ggplot2, lib.loc="/usr/local/lib/R/site-library")
library("reshape2", lib.loc="/usr/local/lib/R/site-library")
source('/home/momal/R/TreeMixture-RML.R')
source('/home/momal/R/fonctions.R')
# source('/home/momal/Git/these/pack1/R/TreeMixture-RML.R')
# source('/home/momal/Git/these/pack1/R/fonctions.R')
mvrnorm_rml<-function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE){
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  if (EISPACK) 
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) 
    drop(X)
  else t(X)
}

###########################
# functions of generation #
###########################
generator_graph<-function(d = 20, graph = "tree", g = NULL, prob = NULL, vis = FALSE,
                          verbose = TRUE){
  gcinfo(FALSE)
  if (verbose)
    cat("Generating data from the multivariate normal distribution with the",
        graph, "graph structure....")
  if (is.null(g)) {
    g = 1
    if (graph == "hub" || graph == "cluster") {
      if (d > 40)
        g = ceiling(d/20)
      if (d <= 40)
        g = 2
    }
  }
  if (graph == "cluster") {
    if (is.null(prob)) {
      if (d/g > 30)
        prob = 0.3
      if (d/g <= 30)
        prob = min(1, 6 * g/d)
    }
    prob = sqrt(prob/2) * (prob < 0.5) + (1 - sqrt(0.5 -0.5 * prob)) * (prob >= 0.5)
  }
  g.large = d%%g
  g.small = g - g.large
  n.small = floor(d/g)
  n.large = n.small + 1
  g.list = c(rep(n.small, g.small), rep(n.large, g.large))
  g.ind = rep(c(1:g), g.list)
  rm(g.large, g.small, n.small, n.large, g.list)
  gc()
  theta = matrix(0, d, d)
  if (graph == "cluster") {
   
    for (i in 1:g) {
      tmp = which(g.ind == i)
      tmp2 = matrix(runif(length(tmp)^2, 0, 0.5), length(tmp),
                    length(tmp))
      tmp2 = tmp2 + t(tmp2)
      theta[tmp, tmp][tmp2 < prob] = 1
      rm(tmp, tmp2)
      gc()
    }
  }
  if (graph == "scale-free") {
  
    out = .C("SFGen", dd0 = as.integer(2), dd = as.integer(d),
             G = as.integer(theta), PACKAGE = "huge")
    theta = matrix(as.numeric(out$G), d, d)
  }
  if(graph=="tree"){
    theta<-SpannTree(d,prob)
  }
  if(graph=="erdos"){
    theta<- erdos(d=d,p=prob)
  }

  #browser()
  if (verbose)
    cat("done.\n")
  rm(vis, verbose)
  gc()
  return(theta = Matrix(theta, sparse = TRUE))
}
generator_param<-function(theta,v=0.3,u){
  diag(theta) = 0
  omega = theta * v
  diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
  sigma = cov2cor(solve(omega))
  omega = solve(sigma)
  sim=list(sigma=sigma,omega=omega)
  return(sim)
}
generator_data<-function(n,sigma){
  d<-ncol(sigma)
  x = scale(mvrnorm_rml(n, rep(0, d), sigma))
  sigmahat = cor(x)
  return(x)
}

############################
# functions of diagnostics #
############################

graph<-function(type,param,path){
  pdf(paste0(path,"/R/Simu/images/",type,"_",param,".pdf"),
      width=6,
      height=4,onefile=TRUE)
  print(diagnostics(paste0(path,"/R/Simu/",type,"/",param,"/spec.rds")))
  print(diagnostics(paste0(path,"/R/Simu/",type,"/",param,"/sens.rds")))
  print(diagnostics(paste0(path,"/R/Simu/",type,"/",param,"/auc.rds")))
  dev.off()
}
diagnostics<-function(file){
  tab<-data.frame(readRDS(file))
  elmts<-unlist(strsplit(file,"/"))
  variable<-elmts[length(elmts)]
  param<-elmts[length(elmts)-1]
  type<-elmts[length(elmts)-2]
  variable<-substr(variable,1,nchar(variable)-4)
  variable<-switch(variable,"auc"="AUC","sens"="Sensitivity","spec"="Specificity")
  
  ggplot(tab, aes(y=mns,x=as.numeric(as.character(tab[,5])),variable = method,group=method))+
    geom_line(aes(y=mns,color=method),size=1)+
    geom_ribbon(aes(ymin=inf, ymax=sup,fill=method), alpha=0.2)+
    labs(title = paste0(type,": effect of ",param," on ",variable),y=variable,x=param)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  
}

####################
###   Main code  ###
####################

save_params<-function(x,variable,type,path,graph,prob,u){ 
  if(variable!="u"){
    graph<-switch(variable, "d"=generator_graph(graph=type,d=x,prob=prob),
                  "prob"=generator_graph(graph=type,prob=x))
    param<-generator_param(as.matrix(graph),u=u)
  }else{
    param<-generator_param(as.matrix(graph),u=x)
  }
  file<-paste0(path,"/R/Simu/",type,"/",variable,"/param_",x,".rds")
  saveRDS(param,file)
}
compare_methods<-function(x,n,sigma,K,criterion){
  X<-generator_data(n,sigma)
  print(paste0("in sapply B = ",x))
  try({inf_treeggm<-TreeGGM(X,"FALSE")$P },silent=TRUE)
  try({inf_treeggm1step<-TreeGGM(X,"TRUE")$P},silent=TRUE)
  rho<-seq(0.0005,1.1,by=0.0002)
  tab<-tableau3D(X,rho)
  inf_glasso<-mat_rho(tab,rho,"max")
  if(exists("inf_treeggm1step") && exists("inf_treeggm")){
    diagnost<-c(diagnostic.auc.sens.spe(pred=inf_treeggm,obs=K,criterion),
                diagnostic.auc.sens.spe(pred=inf_treeggm1step,obs=K,criterion),
                diagnostic.auc.sens.spe(pred=inf_glasso,obs=K,criterion))
  }else{
    if(!exists("inf_treeggm")){
      diagnost<-c(NA,
                  diagnostic.auc.sens.spe(pred=inf_treeggm1step,obs=K,criterion),
                  diagnostic.auc.sens.spe(pred=inf_glasso,obs=K,criterion))
    }else{
      diagnost<-c(NA,NA,diagnostic.auc.sens.spe(pred=inf_glasso,obs=K,criterion))
    }
  }
  return(diagnost)
}
bootstrap_summary<-function(x,type,variable,B,path,n,criterion){
  print(paste0("seq = ",x))
  if(variable!="n"){
    param<-readRDS(paste0(path,"/R/Simu/",type,"/",variable,"/param_",x,".rds"))
    sigma<-param$sigma
    K<-param$omega
  # on génère les données B fois pour B inférences, à partir des param sauvés
    res2<-sapply(1:B,function(x) compare_methods(x,n=n,sigma=sigma,K=K,criterion=criterion))
  }else{
    param<-readRDS(paste0(path,"/R/Simu/",type,"/n/param_n.rds"))
    sigma<-param$sigma
    K<-param$omega
    res2<-sapply(1:B,function(y) compare_methods(y,n=x,sigma=sigma,K=K,criterion=criterion))
  }
  method<-c("treeggm","ggm1step","glasso")
  mns<-rowMeans(res2,na.rm = TRUE)
  sds <- apply(res2,1,function(x)sd(x,na.rm = TRUE))
  inf <- mns-qt(0.975,ncol(res2)-1)*sds/sqrt(ncol(res2))#IC95 autour de la moyenne
  sup <- mns+qt(0.975,ncol(res2)-1)*sds/sqrt(ncol(res2))
  u<-rep(x,length(mns))
  return(list(cbind(mns,inf,sup,method,u),sum(is.na(res2[,1]))/B))
}

simu<-function(type,variable,seq,n,u,prob,B,path){
  if( variable=="u") graph<-generator_graph(graph=type,prob=prob)
  if( variable=="n"){
    graph<-generator_graph(graph=type,prob=prob)
    param<-generator_param(as.matrix(graph),u=u)
    file<-paste0(path,"/R/Simu/",type,"/",variable,"/param_n.rds")
    saveRDS(param,file)
  }else{
    lapply(seq,function(x) save_params(x,type=type,variable=variable,path=path,graph=graph,prob=prob,u=u))
  }
  reussite<-data.frame(matrix(0,ncol=3,nrow=length(seq)))
  colnames(reussite)<-c("auc","sens","spec")
  for(criterion in c("auc","sens","spec")){
    res<-lapply(seq,function(x) bootstrap_summary(x,type=type,variable=variable,B=B,path=path,n=n,criterion=criterion))
    reussite[,which(colnames(reussite)==criterion)]<-unlist(lapply(res,function(x) x[[2]]))
    bla<-data.frame(do.call("rbind",lapply(res,function(x) x[[1]])))
    lignes<-which(bla[,1]=="NaN")
    if (length(lignes)!=0) bla<-bla[-lignes,]
    bla[,1:3]<-apply(bla[,1:3],2,function(x) as.numeric(x))
    saveRDS(bla,paste0(path,"/R/Simu/",type,"/",variable,"/",criterion,".rds"))
  }#sortie for
  saveRDS(reussite,paste0(path,"/R/Simu/",type,"/",variable,"/fail_rate.rds"))
}

#############
#**  RUN  **#
#############
B<-35
path<-getwd()#path =getwd() || "/home/momal/Git/these/pack1/"
parameters<-list(c(seq(7,30,2)),c(seq(20,100,10)),c(seq(0,1.5,0.2)),c(seq(0.1,0.9,0.1)))
names(parameters)<-c("d","n","u","prob")

for(type in c("scale-free","tree","erdos","cluster")){
  for(param in names(parameters)){
    simu(type,variable=param,seq=parameters[[param]],n=100,u=0.2,prob=0.4,B=B,path=path)  
    graph(type,param,path=path)
  }
}


# SMALL RUN for local tests :
# type<-"tree"
# B<-3
# simu(type,variable="d",seq(7,9,2),n=100,u=0.2,prob=0.4,B=B,path="/home/momal/Git/these/pack1/")
# param<-"n"
# graph(type,param=param,path="/home/momal/Git/these/pack1/")

# readRDS("/home/momal/Git/these/pack1/R/Simu/tree/d/taux_reussite.rds")


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

set.seed(2)
# install.packages('/home/momal/lib/igraph_1.1.2.tar.gz',repos=NULL,lib="/home/momal/R/lib/")
# install.packages('/home/momal/lib/igraph_1.0.1.tar.gz',lib='/home/momal/R/lib/',
#                  repos = c(CRAN = "https://cran.rstudio.com"))
library(igraph, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library(huge, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library(glasso)
library(ROCR)
library(dplyr)
library(tidyr)
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
  theta = matrix(0, d, d)
  if (graph == "cluster") {
   
    theta<-SimCluster(d,3,5/d,10)

  }
  if (graph == "scale-free") {
  
    out = .C("SFGen", dd0 = as.integer(2), dd = as.integer(d),
             G = as.integer(theta), PACKAGE = "huge")
    theta = matrix(as.numeric(out$G), d, d)
  }
  if(graph=="tree"){
    theta<-SpannTree(d)
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
# generator_param<-function(theta,v=0.3,u){
#   diag(theta) = 0
#   omega = theta * v
#   diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
#   sigma = cov2cor(solve(omega))
#   omega = solve(sigma)
#   sim=list(sigma=sigma,omega=omega)
#   return(sim)
# }
generator_param<-function(G,v=1){
  lambda = 1
  omega = diag(rep(lambda, ncol(G))) + G*v
  while (min(eigen(omega)$values) < 0){
    lambda = 1.1*lambda
    omega = diag(rep(lambda, ncol(G))) + G*v
  }

  sigma = cov2cor(solve(omega))
  #omega = solve(sigma)
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
  # print(diagnostics(paste0(path,"/R/Simu/",type,"/",param,"/spec.rds")))
  # print(diagnostics(paste0(path,"/R/Simu/",type,"/",param,"/sens.rds")))
  print(diagnostics(paste0(path,"/R/Simu/",type,"/",param,"/auc.rds")))
  dev.off()
}

diagnostics<-function(file){
  tab<-data.frame(readRDS(file))
  lignes<-which(is.na(tab[,1]))
  if (length(lignes)!=0) tab<-tab[-lignes,]
  tab<- gather(tab,key=method,value=value,treeggm,ggm1step,glasso)
  tab<-summarise(group_by(tab,var,method),mns=mean(value),inf=quantile(value,0.25),sup=quantile(value,0.75))
  elmts<-unlist(strsplit(file,"/"))
  variable<-elmts[length(elmts)]
  param<-elmts[length(elmts)-1]
  type<-elmts[length(elmts)-2]
  variable<-substr(variable,1,nchar(variable)-4)
  variable<-switch(variable,"auc"="AUC","sens"="Sensitivity","spec"="Specificity")
   # geom_line(aes(y=mns,color=method),size=1)+
    # geom_ribbon(aes(ymin=inf, ymax=sup,fill=method), alpha=0.2)+
  
  ggplot(tab, aes(y=mns,x=as.numeric(as.character(var)),group=method,color=method))+
    geom_errorbar(aes(ymin=inf, ymax=sup), width=.3,position=position_dodge(0.05))+
    geom_smooth(se=FALSE,size=1)+
    # geom_linerange(aes(ymin = quantile(value,0.25), ymax = quantile(value,0.75)),group=tab$method)+
    labs(title = paste0("Graph of type ",type,": effect of ",param," on ",variable,".\n Cruves of means, and 1rst and 3rd quartiles."),y=variable,x=param)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))

}

####################
###   Main code  ###
####################

save_params<-function(x,variable,type,path,graph,prob,u,nbgraph=nbgraph){ 
  if(variable!="u"){
    graph<-switch(variable, "d"=generator_graph(graph=type,d=x,prob=prob),
                  "prob"=generator_graph(graph=type,prob=x))
    
    param<-generator_param(as.matrix(graph))
  }else{
    param<-generator_param(as.matrix(graph),prob=0.1)
  }
 
  file<-paste0(path,"/R/Simu/",type,"/",variable,"/",nbgraph,"/param_",x,".rds")
  saveRDS(param,file)
  
  nb_triangles = sum(adjacent.triangles(graph_from_adjacency_matrix(as.matrix(graph),mode="max")))/3
  record(nb_triangles,x,c("nb_triangles"),paste0(path,"/R/Simu/",type,"/",variable,"/",nbgraph,"/"),rep=FALSE)
}
save_scores<-function(list_scores,save_file,param){
  names<-c("scores_treeggm_","scores_one_step_","scores_glasso_")
  for(i in 1:3){
   inf<-lapply(list_scores,function(x){x[[i]] } )
   L<-length(na.omit.list(inf))
   inf<-lapply(inf,function(x) { replace(x, is.na(x), 0)})
  saveRDS(Reduce("+", inf) / L,paste0(save_file,names[i],param,".rds"))
  }
}
na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
compare_methods<-function(x,n,sigma,K,criterion){
  X<-generator_data(n,sigma)
  
  print(paste0("in sapply B = ",x))
  temps<-rep(NA,3)
  try({T1<-Sys.time()
  inf_treeggm<-TreeGGM(X,"FALSE")$P
  T2<-Sys.time()
  temps[1]<- difftime(T2, T1) },silent=TRUE)
  try({T1<-Sys.time()
  inf_treeggm1step<-TreeGGM(X,"TRUE")$P
  T2<-Sys.time()
  temps[2]<-difftime(T2, T1) },silent=TRUE)
  T1<-Sys.time()
  # rho<-seq(0.0005,1.1,by=0.0002)
  # tab<-tableau3D(X,rho)
  # inf_glasso<-mat_rho(tab,rho,"max")
  inf_glasso<-inf_glasso_MB(X)
  diag(inf_glasso)<-0
  T2<-Sys.time()
  temps[3]<- difftime(T2, T1)
  #if(!exists("inf_treeggm1step") || !exists("inf_treeggm"))  browser()
  if(exists("inf_treeggm1step") && exists("inf_treeggm")){
    diagnost<-c(diagnostic.auc.sens.spe(pred=inf_treeggm,obs=K,criterion),
                diagnostic.auc.sens.spe(pred=inf_treeggm1step,obs=K,criterion),
                diagnostic.auc.sens.spe(pred=inf_glasso,obs=K,criterion))
    
    inferences<-list(inf_treeggm,inf_treeggm1step,inf_glasso)
  }else{
      diagnost<-c(NA, NA,diagnostic.auc.sens.spe(pred=inf_glasso,obs=K,criterion))
      inferences<-list(inf_glasso*NA,inf_glasso*NA,inf_glasso)
  }
  
  return(list(diagnost,temps,inferences,estim_reg(X,K)))
}

record <- function(var, x, col_names, path2,rep=TRUE) {
  
  if (rep){
    temps <- data.frame(var, "B" = rep(1:B,nrow(var)/B), "param" = rep(x, nrow(var)))
  }else{
    temps <- data.frame(var,  "param" = x)
  } 
  colnames(temps)[1:length(col_names)] <- col_names
 
  save_path <- paste0(path2, deparse(substitute(var)), ".rds")
  if (file.exists(save_path)) {
    tmp <- rbind(readRDS(save_path), temps)
    saveRDS(tmp, save_path)
  } else{
    saveRDS(temps, save_path)
  }
}

bootstrap_summary<-function(x,type,variable,B,path,n,criterion,nbgraph=nbgraph){
  print(paste0("seq = ",x))
  if(variable!="n"){
    param<-readRDS(paste0(path,"/R/Simu/",type,"/",variable,"/",nbgraph,"/param_",x,".rds"))
    sigma<-param$sigma
    K<-param$omega
  # on génère les données B fois pour B inférences, à partir des param sauvés
    obj<-lapply(1:B,function(x) compare_methods(x,n=n,sigma=sigma,K=K,criterion=criterion))
  }else{
    param<-readRDS(paste0(path,"/R/Simu/",type,"/n/param_n.rds"))
    sigma<-param$sigma
    K<-param$omega
    obj<-sapply(1:B,function(y) compare_methods(y,n=x,sigma=sigma,K=K,criterion=criterion))
    
  }
  res2<-do.call(rbind,lapply(obj, function(x){x[[1]]}))
  res2<-cbind(res2,rep(x,nrow(res2)))
  temps<-do.call(rbind,lapply(obj, function(x){x[[2]]}))
  scores<-do.call(list,lapply(obj, function(x){x[[3]]}))
  estim_nb_edges<-do.call(rbind,lapply(obj, function(x){x[[4]]}))
  
  
  method<-c("treeggm","ggm1step","glasso")
  path2<-paste0(path,"/R/Simu/",type,"/",variable,"/",nbgraph,"/")
  record(temps,x,method,path2)
  record(estim_nb_edges,x,c("nb_pred","nb_obs"),path2)
  save_scores(scores,path2,x)
  # mns<-colMeans(res2,na.rm = TRUE)
  # sds <- apply(res2,2,function(x)sd(x,na.rm = TRUE))
  # inf <- mns-qt(0.975,nrow(res2)-1)*sds/sqrt(nrow(res2))#IC95 autour de la moyenne
  # sup <- mns+qt(0.975,nrow(res2)-1)*sds/sqrt(nrow(res2))
  #cbind(mns,inf,sup,method,var)
  # var<-rep(x,length(mns))
  print(res2)
  #colnames(res2)<-method
  #res2$param<-rep(x,nrow(res2))
  return(res2)
}

simu<-function(type,variable,seq,n,B,prob=0.4,path,Bgraph){
  for(nbgraph in 1:Bgraph){
  print(paste0("graph number ",nbgraph))    
    if( variable=="u") graph<-generator_graph(graph=type,prob=prob)
    if( variable=="n"){
      graph<-generator_graph(graph=type,prob=0.1)
      param<-generator_param(as.matrix(graph))
      file<-paste0(path,"/R/Simu/",type,"/",variable,"/",nbgraph,"/param_n.rds")
      saveRDS(param,file)
    }else{
      lapply(seq,function(x) save_params(x,type=type,variable=variable,path=path,graph=graph,
                                         prob=(log(x)/2)/x,nbgraph=nbgraph))
    }
    for(criterion in c("auc")){
      res<-lapply(seq,function(x) bootstrap_summary(x,type=type,variable=variable,
                                                    B=B,path=path,n=n,
                                                    criterion=criterion,nbgraph=nbgraph))
      print(res)
      #auc<-data.frame(do.call("rbind",lapply(res,function(x) x[[1]])))
      auc<-data.frame(do.call("rbind",res))
      print(auc)
  
      auc[,1:3]<-apply(auc[,1:3],2,function(x) as.numeric(x))
      record(auc,nbgraph,c("treeggm","ggm1step","glasso","var"),paste0(path,"/R/Simu/",type,"/",variable,"/"))
      #reussite[,which(colnames(reussite)==criterion)]<-unlist(lapply(res,function(x) x[[2]]))
    }#sortie for
  }  
  auc<-readRDS(paste0(path,"/R/Simu/",type,"/",variable,"/auc.rds"))
    fail_rate<-auc %>%
      group_by(var) %>%
      summarise (fail=mean(is.na(treeggm)))
  saveRDS(fail_rate,paste0(path,"/R/Simu/",type,"/",variable,"/fail_rate.rds"))

}

#############
#**  RUN  **#
#############
B<-6
path<-getwd()#path =getwd() || "/home/momal/Git/these/pack1/"
parameters<-list(c(seq(7,30,2)),c(seq(20,100,10)),c(seq(0,1.5,0.2)),c(seq(0.5,5,0.5)/20))
names(parameters)<-c("d","n","u","prob")

for(type in c("tree","erdos","cluster","scale-free")){
  for(param in c("d","prob","n")){
    simu(type,variable=param,seq=parameters[[param]],n=100,B=B,path=path,Bgraph=6)  
    graph(type,param,path=path)
  }
}
# type<-"tree"
# path<-"/home/momal/Git/these/pack1/"



# # SMALL RUN for local tests :
# type<-"tree"
# B<-5
# simu(type,variable="d",c(seq(8,12,2)),n=100,B=B,path=path,Bgraph=2)
# 
# # # param<-"n"
#  graph(type,param="d",path="/home/momal/Git/these/pack1/")
#
# list_scores<-readRDS("/home/momal/Git/these/pack1/R/Simu/erdos/d/scores.rds")


# #### essai rbind auc
# tab %>%
#   group_by(u,method) %>%
#   summarise(avg = mean(mns),inf=means(mns)-)
# mns-qt(0.975,nrow(res2)-1)*sds/sqrt(nrow(res2))
# tab<-rbind(auc1,auc2)
# ggplot(tab, aes(y=mns,x=as.numeric(as.character(tab[,5])),variable = method,group=method))+
#   #geom_line(aes(y=mns,color=method),size=1)+
#   geom_smooth(aes(y=mns,color=method),size=1)+
# #  geom_ribbon(aes(ymin=inf, ymax=sup,fill=method), alpha=0.2)+
#   labs(title = paste0(type,": effect of ",param," on ",variable),y=variable,x=param)+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5))
# tab<-rbind(auc1,auc1)


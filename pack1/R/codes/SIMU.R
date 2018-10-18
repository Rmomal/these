
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

set.seed(3)
# install.packages('/home/momal/lib/igraph_1.1.2.tar.gz',repos=NULL,lib="/home/momal/R/lib/")
# install.packages('/home/momal/lib/igraph_1.0.1.tar.gz',lib='/home/momal/R/lib/',
#                  repos = c(CRAN = "https://cran.rstudio.com"))
library(igraph)
#, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library(huge)
#, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library(glasso)
library(ROCR)
library(SpiecEasi)
library(Matrix)
library(dplyr)
library(tidyr)
library(parallel)
library(PLNmodels)
library(mvtnorm)
library(ggplot2)
#, lib.loc="/usr/local/lib/R/site-library")
library(reshape2)
#, lib.loc="/usr/local/lib/R/site-library")
# source('/home/momal/R/TreeMixture-RML.R')
# source('/home/momal/R/fonctions.R')
source('/home/momal/Git/these/pack1/R/TreeMixture-RML.R')
source('/home/momal/Git/these/pack1/R/fonctions.R')
mvrnorm_rml<-function (n = 1, mu=0, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE){
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
                          verbose = TRUE,r=10){
  gcinfo(FALSE)
  if (verbose)
    #cat("Generating data from the multivariate normal distribution with the",
    #     graph, "graph structure....")
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
    #browser()
    theta<-SimCluster(d,3,prob,r) #prob=5/d ?

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
    if(sum(theta)<4){
      while(sum(theta)<4){
        theta<- erdos(d=d,p=prob)
      }
    }
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
  cste = 1
  omega = diag(rep(cste, ncol(G))) + G*v
  while (min(eigen(omega)$values) < 1e-6){
    cste = 1.1*cste
    omega = diag(rep(cste, ncol(G))) + G*v
  }
  #browser()
  sigma = cov2cor(solve(omega))
  #omega = solve(sigma)
  sim=list(sigma=sigma,omega=omega,cste=cste)
  return(sim)
}

generator_data<-function(n,sigma){
  d<-ncol(sigma)
  x = scale(mvrnorm_rml(n, rep(0, d), sigma))
  sigmahat = cor(x)
  return(x)
}
generator_composi_data<-function(Sigma,offset,covariates){#c = nb covaiables
  # vraies abondances, log normales
  n<-nrow(covariates)
  c<-ncol(covariates)# nb covariables
  p<-ncol(Sigma) # nb esp??ces
  beta<-matrix(runif(c*p),c,p)
  real_counts <- 1 + covariates %*% beta + rmvnorm(n, rep(0,nrow(Sigma)), Sigma)

  # transformation logistique en proportions
  pi<-exp(real_counts)/rowSums(exp(real_counts))
  #observed counts multinomiaux
  Y<-matrix(0,n,p)
  for(i in 1:n){
    Y[i,]<-rmultinom(1,offset[i],pi[i,]) #offset est un vecteur de taille n
  }
  return(Y)
}

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
generator_PLN_nocov<-function(Sigma,n){

  # vraies abondances, log normales
  p<-ncol(Sigma) # nb esp??ces
  # browser()
  Z<- rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
  Y = matrix(rpois(n*p, exp(Z)), n, p)
  return(list(Y,cor(Z)))
}
generator_ZpZhat<-function(Sigma,O,covariables,n=200){#O=NULL,covariables=NULL, pb avec ces argms
  Z = rmvnorm(n, sigma=Sigma);
  p<-ncol(Sigma)
  Y = matrix(rpois(n*p, exp(O+Z)), n, p)
  PLN = PLN(Y ~  offset(O)+covariables)
  ZpZ.hat = PLN$model_par$Sigma
  return(ZpZ.hat)
}
############################
# functions of diagnostics #
############################

graph<-function(type,param,path){
  pdf(paste0(path,"images/",type,"_",param,".pdf"),
      width=6,
      height=4,onefile=TRUE)
  # print(diagnostics(paste0(path,"/R/Simu/",type,"/",param,"/spec.rds")))
  # print(diagnostics(paste0(path,"/R/Simu/",type,"/",param,"/sens.rds")))
  print(diagnostics(paste0(path,type,"/",param,"/auc.rds")))
  dev.off()
}
#file<-paste0(path,type,"/",param,"/auc.rds")
diagnostics<-function(file){

  tab<-data.frame(readRDS(file))
  lignes<-which(is.na(tab[,1]))
  if (length(lignes)!=0) tab<-tab[-lignes,]
  tab<- gather(tab,key=method,value=value,treeggm,ggm1step,oracle,glasso)
  tab<-summarise(group_by(tab,var,method),mns=median(value),inf=quantile(value,0.25),sup=quantile(value,0.75))
  elmts<-unlist(strsplit(file,"/"))
  variable<-elmts[length(elmts)]
  param<-elmts[length(elmts)-1]
  type<-elmts[length(elmts)-2]
  variable<-substr(variable,1,nchar(variable)-4)
  variable<-switch(variable,"auc"="AUC","sens"="Sensitivity","spec"="Specificity")

  # geom_line(aes(y=mns,color=method),size=1)+
  # geom_ribbon(aes(ymin=inf, ymax=sup,fill=method), alpha=0.2)+
  tab$var<-as.numeric(as.character(tab$var))
  ggplot(tab, aes(y=mns,x=as.numeric(as.character(var)),group=method,color=method))+
    geom_errorbar(aes(ymin=inf, ymax=sup), width=0,position=position_dodge((max(tab$var)-min(tab$var))/100))+
    #geom_smooth(se=FALSE,size=0.3)+
    geom_point()+
    geom_line(size=0.2)+
    # geom_linerange(aes(ymin = quantile(value,0.25), ymax = quantile(value,0.75)),group=tab$method)+
    labs(y=variable,x=param)+
    scale_color_manual(values=c("#076443", "#56B4E9","#fc5f94","#E69F00" ),name="Method:",
                       breaks=c("treeggm","ggm1step","oracle", "glasso" ),
                       labels=c("EM ","1 step", "oracle","SpiecEasi" ))+
    scale_y_continuous(limits = c(0.5,1))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())
}
#title = paste0("Graph of type ",type,": effect of ",param," on ",variable,".\n Cruves of medians, and 1rst and 3rd quartiles.")
####################
###   Main code  ###
####################

save_params<-function(x,variable,type,path,graph,prob,nbgraph=nbgraph){
  if(variable!="u"){
    graph<-switch(variable, "d"=generator_graph(graph=type,d=x,prob=2/x),
                  "prob"=generator_graph(graph=type,prob=x),
                  "r"=generator_graph(graph=type,r=x))

    param<-generator_param(as.matrix(graph))
  }else{
    param<-generator_param(as.matrix(graph),prob=prob)
  }
  file<-paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",x,".rds")
  saveRDS(param,file)
  nb_triangles = sum(adjacent.triangles(graph_from_adjacency_matrix(as.matrix(graph),mode="max")))/3
  record(nb_triangles,x,c("nb_triangles"),paste0(path,type,"/",variable,"/Graphs_characteristics/Graph",nbgraph),rep=FALSE)
}
save_scores<-function(list_scores,save_file,val,nbgraph){
  names<-c("treeggm_","one_step_","oracle","glasso_")
  for(i in 1:4){
    inf<-lapply(list_scores,function(x){x[[i]] } )
    L<-length(na.omit.list(inf))
    inf<-lapply(inf,function(x) { replace(x, is.na(x), 0)})
    saveRDS(Reduce("+", inf) / L,paste0(save_file,"Graph",nbgraph,"_",names[i],val,".rds"))
  }
}
na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

compare_methods<-function(x,n,sigma,K,criterion){
  X<-generator_data(n,sigma)

  print(paste0("in sapply B = ",x))
  temps<-rep(NA,3)
  # try({
  T1<-Sys.time()

  inf_treeggm<-TreeGGM(X,"FALSE",FALSE)$P
  T2<-Sys.time()
  temps[1]<- difftime(T2, T1)
  #  },silent=TRUE)
  try({T1<-Sys.time()
  inf_treeggm1step<-TreeGGM(X,"TRUE",FALSE)$P
  T2<-Sys.time()
  temps[2]<-difftime(T2, T1) },silent=TRUE)
  T1<-Sys.time()
  # rho<-seq(0.0005,1.1,by=0.0002)
  # tab<-tableau3D(X,rho)
  # inf_glasso<-mat_rho(tab,rho,"max")
  #inf_glasso<-inf_glasso_MB(X)
  inf_glasso<-inf_spieceasi(X) # /!\ X doit ??tre comptages
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

compare_methods2<-function(x,n,Y,K,covariables,estim_cov=TRUE,corZ){
  browser()
  #offset<-matrix(offset,n,ncol(Y))

  if(estim_cov){
    PLN = PLN(Y ~ -1+covariables)
    Sigma<-PLN$model_par$Sigma
    corVEM<-cov2cor(Sigma)
  }else{
    PLN = PLN(Y ~ 1)
    Sigma<-PLN$model_par$Sigma
    corVEM<-cov2cor(Sigma)
  }
  print(paste0("in sapply B = ",x))
  temps<-rep(NA,4)
  # try({
  T1<-Sys.time()

  inf_treeggm<-TreeGGM(corVEM,"FALSE",FALSE)$P
  T2<-Sys.time()
  temps[1]<- difftime(T2, T1)
  #  },silent=TRUE)
  try({T1<-Sys.time()
  inf_treeggm1step<-TreeGGM(corVEM,"TRUE",FALSE)$P
  T2<-Sys.time()
  temps[2]<-difftime(T2, T1) },silent=TRUE)

  T1<-Sys.time()
  inf_treeggmOracle<-TreeGGM(corZ,"FALSE",FALSE)$P
  T2<-Sys.time()
  temps[3]<- difftime(T2, T1)

  # rho<-seq(0.0005,1.1,by=0.0002)
  # tab<-tableau3D(X,rho)
  # inf_glasso<-mat_rho(tab,rho,"max")
  #inf_glasso<-inf_glasso_MB(X)

  try({T1<-Sys.time()
  inf_glasso<-inf_spieceasi(Y) # /!\ Y doit ??tre comptages
  diag(inf_glasso)<-0
  T2<-Sys.time()
  temps[4]<-difftime(T2, T1) },silent=TRUE)

  #if(!exists("inf_treeggm1step") || !exists("inf_treeggm"))  browser()
  if(exists("inf_treeggm1step") && exists("inf_treeggm") && exists("inf_treeggmOracle") && exists("inf_glasso")){
    #browser()
    diagnost_total<-list(diagnostic.auc.sens.spe(pred=inf_treeggm,obs=K,method="treeggm"),
                         diagnostic.auc.sens.spe(pred=inf_treeggm1step,obs=K,method="1step"),
                         diagnostic.auc.sens.spe(pred=inf_treeggmOracle,obs=K,method="oracle"),
                         diagnostic.auc.sens.spe(pred=inf_glasso,obs=K,method="SpiecEasi"))
    diagnost<-t(do.call(rbind,lapply(diagnost_total, function(x){x[[1]]})))
    precrec<-do.call(rbind,lapply(diagnost_total, function(x){x[[2]]}))

    inferences<-list(inf_treeggm,inf_treeggm1step,inf_treeggmOracle,inf_glasso)
  }else{
    diagnost<-c(NA, NA,NA,NA)
    inferences<-list(matrix(NA,ncol(K),ncol(K)),matrix(NA,ncol(K),ncol(K)),matrix(NA,ncol(K),ncol(K)),matrix(NA,ncol(K),ncol(K)))
  }

  return(list(diagnost,temps,inferences,estim_reg(generator_PLN_nocov(Sigma,n)[[1]],K),precrec=precrec))
}

record <- function(var, x, col_names, path2, B=1, rep = TRUE) {
  if (rep) {
    frame <-
      data.frame(var,
                 "B" = rep(1:B, nrow(var) / B),
                 "param" = rep(x, nrow(var)))
  } else{
    frame <- data.frame(var,  "param" = x)
  }
  colnames(frame)[1:length(col_names)] <- col_names

  save_path <- paste0(path2, deparse(substitute(var)), ".rds")
  if (file.exists(save_path)) {
    tmp <- rbind(readRDS(save_path), frame)
    saveRDS(tmp, save_path)
  } else{
    saveRDS(frame, save_path)
  }
}

bootstrap_summary<-function(x,type,variable,B,path,n,criterion,nbgraph=nbgraph,PLN,covariables,
                            cov,estim_cov){

  print(paste0("seq = ",x))
  if(variable!="n"){
    param<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",x,".rds"))
    if( PLN){
      K<-param$omega
      #  Y<-generator_composi_data(param$sigma,offset,covariables)
      # offset<-round(matrix(runif(n*ncol(K)),n,ncol(K))*100+1)
      if(cov){

        Y<-generator_PLN(param$sigma,covariables)[[1]]
        corZ<-generator_PLN(param$sigma,covariables)[[2]]
        obj<-lapply(1:B,function(x) compare_methods2(x,n=n,Y=Y,K=K,covariables,
                                                     estim_cov=estim_cov,corZ=corZ))
      }else{
        Y<-generator_PLN_nocov(param$sigma,n)[[1]]
        corZ<-generator_PLN_nocov(param$sigma,n)[[2]]
        obj<-lapply(1:B,function(x) compare_methods2(x,n=n,Y=Y,K=K,covariables,
                                                     estim_cov=FALSE,corZ=corZ))
      }
    }else{
      K<-param$omega
      sigma<-param$sigma
      obj<-lapply(1:B,function(x) compare_methods(x,n=n,sigma=sigma,K=K,criterion=criterion))
    }

    # on g??n??re les donn??es B fois pour B inf??rences, ?? partir des param sauv??s

  }else{
    param<-readRDS(paste0(path,type,"/n/Sets_param/Graph",nbgraph,".rds"))
    if( PLN){
      K<-param$omega
      #Y<-generator_composi_data(param$sigma,offset,covariables)
      if(cov){
        covariables<-matrix(c(rep(c(1,0,1),each=round(x/3)),rep(1,max(x-3*round(x/3),0))),x,1)

        Y<-generator_PLN(param$sigma,covariables)[[1]]
        corZ<-generator_PLN(param$sigma,covariables)[[2]]
        obj<-lapply(1:B,function(y) compare_methods2(y,n=x,Y=Y,K=K,covariables,
                                                     estim_cov=estim_cov,corZ=corZ))
      }else{
        Y<-generator_PLN_nocov(param$sigma,n=x)[[1]]
        corZ<-generator_PLN_nocov(param$sigma,n=x)[[2]]
        obj<-lapply(1:B,function(y) compare_methods2(y,n=x,Y=Y,K=K,covariables,
                                                     estim_cov=FALSE,corZ=corZ))
      }
    }else{
      K<-param$omega
      sigma<-param$sigma
      obj<-lapply(1:B,function(y) compare_methods(y,n=x,sigma=sigma,K=K,criterion=criterion))
    }
  }
  #browser()
  res_auc<-do.call(rbind,lapply(obj, function(x){x[[1]]}))
  res_auc<-cbind(res_auc,rep(x,nrow(res_auc)))

  res_precrec<-do.call(rbind,lapply(obj, function(x){x[["precrec"]]}))
  res_precrec<-cbind(res_precrec,rep(x,nrow(res_precrec)))

  temps<-do.call(rbind,lapply(obj, function(x){x[[2]]}))
  scores<-do.call(list,lapply(obj, function(x){x[[3]]}))
  estim_nb_edges<-do.call(rbind,lapply(obj, function(x){x[[4]]}))

  method<-c("treeggm","ggm1step","treggmOracle","glasso")
  path2<-paste0(path,type,"/",variable,"/")
  record(temps,x,method,paste0(path2,"temps/Graph",nbgraph),B)
  record(estim_nb_edges,x,c("nb_pred","nb_obs"),paste0(path2,"Graphs_characteristics/Graph",nbgraph),B)
  save_scores(scores,paste0(path2,"Scores/"),val=x,nbgraph=nbgraph)

  #print(res_auc)
  return(list(res_auc,res_precrec))
}

simu<-function(type,variable,seq,n,B,prob=0.1,path,Bgraph,PLN=FALSE,covariables,cov=TRUE,estim_cov=TRUE,cores){
  T1<-Sys.time()
  for(nbgraph in 1:Bgraph){
    print(paste0("graph number: ",nbgraph,"// type: ", type," // var: ", variable))
    if( variable=="u") graph<-generator_graph(graph=type,prob=prob)
    if( variable=="n"){
      graph<-generator_graph(graph=type,prob=prob)
      param<-generator_param(as.matrix(graph))

      file<-paste0(paste0(path,type,"/n/Sets_param/Graph",nbgraph,".rds"))
      saveRDS(param,file)
    }else{
      lapply(seq,function(x) save_params(x,type=type,variable=variable,path=path,graph=graph,
                                         prob=prob,nbgraph=nbgraph))
    }
    for(criterion in c("auc")){

      res<-mclapply(seq,function(x) bootstrap_summary(x,type=type,variable=variable,
                                                      B=B,path=path,n=n,
                                                      criterion=criterion,nbgraph=nbgraph,PLN,covariables,
                                                      cov=cov,estim_cov=estim_cov),mc.cores=cores)
      # print(res[[1]][[1]])
      auc<-data.frame(do.call(rbind,lapply(res,function(x){x[[1]]})))
      precrec<-data.frame(do.call(rbind,lapply(res,function(x){x[[2]]})))
      #auc<-data.frame(do.call("rbind",res))
      # browser()
      auc[,1:4]<-apply(auc[,1:4],2,function(x) as.numeric(x))
      record(auc,nbgraph,c("treeggm","ggm1step","oracle","glasso","var"),paste0(path,type,"/",variable,"/"),B)
      record(precrec,nbgraph,c("prec","rec","method","var"),paste0(path,type,"/",variable,"/"),B)

    }#sortie for
  }
  auc<-readRDS(paste0(path,type,"/",variable,"/auc.rds"))
  fail_rate<-auc %>%
    group_by(var) %>%
    summarise (fail=mean(is.na(treeggm)))
  saveRDS(fail_rate,paste0(path,type,"/",variable,"/fail_rate.rds"))
  T2<-Sys.time()
  difftime(T2, T1)
}

#############
#**  RUN  **#
#############

path<-"/home/momal/Git/these/pack1/R/Simu/test/"#path =paste0(getwd(),"/R/Simu/") || "/home/momal/Git/these/pack1/R/Simu/"
parameters<-list(c(seq(10,30,2)),c(seq(10,120,10)),c(seq(0,1.5,0.2)),c(seq(0.5,5,0.5)/20),c(seq(1,50,10)))
names(parameters)<-c("d","n","u","prob","r")


cparam<-c("d","prob","n")
type=c("erdos","tree","scale-free","cluster")
n<-100
covariables<-cbind(rep(c(0,1),each=n/2),rnorm(n,8,0.5),
                   c(rep(c(1,0,1),each=round(n/3)),rep(1,n-3*round(n/3))),
                   round(runif(n)*10))
covariables2<-matrix(c(rep(c(1,0,1),each=round(n/3)),rep(1,n-3*round(n/3))),n,1)
#offset<-round(matrix(runif(n*)*10000)
# for(type in type) {
#   cparam <- ifelse(type == "tree", c("d","n"), cparam)
#   for (param in cparam) {
type<-"tree"
for(type in type) {
  for(param in c("d","n")){
    simu(
      type,
      variable = param,
      seq = parameters[[param]],
      n = n,
      B = 2,
      path = path,
      Bgraph = 100,
      PLN = TRUE,
      covariables = covariables,cov=FALSE,estim_cov=FALSE,
      cores = 3
    )
  }
}
#difftime(T2, T1)
graph(type, param, path = path)
# }
# }
##### avec covariables, stockages dans Simu/type/
for(type in c("erdos")) {
  #  cparam <- ifelse(type == "tree", c("d","n"), cparam)
  for (param in c("prob")) {

    simu(
      type,
      variable = param,
      seq = parameters[[param]],
      n = n,
      B = 2,
      path = path,
      Bgraph = 3,
      PLN = TRUE,
      covariables = covariables2,cov=TRUE,estim_cov=TRUE,
      cores = 1
    )
     graph(type, param, path = path)
  }
}






# # SMALL RUN for local tests :
# type<-"tree"
#  B<-2
# # simu(type,variable="n",c(seq(20,30,5)),n=100,B=B,path=path,Bgraph=2)
#  simu(type,variable="d",c(seq(10,14,2)),n=100,B=2,path=path,Bgraph=3,PLN=TRUE,cores=12)
# #
# # # param<-"d"
#  graph(type,param="d",path="/home/momal/Git/these/pack1")

######################
##### PLN MODEL ######
######################
# # Model parms
# lambda = 1; Omega = diag(rep(lambda, p)) + G;
# while (min(eigen(Omega)$values) < 0){lambda = 1.1*lambda; Omega = diag(rep(lambda, p)) + G}
# Sigma = solve(Omega)
# O = matrix(rnorm(n*p), n, p)
#
# # Data
# Z = rmvnorm(n, sigma=Sigma);
# Y = matrix(rpois(n*p, exp(O+Z)), n, p)
# PLN = PLN(Y ~ -1 + offset(O))
# ZpZ.hat = PLN$model_par$Sigma


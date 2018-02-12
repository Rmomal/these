rm(list=ls())
source('/home/momal/Git/these/pack1/R/TreeMixture-RML.R')
source('/home/momal/Git/these/pack1/R/fonctions.R')
library(huge)
library(LITree)
library(glasso)
library(ggplot2)
library(ROCR)
library(RColorBrewer)
library(reshape2)
library(igraph)
library(mixtools)


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
    if (is.null(u))
      u = 0.1
    if (is.null(v))
      v = 0.3
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
    if (is.null(u))
      u = 0.1
    if (is.null(v))
      v = 0.3
    out = .C("SFGen", dd0 = as.integer(2), dd = as.integer(d),
             G = as.integer(theta), PACKAGE = "huge")
    theta = matrix(as.numeric(out$G), d, d)
  }
  if(graph=="tree"){
    theta<-SpannTree(d,prob)
  }
  if(graph=="erdos"){
    theta<- graphModel$new(type =graph,size=d, p.or.m =prob)$adjmat
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
  x = scale(mvrnorm(n, rep(0, d), sigma))
  sigmahat = cor(x)
  return(x)
}

simu_u<-function(type,u,n,B){
  graph<-generator_graph(graph=type,prob=0.4)
  lapply(u,function(x){
    param<-generator_param(as.matrix(graph),u=x)
    sigma<-param$sigma
    K<-param$omega
    file<-paste0(getwd(),"/",type,"/u/param_",as.character(x),".rds")
    saveRDS(param,file)})

  for(criterion in c("auc","sens","spec")){
    res<-lapply(u,function(x){
      param<-readRDS(file)
      sigma<-param$sigma
      K<-param$omega
      res2<-sapply(1:B,function(x){
        X<-generator_data(n,param$sigma)
        inf_treeggm<-TreeGGM(X,"FALSE")$P
        inf_treeggm1step<-TreeGGM(X,"TRUE")$P
        rho<-seq(0.0005,1.1,by=0.0002)
        tab<-tableau3D(X,rho)
        inf_glasso<-mat_rho(tab,rho,"max")

        c(diagnostic.auc.sens.spe(pred=inf_treeggm,obs=K,criterion),
          diagnostic.auc.sens.spe(pred=inf_treeggm1step,obs=K,criterion),
          diagnostic.auc.sens.spe(pred=inf_glasso,obs=K,criterion))
      })
      method<-c("treeggm","ggm1step","glasso")
      ####### pb rowMeans at least two dim
      #browser()
      mns<-rowMeans(res2)
      sds <- apply(res2,1,sd)
      inf <- mns-qt(0.975,ncol(res2)-1)*sds/sqrt(ncol(res2))
      sup <- mns+qt(0.975,ncol(res2)-1)*sds/sqrt(ncol(res2))
      u<-rep(x,length(mns))
      return(cbind(mns,inf,sup,method,u))
    }
    )#sortie lapply
    bla<-data.frame(do.call("rbind",res))
    bla[,1:3]<-apply(bla[,1:3],2,function(x) as.numeric(x))
    saveRDS(bla,paste0(getwd(),"/",type,"/u/",criterion,".rds"))
  }#sortie for
}

diagnostics<-function(file){
  tab<-data.frame(readRDS(file))
  elmts<-unlist(strsplit(file,"/"))
  variable<-elmts[length(elmts)]
  param<-elmts[length(elmts)-1]
  type<-elmts[length(elmts)-2]
  variable<-substr(variable,1,nchar(variable)-4)
  variable<-switch(variable,"auc"="AUC","sens"="Sensitivity","spec"="Specificity")

ggplot(tab, aes(y=mns,x=tab[,5],variable = method,group=method))+
    geom_line(aes(y=mns,color=method),size=1)+
    geom_ribbon(aes(ymin=inf, ymax=sup,fill=method), alpha=0.2)+
    labs(title = paste0(type,": effect of ",param," on ",variable),y=variable,x=param)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))


}

##### RUN
#setwd("/home/momal/Git/these/pack1/R")
setwd("/mnt/mmip/momal/Simu/")
for(type in c("tree","erdos","cluster","scale-free")){
  simu_u(type,c(0,1),100,3)
  param<-"u"
  pdf(paste0(getwd(),"/images/",type,"_",param),
      width=6,
      height=4,onefile=TRUE)
  diagnostics(paste0(getwd(),"/",type,"/",param,"/spec.rds"))
  diagnostics(paste0(getwd(),"/",type,"/",param,"/sens.rds"))
  diagnostics(paste0(getwd(),"/",type,"/",param,"/auc.rds"))
  dev.off()
}






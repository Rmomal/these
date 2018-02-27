setwd("/home/momal/Git/these/pack1/R")
rm(list=ls()); par(pch=20, mfrow=c(2, 2), mex=3/4);
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
###############
# I. Graph generation & Simuler données
##############

p<-15 #nb dimensions

scalef<-huge.generator(d=p,graph="scale-free")$data
cluster<-huge.generator(d=p,graph="cluster")$data
erdos.graph <- graphModel$new(type = "erdos",size=p, p.or.m = 0.5) # pb type="tree"
model <- GGMmodel$new(graph=erdos.graph)
model$randomSample(n=200)
erdos<-model$getX()
#K=model$K; Sigma=model$Sigma ;save(X,K,Sigma,file="Erdos20ind5var.Rdata")


#tree
omega<-Omega(SpannTree,c(10,0.4),seq(2,5,by=0.2))
data<-simu(omega,200)
###############
# II. Inférences
##############
erdos.graph <- graphModel$new(type ="erdos",size=10, p.or.m =0.5)$adjmat

data_omega<-function(huge,type,p,n=200,prob=0.5){
  if(huge){
  data<-huge.generator(d=p,graph=type)
 # browser()
  X<-data$data
  K<-data$omega
  }else{
    if(type=="tree"){
      K<-Omega(SpannTree,c(p,prob),seq(2,5,by=0.2))
      X<-simu(K,n)
    }else{
      erdos.graph <- graphModel$new(type =type,size=p, p.or.m =prob) # pb type="tree"
  model <- GGMmodel$new(graph=erdos.graph)
  model$randomSample(n=n)
  X<-model$getX()
  K<-model$K
    }
  }
  return(list(X=X,K=K))
}

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

#tests data
data<-data_omega(FALSE,"tree",15,prob=0.3)
X<-data$X
K<-data$K
##
graph<-generator_graph(graph="erdos",prob=0.2)
param<-generator_param(as.matrix(graph),u=0)
sigma<-param$sigma
K<-param$omega
X<-generator_data(100,sigma)

#tests inference
inf_treeggm<-TreeGGM(X,"FALSE")$P
inf_treeggm1step<-TreeGGM(X,"TRUE")$P #true pour 1step
rho<-seq(0.0005,1.1,by=0.0002)
tab<-tableau3D(X,rho)
inf_glasso<-mat_rho(tab,rho,"max")


###############
# III. Comparaison
##############
### courbes ROC
# résultats variables
dev.off()
q<-fun.auc.ggplot(pred=inf_treeggm,obs=K,title="Mixture Tree GGM",rho)
g<-fun.auc.ggplot(pred=inf_glasso,obs=K,title="Glasso",rho)
h<-fun.auc.ggplot(pred=inf_treeggm1step,obs=K,title="Mixture Tree GGM 1 step",rho)
grid.arrange(q,h,g, ncol=3, clip=TRUE)


##### fonction solide
# 200 observations

d<-8:30
prob<-seq(0.25,0.8,0.05)
u<-c(0,1)

simu_u<-function(u,n,B){
  graph<-generator_graph(graph="tree",prob=0.4)
  lapply(u,function(x){
    param<-generator_param(as.matrix(graph),u=x)
    sigma<-param$sigma
    K<-param$omega
    saveRDS(param,paste0(getwd(),"/u/param_",as.character(x),".rds"))})

  for(criterion in c("auc","sens","spec")){
    res<-lapply(u,function(x){
    param<-readRDS(paste0(getwd(),"/u/param_",as.character(x),".rds"))
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
  saveRDS(bla,paste0(getwd(),"/u/",criterion,".rds"))
  }#sortie for
}

simu_u(c(0,1),100,3)

#### fonction graphique à partir de auc spec et sens
diagnostics<-function(file){
  tab<-data.frame(readRDS(file))

g<-ggplot(tab, aes(y=mns,x=u,variable = method,group=method))+
geom_line(aes(y=mns,color=method))+
geom_ribbon(aes(ymin=inf, ymax=sup,fill=method), alpha=0.2)+
  theme_bw()
return(g)
}
diagnostics(paste0(getwd(),"/u/auc.rds"))

res<-data.frame(t(res))[-1,]
colnames(res)<-c("treeggm","ggm1step","glasso")
res$u<-u[1:nrow(res)]
saveRDS(res,"u_tree.rds")


#res<-readRDS(file="nods_erdos.rds")

meltdf <- melt(res,id="u")
ggplot(meltdf,aes(x=u,y=value,colour=variable,group=variable)) +
  geom_point(size=2)+
  geom_line()+
  labs(x="Diagonal penalty",y="AUC",title="Effect of u on spanning tree graph \n with 20 nods")+
  theme(plot.title = element_text(hjust = 0.5))



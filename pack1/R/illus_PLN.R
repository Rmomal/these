rm(list=ls()); par(pch=20);
#devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
library(PLNmodels); library(sna);
library(igraph)
library(RColorBrewer)
source('/home/momal/Git/these/pack1/R/FunctionsMatVec.R')
source('/home/momal/Git/these/pack1/R/FunctionsTree.R')
source('/home/momal/Git/these/pack1/R/FunctionsInference.R')
source('/home/momal/Git/these/pack1/R/TreeMixture-RML.R')
source('/home/momal/Git/these/pack1/R/fonctions.R')
# Data
data.dir = '/home/momal/Git/these/Data/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))

# Parms
Y = as.matrix(Data$count); n = nrow(Y); p = ncol(Y)
O = Data$offset; X = Data$covariates
# dev.off()
# plot(X$distTObase^2,X$distTOground^2+X$distTOtrunk^2)
# abline(0,1,col="red")
# model_quadra<-lm(X$distTObase^2 ~ X$distTOground^2+X$distTOtrunk^2)
# model_lin<-lm(X$distTObase ~ X$distTOground+X$distTOtrunk)
# par(mfrow=c(1,2))
# qqnorm(model_lin$residuals/sd(model_lin$residuals),col="blue")
# abline(0,1)
# qqnorm(model_quadra$residuals/sd(model_quadra$residuals),col="red")
# abline(0,1)
# plot(model_lin$fitted.values,model_lin$residuals/sd(model_lin$residuals),col="blue")
# abline(0,0)
# plot(model_lin$fitted.values,model_quadra$residuals/sd(model_quadra$residuals),col="red")
# abline(0,0)
# summary(model_lin)
# summary(model_quadra)
# PLN models
PLN.offset = PLN(Y ~ 1 + offset(log(O)))
PLN.tree = PLN(Y ~ 1 + X$tree + offset(log(O)))
PLN.tree.base = PLN(Y ~ 1 + X$tree + X$distTObase + offset(log(O)))
PLN.tree.base.infect = PLN(Y ~ 1 + X$tree + X$distTObase + X$pmInfection + offset(log(O)))

# BIC
Crit = rbind(PLN.offset$criteria, PLN.tree$criteria, PLN.tree.base$criteria, PLN.tree.base.infect$criteria)
Crit
apply(Crit, 2, which.max)

# inférences
Z.offset = PLN.offset$model_par$Sigma
inf.offset<-TreeGGM(cov2cor(Z.offset),print=TRUE,step="FALSE")

Z.tree = PLN.tree$model_par$Sigma
inf.tree<-TreeGGM(cov2cor(Z.tree),print=TRUE,step="FALSE")

Z.tree.base = PLN.tree.base$model_par$Sigma
inf.tree.base<-TreeGGM(cov2cor(Z.tree.base),print=TRUE,step="FALSE")

Z.tree.base.infect = PLN.tree.base.infect$model_par$Sigma
inf.tree.base.infect<-TreeGGM(cov2cor(Z.tree.base.infect),print=TRUE,step="FALSE")

saveRDS(inf.offset,"inf_offset.rds")
saveRDS(inf.tree,"inf_tree.rds")
saveRDS(inf.tree.base,"inf_treebase.rds")
saveRDS(inf.tree.base.infect,"inf_treebaseinfect.rds")

par(mfrow=c(2,2))
hist(offset,breaks=100)
hist(tree,breaks=100)
hist(tree.base,breaks=100)
hist(tree.base.infect,breaks=100)

offset<-readRDS("/home/momal/Git/these/pack1/R/inf_offset.rds")[[1]]
tree<-readRDS("/home/momal/Git/these/pack1/R/inf_tree.rds")[[1]]
tree.base<-readRDS("/home/momal/Git/these/pack1/R/inf_treebase.rds")[[1]]
tree.base.infect<-readRDS("/home/momal/Git/these/pack1/R/inf_treebaseinfect.rds")[[1]]

plot(density(Z.tree.base.infect), main="offset vs tree densities",col="blue")
lines(density(Z.offset),col="red")
lines(density(Z.tree.base),col="darkgreen")
lines(density(Z.tree),col="purple")

# EstimM <- function(Prob){
#   p = ncol(Prob);
#   M = 2*sum(Prob[upper.tri(Prob)]>.5);
#   # hist(as.vector(Prob), breaks=p, main=paste(M, '/', sum(G[upper.tri(G)]==0)))
#   return(M)
# }
# nbNonEdge<-function(OmegaY,n,p){
#   Rpart= -diag(1/sqrt(diag(OmegaY)))%*%OmegaY%*%diag(1/sqrt(diag(OmegaY)))
# Stat = Rpart * sqrt((n-2)/(1-Rpart^2))
# Pval =  matrix(2*pt(abs(Stat), lower.tail=F, df=n-p-2), p, p)
# return(EstimM(Pval))
# }
# nbNonEdge(offset,n,p)
min(tree.base.infect[which(tree.base.infect!=0)])
tree<-tree- 0.01685373
tree.base<-tree.base- 0.01685373
tree.base.infect<-tree.base.infect- 0.01685373

dev.off()
hist(tree.base[which(tree.base>0.015)],breaks=100)
hist(tree.base)
net<-net_from_matrix(tree.base,0.05,FALSE)
degree(net)[48]

pal<-brewer.pal(8, "Spectral")
plotnet<-function(omega,nbedges){
  p<-ncol(omega)
  seuil<-sort(omega[upper.tri(omega)])[p*(p-1)/2-nbedges]
  net<-net_from_matrix(omega,seuil,FALSE)
  V(net)$label=NA
  E(net)$color=pal[7]
  E(net)$curved=.1
  V(net)$color="black"
  V(net)$size=3
  return(net)
}
par(mfrow=c(2,2))
nbedges<-200
coords <- layout_(plotnet(tree.base.infect,nbedges), nicely())
plot(plotnet(tree.base.infect,nbedges), layout = coords,main="tree.base.infect")
plot(plotnet(tree.base,nbedges),layout=coords,main="tree.base")
plot(plotnet(tree,nbedges),layout=coords,main="tree")
plot(plotnet(offset,nbedges),layout=coords,main="offset")

calcdeg<-function(matrix){
  deg<-c()
  for(nbedges in seq(6000,100,-50)){
  net<-plotnet(matrix,nbedges)
  deg<-c(deg,degree(net)[48] )
  }
  return(deg)
}

L<-list(offset,tree,tree.base,tree.base.infect)
listdeg<-lapply(L,function(x) calcdeg(x))
dev.off()
plot(listdeg[[1]],x=seq,main="degree of pathogen",ylab="",xlab="number of edges",pch=20)
points(listdeg[[2]],x=seq,col="red",pch=20)
points(listdeg[[3]],x=seq,col="blue",pch=20)
points(listdeg[[4]],x=seq,col="goldenrod",pch=20)
legend("bottomright",c("offset","tree","tree.base","tree.base.infect"),col=c("black","red","blue","goldenrod"),pch=20)

seq<-seq(6000,100,-50)
axis(1, at=seq(1,120,30),labels=seq, col.axis="black", las=2)

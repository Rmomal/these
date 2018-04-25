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
Y = as.matrix(Data$count);
O = Data$offset; X = Data$covariates

# Selection des especes
YO = Y/O
Rank = rank(colSums(YO))
plot(cumsum(sort(colSums(YO))))
Seuil = 20
Ycum = colSums(Y); Order = order(Ycum)
plot(cumsum(Ycum[Order]), col = 1+(Rank[Order]<Seuil))
Y = Y[, Rank > Seuil]; O = O[, Rank >Seuil];  n = nrow(Y); p = ncol(Y)

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
offset<-readRDS("/home/momal/Git/these/pack1/R/inf_offset.rds")[[1]]
tree<-readRDS("/home/momal/Git/these/pack1/R/inf_tree.rds")[[1]]
tree.base<-readRDS("/home/momal/Git/these/pack1/R/inf_treebase.rds")[[1]]
tree.base.infect<-readRDS("/home/momal/Git/these/pack1/R/inf_treebaseinfect.rds")[[1]]

par(mfrow=c(2,2))
hist(offset,breaks=100)
hist(tree,breaks=100)
hist(tree.base,breaks=100)
hist(tree.base.infect,breaks=100)

dev.off()
plot(density(tree.base.infect), main="offset vs tree densities",col="blue",xlim=c(0,0.1))
lines(density(offset),col="red")
lines(density(tree),col="darkgreen")
lines(density(tree.base),col="purple")

EstimM <- function(Prob){
  p = ncol(Prob);
  M = 2*sum(Prob[upper.tri(Prob)]>.5);
  # hist(as.vector(Prob), breaks=p, main=paste(M, '/', sum(G[upper.tri(G)]==0)))
  return(M)
}
nbNonEdge<-function(OmegaY,n,p){
  Rpart= -diag(1/sqrt(diag(OmegaY)))%*%OmegaY%*%diag(1/sqrt(diag(OmegaY)))
Stat = Rpart * sqrt((n-2)/(1-Rpart^2))
Pval =  matrix(2*pt(abs(Stat), lower.tail=F, df=n-p-2), p, p)
hist(Pval)
return(EstimM(Pval))
}
noffset<-nbNonEdge(Z.offset,n,p)
ntree<-nbNonEdge(Z.tree,n,p)
ntreebase<-nbNonEdge(Z.tree.base,n,p)
ntreebaseinfect<-nbNonEdge(Z.tree.base.infect,n,p)


pal<-brewer.pal(8, "Spectral")
plotnet<-function(omega,nbedges){
  p<-ncol(omega)
  seuil<-sort(omega[upper.tri(omega)])[p*(p-1)/2-nbedges/2]
  net<-net_from_matrix(omega,seuil,FALSE)
  V(net)$label=NA
  E(net)$color=pal[7]
  E(net)$curved=.1
  V(net)$color="black"
  V(net)$size=3
  return(net)
}
par(mfrow=c(2,2))

coords <- layout_(plotnet(tree.base.infect,ntreebaseinfect), nicely())
plot(plotnet(tree.base.infect,ntreebaseinfect), layout = coords,main="tree.base.infect")
plot(plotnet(tree.base,ntreebase),layout=coords,main="tree.base")
plot(plotnet(tree,ntree),layout=coords,main="tree")
plot(plotnet(offset,noffset),layout=coords,main="offset")

net<-plotnet(offset,noffset)
degree(net)[48]
gsize(net)
net<-plotnet(tree,ntree)
degree(net)[48]
gsize(net)
net<-plotnet(tree.base,ntreebase)
degree(net)[48]
gsize(net)
net<-plotnet(tree.base.infect,ntreebaseinfect)
degree(net)[48]
gsize(net)

# calcdeg<-function(matrix){
#   deg<-c()
#   for(nbedges in seq(6000,100,-50)){
#   net<-plotnet(matrix,nbedges)
#   deg<-c(deg,degree(net)[48] )
#   }
#   return(deg)
# }
#
# L<-list(offset,tree,tree.base,tree.base.infect)
# listdeg<-lapply(L,function(x) calcdeg(x))
# dev.off()
# plot(listdeg[[1]],x=seq,main="degree of pathogen",ylab="",xlab="number of edges",pch=20)
# points(listdeg[[2]],x=seq,col="red",pch=20)
# points(listdeg[[3]],x=seq,col="blue",pch=20)
# points(listdeg[[4]],x=seq,col="goldenrod",pch=20)
# legend("bottomright",c("offset","tree","tree.base","tree.base.infect"),col=c("black","red","blue","goldenrod"),pch=20)

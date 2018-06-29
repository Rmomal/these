rm(list=ls()); par(pch=20);
#devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
library(PLNmodels); library(sna);
library(igraph)
library(RColorBrewer)
library(ggplot2)
source('/home/momal/Git/these/pack1/R/FunctionsMatVec.R')
source('/home/momal/Git/these/pack1/R/FunctionsTree.R')
source('/home/momal/Git/these/pack1/R/FunctionsInference.R')
source('/home/momal/Git/these/pack1/R/TreeMixture-RML.R')
source('/home/momal/Git/these/pack1/R/fonctions.R')
# Data
data.dir = '/home/momal/Git/these/Data/Oaks-CVacher/'
data.dir = '../../Data/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))
theme_set(theme_bw())
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

Y = Y[X$tree=="2", Rank > Seuil];
O = O[X$tree=="2", Rank >Seuil];
n = nrow(Y); p = ncol(Y)

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

######
#Spieceasi
#####
inf_spiec<-inf_spieceasi(Y)
# PLN models
PLN.offset = PLN(Y ~ 1 + offset(log(O)))
PLN.dist = PLN(Y ~ 1 + X$distTObase[X$tree=="2"]+ X$distTOtrunk[X$tree=="2"] + X$distTOground[X$tree=="2"]  + offset(log(O)))
PLN.orient = PLN(Y ~ 1 + X$orientation[X$tree=="2"] + offset(log(O)))
PLN.dist.orient = PLN(Y ~ 1 + X$distTObase[X$tree=="2"] + X$orientation[X$tree=="2"] + offset(log(O)))
# BIC
Crit = rbind(PLN.offset$criteria, PLN.dist$criteria, PLN.orient$criteria,PLN.dist.orient$criteria)
apply(Crit, 2, which.max)

# inférences
Z.offset = PLN.offset$model_par$Sigma
inf.offset<-TreeGGM(cov2cor(Z.offset),print=TRUE,step="FALSE")
Z.dist = PLN.dist$model_par$Sigma
inf.dist<-TreeGGM(cov2cor(Z.dist),print=TRUE,step="FALSE")
Z.orient = PLN.orient$model_par$Sigma
inf.orient<-TreeGGM(cov2cor(Z.orient),print=TRUE,step="FALSE")
Z.dist.orient = PLN.dist.orient$model_par$Sigma
inf.dist.orient<-TreeGGM(cov2cor(Z.dist.orient),print=TRUE,step="FALSE")

heatmap(solve(Z.offset),Rowv = NA,Colv = NA)
heatmap(inf.offset$P,Rowv = NA,Colv = NA)
heatmap(solve(Z.dist),Rowv = NA,Colv = NA)
heatmap(inf.dist$P,Rowv = NA,Colv = NA)

saveRDS(inf.offset,"inf_offset.rds")
saveRDS(Z.offset,"Z.offset.rds")
saveRDS(inf.dist,"inf_dist.rds")
saveRDS(Z.dist,"Z.dist.rds")
saveRDS(inf.orient,"inf_orient.rds")
saveRDS(Z.orient,"Z.orient.rds")
saveRDS(inf_spiec,"spiec.rds")
saveRDS(Z.dist.orient,"Z.dist.orient.rds")
saveRDS(inf.dist.orient,"inf_dist_orient.rds")
########
offset<-readRDS("/home/momal/Git/these/pack1/R/inf_offset.rds")[[1]]
dist<-readRDS("/home/momal/Git/these/pack1/R/inf_dist.rds")[[1]]
orient<-readRDS("/home/momal/Git/these/pack1/R/inf_orient.rds")[[1]]
Z.orient<-readRDS("/home/momal/Git/these/pack1/R/Z.orient.rds")
Z.offset<-readRDS("/home/momal/Git/these/pack1/R/Z.offset.rds")
Z.dist<-readRDS("/home/momal/Git/these/pack1/R/Z.dist.rds")
spiec<-readRDS("/home/momal/Git/these/pack1/R/spiec.rds")
Z.dist.orient<-readRDS("/home/momal/Git/these/pack1/R/Z.dist.orient.rds")
dist_orient<-readRDS("/home/momal/Git/these/pack1/R/inf_dist_orient.rds")[[1]]

# par(mfrow=c(2,2))
# hist(offset,breaks=100)
# hist(tree,breaks=100)
# hist(tree.base,breaks=100)
# hist(tree.base.infect,breaks=100)
#
# dev.off()
# plot(density(tree.base.infect), main="offset vs tree densities",col="blue",xlim=c(0,0.1))
# lines(density(offset),col="red")
# lines(density(tree),col="darkgreen")
# lines(density(tree.base),col="purple")

#@ nb non edges from pvalue matrix
# EstimM <- function(Prob){
#   p = ncol(Prob);
#   M = 2*sum(Prob[upper.tri(Prob)]>.5);
#   return(M)
# }
# #@ build pvalue matrix and return nb non edges using estimM
# nbNonEdge<-function(OmegaY,n,p){
#   Rpart= -diag(1/sqrt(diag(OmegaY)))%*%OmegaY%*%diag(1/sqrt(diag(OmegaY)))
#   Stat = Rpart * sqrt((n-2)/(1-Rpart^2))
#   Pval =  matrix(2*pt(abs(Stat), lower.tail=F, df=n-p-2), p, p)
#   hist(Pval)
#   return(EstimM(Pval))
# }
# par(mfrow=c(2,2))
# n<-nrow(Y)
# p<-ncol(Y)
# noffset<-nbNonEdge(Z.offset,n,p)
# ntree<-nbNonEdge(Z.tree,n,p)
# ntreebase<-nbNonEdge(Z.tree.base,n,p)
# ntreebaseinfect<-nbNonEdge(Z.tree.base.infect,n,p)
#
# #@ build net from precision matrix with specified number of non edges by default,
# #@ edges if boolean FALSE
# net_precision_density<-function(omega,nbnonedges, absence=TRUE){
#   p<-ncol(omega)
#   nb<-ifelse(absence,nbnonedges, p*(p-1)/2)
#   seuil<-sort(omega[upper.tri(omega)])[nb]
#   net<-net_from_matrix(omega,seuil,FALSE)
#   V(net)$label=NA;  E(net)$color="darkolivegreen3";  E(net)$curved=.1;  V(net)$color="black";  V(net)$size=3
#   return(net)
# }
#
# net1<-net_precision_density(offset,noffset)
# net2<-net_precision_density(tree,ntree)
# net3<-net_precision_density(tree.base,ntreebase)
# net4<-net_precision_density(tree.base.infect,ntreebaseinfect)
# c(degree(net1)[44],degree(net2)[44],degree(net3)[44],degree(net4)[44])
# c(gsize(net1),gsize(net2),gsize(net3),gsize(net4))

net_seuil<-function(omega,seuil){
  net<-net_from_matrix(omega,seuil,FALSE)
  V(net)$label = NA
  V(net)[44]$label = colnames(omega)[44]
  E(net)$color = "black"
  E(net)$curved = .1
  E(net)$width=3
  deg <- degree(net, mode="out")
  V(net)$size <- deg*2+2
  V(net)$color = "darkolivegreen3"
  return(net)
}
net_nbedges<-function(omega,size){
  pal<-brewer.pal(10, "Spectral")
 # seuil<-sort(omega[upper.tri(omega)],decreasing=TRUE)[nbedges]
  net<-net_from_weight(plotG(omega),seuil,FALSE)
  E(net)$width=E(net)$weight*size
  V(net)$label = NA
 V(net)[44]$label = "EA"
 V(net)[14]$label = "F19"
  V(net)$size =1
  vcol <- rep("grey40", vcount(net))
  vcol[V(net)[44]] <- "gold"
  vcol[V(net)[14]] <- "deepskyblue"
  ecol <- rep("gray80", ecount(net))
  inc.edges44 <- incident(net,  V(net)[44], mode="all")
  inc.edges14 <- incident(net,  V(net)[14], mode="all")
   ecol[inc.edges44] <- "orange"
   ecol[inc.edges14] <- "#0086CC"

 # E(net)$curved = .1
  #E(net)$width=3
  # deg <- degree(net, mode="out")
  # V(net)$size <- deg*2+1
 V(net)$color = "darkolivegreen3"
  # V(net)$color=pal

  # ecol <- rep("gray80", ecount(net))
  # ecol[inc.edges] <- "orange"
  # ew <- rep(2, ecount(net))
  # ew[inc.edges] <- 4
  # vcol <- rep("grey40", vcount(net))
  # vcol[V(net)[44]] <- "gold"
  # neigh.nodes <- neighbors(net, V(net)[44], mode="out")
  # # Set colors to plot the neighbors:
  # vcol[neigh.nodes] <- "#ff9d00"
  V(net)$color=vcol
  E(net)$color=ecol
  # E(net)$width=ew
  return(net)
}

par(mfrow=c(2,2))


colnames(offset)<-colnames(Y)
colnames(dist)<-colnames(Y)
colnames(orient)<-colnames(Y)
colnames(spiec)<-colnames(Y)
colnames(dist_orient)<-colnames(Y)


hist(offset,breaks=100)
hist(dist,breaks=100)
hist(orient,breaks=100)
hist(dist_orient,breaks=100)
summary(c(offset))
summary(c(dist))
summary(c(orient))
summary(c(dist_orient))
#seuil<-0.11

coords1 <- layout_(net_spiec, nicely())
coords<-layout_(net1, as_star(center = V(net1)[44]))
par(mfrow=c(1,1))
plot(net1,layout=coords,main="")
plot(net2, layout = coords,main="")
plot(net3, layout = coords,main="orient")
plot(net_spiec, layout = coords,main="Spieceasi")
plot(net4, layout = coords,main="DistOrient")

############# NETS ##############
size<-1.5
net1<-net_nbedges(offset,size)
net2<-net_nbedges(dist,size)
net3<-net_nbedges(orient,size)
spiec2<-spiec
spiec2[which(spiec2==0)]<-1e-4
net_spiec<-net_nbedges(spiec2,size)
net4<-net_nbedges(dist_orient,size)
############ PLOTS ############

par(mfrow=c(1,3))
coords<-coords
plot(net_spiec, layout = coords,main="Spieceasi")
plot(net1,layout=coords,main="")
plot(net2, layout = coords,main="dist")
#plot(net3, layout = coords,main="orient")
plot(net4, layout = coords,main="DistOrient")
########### BOXES #############
data<-data.frame(degree=c(degree(net1),degree(net2),degree(net3),degree(net_spiec),degree(net4)),
                model=rep(c("Offset","Distance","Orientation","Spieceasi","DistOrient"),each=length(degree(net1))))
p<-ggplot(data, aes(x=model, y=degree,fill=model)) +
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
  theme(legend.position="none")
par(mfrow=c(1,1))
p+geom_boxplot(notch=TRUE,width=0.5)
p+geom_violin(alpha = 0.5)  + stat_summary(fun.y=median, geom="point", shape=16, size=2)
###############################
p<-ncol(Y)
color<-matrix(1,p,p)
color[44,]<-2

color[,44]<-2
color[14,]<-3

color[,14]<-3
plot(offset,dist_orient,col=color)




# M1<-degree(net1)
# M4<-degree(net4)
# degree_table<-data.frame(nods=1:ncol(Y),M1=M1,M4=M4)
# attach(degree_table)
# dev.off()
# plot(nods,M1-M4)
# abline(v=48,h=c(-20,0,20))
# plot(density(M1-M4))
# hist(M1-M4,breaks=50)
# summary(M1-M4)
# degree_table2<-degree_table[order(M1),]
# attach(degree_table2)
# plot(M1,M1-M4)

###################
calcdeg<-function(matrix){
  deg<-c()
  for(nbedges in seq(4000,100,-50)){
  net<-net_precision_density(matrix,nbedges,FALSE)
  deg<-c(deg,centralization.degree(net)$res[44] )
  }
  return(deg)
}
seq<-seq(4000,100,-50)
matrix<-offset
L<-list(offset,tree,tree.base,tree.base.infect)
listdeg<-lapply(L,function(x) calcdeg(x))
dev.off()
plot(listdeg[[1]],x=seq,main="degree of pathogen",ylab="",xlab="number of edges",pch=20)
points(listdeg[[2]],x=seq,col="red",pch=20)
points(listdeg[[3]],x=seq,col="blue",pch=20)
points(listdeg[[4]],x=seq,col="goldenrod",pch=20)
legend("bottomright",c("offset","tree","tree.base","tree.base.infect"),col=c("black","red","blue","goldenrod"),pch=20)

####################
# comparaison des graphs à densité égale
# inf_glasso_MBbis<-function(X){
#   S<-X
#   d<-ncol(X)
#   log.lambda.min <- -5
#   log.lambda.max <- log(get.lambda.l1(S))
#   log.lambda <- seq(log.lambda.min, log.lambda.max, length.out = d*2)
#   MB.res <- lapply(exp(log.lambda), function(lambda) glasso(S, lambda, trace = FALSE, approx = FALSE, penalize.diagonal = FALSE))
#   adjmat.array <- simplify2array(Map("*",exp(log.lambda),lapply(MB.res, function(x){ (abs(x$wi)>0)*1})))
#   # Let us replace each edge by the  largest Glasso lambda where it disappears (or a sum related to this)
#   K.score <- apply(adjmat.array,c(1,2),sum)
#   K.score <- K.score / max(K.score)
#   return(K.score)
# }


in_offset_glass<-inf_glasso_MBbis(Z.offset)
in_treebaseinfect_glass<-inf_glasso_MBbis(Z.tree.base.infect)
saveRDS(in_offset_glass,"inf_glass_offset.rds")
saveRDS(in_treebaseinfect_glass,"inf_glass_treebaseinfect.rds")

omega<-offset
matrice_adj_from_density<-function(omega,nbedges){
  p<-ncol(omega)
  seuil<-sort(omega[upper.tri(omega)])[p*(p-1)/2-nbedges]
  M<-matrice_adj(omega,seuil)
  return(M)
}
correspondances<-tibble(nbedges=seq(50,4350,50),gl00mt=NA,gl10mt=NA,gl01mt=NA,gl11mt=NA)
for(i in 1:nrow(correspondances)){
  nbedges<-correspondances$nbedges[i]
  offset_glasso<-matrice_adj_from_density(in_offset_glass,nbedges)
  offset_MixTreeGGM<-matrice_adj_from_density(offset,nbedges)
  if(nbedges<3750){
     correspondances[i,2:5]<-c(table(offset_glasso[upper.tri(offset_glasso)],
                                  offset_MixTreeGGM[upper.tri(offset_MixTreeGGM)]))
  }else{
    tab<-table(offset_glasso[upper.tri(offset_glasso)],
               offset_MixTreeGGM[upper.tri(offset_MixTreeGGM)])
    correspondances[i,2:5]<-c(NA,tab[1],NA,tab[2])
  }

}

long<-correspondances %>%
  gather(decision,count,-nbedges)

p1<-ggplot(long,aes(nbedges,count,colour=factor(decision)))+
  geom_point()+
  labs(title="modele offset")
############################
############################

correspondances2<-tibble(nbedges=seq(50,4350,50),gl00mt=NA,gl10mt=NA,gl01mt=NA,gl11mt=NA)
for(i in 1:nrow(correspondances2)){
  nbedges<-correspondances2$nbedges[i]
  complet_glasso<-matrice_adj_from_density(in_treebaseinfect_glass,nbedges)
  complet_MixTreeGGM<-matrice_adj_from_density(tree.base.infect,nbedges)
  if(nbedges<3650){
    correspondances2[i,2:5]<-c(table(complet_glasso[upper.tri(complet_MixTreeGGM)],
                                     complet_MixTreeGGM[upper.tri(complet_MixTreeGGM)]))
  }else{
    tab<-table(complet_glasso[upper.tri(complet_glasso)],
               complet_MixTreeGGM[upper.tri(complet_MixTreeGGM)])
    correspondances2[i,2:5]<-c(NA,tab[1],NA,tab[2])
  }

}

long2<-correspondances2 %>%
  gather(decision,count,-nbedges)

p2<-ggplot(long2,aes(nbedges,count,colour=factor(decision)))+
  geom_point()+
  labs(title="modele complet")


devtools::install_github("thomasp85/patchwork")
library(patchwork)


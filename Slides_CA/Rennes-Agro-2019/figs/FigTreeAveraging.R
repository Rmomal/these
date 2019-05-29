# Plots for tree averaging principle

rm(list=ls())
library(sna); library(vegan)
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('/home/robin/PgmR/Network/FunctionsSimul.R')
source('/home/robin/PgmR/Network/FunctionsTree.R')
FigDir = ''

# Parms
p = 4; P = p*(p-1)/2
beta.parm = 1/3
edge.width = 20; tree.edge.width = sqrt(p^(p-2)); node.width = 10; label.width = 5; node.col = 1; edge.col = 8
T.nb = 5
seed = 2; set.seed(seed) # 1
par(mfrow=c(1, 1))
export = T

# Functions
F_PlotWeightedGraph <- function(W, xy, edge.prob=F){
   p = nrow(W)
   xlim = 1.1*c(min(xy[, 1]), max(xy[, 1]))
   ylim = 1.1*c(min(xy[, 2]), max(xy[, 2]))
   plot(xy[, 1], xy[, 2], type='p', col=0, pch=20, cex=0, axes=0, 
        xlab='', ylab='', main='', xlim=xlim, ylim=ylim)
   invisible(sapply(1:(p-1), function(i){sapply((i+1):p, function(j){
      if(W[i, j] > 0){
         cat(i, j, '')
         lines(xy[c(i, j), 1], xy[c(i, j), 2], lwd=edge.width*W[i, j], col=edge.col)
         # plot(xy[c(i, j), 1], xy[c(i, j), 2], lwd=edge.width**W[i, j], col=edge.col, type='l')
         lines(xy[c(i, j), 1], xy[c(i, j), 2], lwd=edge.width*W[i, j], col=edge.col)
         if(edge.prob){text(mean(xy[c(i, j), 1]), mean(xy[c(i, j), 2]), labels=round(100*W[i, j], 2))}
      }
   })}))
   points(xy[, 1], xy[, 2], type='p', col=node.col, pch=1, cex=node.width)
   text(xy[, 1], xy[, 2], labels = (1:p), cex=label.width)
}

# Node position
theta = -2*acos(-1)*(1.5:(p+.5))/p; xy = cbind(cos(theta), sin(theta))

# Edge proba
beta = F_Vec2Sym(rbeta(P, beta.parm, beta.parm))
L = diag(rowSums(beta)) - beta; B = det(L[-1, -1])
Pedge = Kirshner(beta)$P
F_PlotWeightedGraph(Pedge, xy)

# Some tree
T.list = list(); T.prob = rep(0, T.nb)
for (t in 1:T.nb){
   cat(t, "")
   T.list[[t]] = SimTree(p)
   T.prob[t] = sqrt(prod(beta[which(T.list[[t]]==1)])) / B
   if(export){pdf(paste0(FigDir, 'FigTreeAveraging-p', p, '-tree', t, '-seed', seed, '.pdf')); par(mex=0.1)}
   # F_PlotWeightedGraph(tree.edge.width*T.prob[t]*T.list[[t]], xy)
   F_PlotWeightedGraph(.5*T.list[[t]], xy)
   if(export){dev.off()}
}

if(export){pdf(paste0(FigDir, 'FigTreeAveraging-p', p, '-avgtree-seed', seed, '.pdf')); par(mex=0.1)}
F_PlotWeightedGraph(Pedge, xy)
if(export){dev.off()}

if(export){pdf(paste0(FigDir, 'FigTreeAveraging-p', p, '-graph-seed', seed, '.pdf')); par(mex=0.1)}
F_PlotWeightedGraph(.5*(Pedge>.5), xy)
if(export){dev.off()}

round(100*T.prob, 2)

p^(p-2)


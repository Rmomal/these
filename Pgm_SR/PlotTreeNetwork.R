# Trace des reseaux

rm(list=ls()); 
library(sna); library(mclust); 

# Res
res.dir = '../pack1/R/'
res.offset = readRDS(paste0(res.dir, 'inf_offset.rds')); P.offset = res.offset$P; 
res.orient = readRDS(paste0(res.dir, 'inf_orient.rds')); P.orient = res.orient$P; 

# Function
P = P.offset; P = P.orient
# P = P[1:20, 1:20]
plotG <- function(P){
   p = ncol(P); 
   # Coords
   X = cbind(cos(2*acos(-1)*(1:p)/p), sin(2*acos(-1)*(1:p)/p))
   
   # Colors
   edge.col = matrix(8, p, p)
   edge.col[, 44] = edge.col[44,] = 2
   edge.col[, 14] = edge.col[14,] = 3
   node.col = rep(1, p)
   node.col[44] = 2
   node.col[14] = 3
   
   # Weights
   W = log((p-1)*P)
   GM = Mclust(as.vector(W[upper.tri(W)]), G=2, modelNames='E')
   tau = GM$z[, which.max(GM$parameters$mean)]
   W = matrix(0, p, p); W[upper.tri(W)] = tau; W = W + t(W)

   # Plot
   #hist(W, breaks=p, main='')
   # gplot(W, gmode='graph', edge.lwd=W*10,edge.width=W*10, edge.col=edge.col,
   #       vertex.col=node.col, coord=X)
   return(W)
}

# Plot
par(mfrow=c(2, 2), mex=.3)
plotG(offset)
plotG(P.orient)

plotG(spiec2)

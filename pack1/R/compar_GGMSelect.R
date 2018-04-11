library("glasso")
library(GGMselect)
library(network)
p = 30
n = 100
# ----------------------------------------
# Random graph generator: use of simulateGraph
# ----------------------------------------
eta = 0.11
Gr <- simulateGraph(p, eta)
X <- rmvnorm(n, mean = rep(0, p), sigma = Gr$C)
# ----------------------------------------
# Graph selection with family C01:  use of selectFast
# ----------------------------------------
GRest <- selectFast(X, family = "C01")
omega<-TreeGGM(X,"FALSE")$P
hist(omega[upper.tri(omega)])estimDensity(X,0.1,TRUE)
Ginf<-matrice_adj(omega,0.1)
# ----------------------------------------
# Plot the result with the help of the package network
# ----------------------------------------


MyFamily <- NULL
gV <- network(Gr$G)
g <- network(GRest$C01$G)
glas <- network(GMF$G)
inf<-network(Ginf)
par(mfrow=c(2,2), pty = "s")
 a <- plot(gV, usearrows = FALSE)
 title(sub="Simulated graph")
 plot(g, coord=a, usearrows = FALSE)
title(sub="Graph selected with C01 family")
plot(inf,coord=a,usearrows=FALSE)
title(sub="Graph from mixture trees")
plot(glas,coord=a,usearrows=FALSE)
title(sub="select a graph within MyFamily")
for (j in 1:3){
   MyFamily[[j]] <- abs(sign(glasso(cov(X),rho=j/5)$wi))
   diag(MyFamily[[j]]) <- 0
 }
 # select a graph within MyFamily
 GMF <- selectMyFam(X,MyFamily)
 
 
 ### AUC
 diagnostic.auc.sens.spe(Ginf,solve(Gr$C),"auc")
 diagnostic.auc.sens.spe(GRest$C01$G,solve(Gr$C),"auc")
 diagnostic.auc.sens.spe(GMF$G,solve(Gr$C),"auc")
 
 
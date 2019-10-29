# essais cluster
library(tidygraph)
library(ggraph)

SimCluster2<-function(p, k, dens, r){
  beta = dens / (r / k + (k - 1) / k)
  alpha = r * beta
  while (alpha > 1) {
    r = .9 * r
    beta = dens / (r / k + (k - 1) / k)
    alpha = r * beta
  }
  print(paste0("alpha =",alpha))
  print(paste0("beta =",beta))
  print(paste0("dens =",dens))
  print(paste0("r =",r))
  
  Z = t(rmultinom(p, 1, rep(1 / k, k)))
 groupe=Z%*%1:3
  Z = Z %*% t(Z)
  diag(Z) = 0
  ZZ = F_Sym2Vec(Z)
  G = F_Vec2Sym(rbinom(p * (p - 1) / 2, 1, alpha * ZZ + beta * (1 - ZZ)))
  return(list(G=G, groupes=groupe))
}
p=30
k=3
dens=5/p
r=20

gr<-SimCluster2(p,k,dens,r)
draw_network(gr$G,pal="black", layout="fr", groupes=gr$groupes, curv=0.1)


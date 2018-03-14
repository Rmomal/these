#### forme des graphs
library(RColorBrewer)
library(igraph)
pal<-brewer.pal(8, "Spectral")
omega<-param_11$omega
net<-net_from_matrix(omega,1e-16,FALSE)
clp <- cluster_optimal(net)
V(net)$community <- clp$membership
E(net)$color=pal[7]
V(net)$label = NA
E(net)$curved=.1
V(net)$frame.color=pal[1]
plot(net, vertex.color=pal[V(net)$community])

clust<-function(methode){
  if(methode=="cluster_spinglass"){
    net<-delete_vertices(net,deg<1)
  }else{
    net<-net
  }
  clp <- get(methode)(net)
  plot(clp, net, main="tree graph 11 nods, edge prob=0.4")
}

clust("cluster_optimal")

###### temps d'inférences

temps<-readRDS("/home/momal/Git/these/pack1/R/Simu/erdos/d/temps.rds")

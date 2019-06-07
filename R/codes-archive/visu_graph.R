#### forme des graphs
library(RColorBrewer)
library(igraph)
graphics.off()
pal<-brewer.pal(8, "Spectral")
type<-"cluster"
variable<-"dens"
nbgraph<-2
valeur=0.25
omega<-readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLN.2.0/",type,"/",
                      variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds"))$omega
net<-net_from_matrix(omega,1e-16,FALSE)
#net<-net_from_matrix(G,1e-16,FALSE)
# clp <- cluster_optimal(net)
# V(net)$community <- clp$membership
E(net)$color="black"
#V(net)$label = colnames(offset)
E(net)$curved=.1
V(net)$color="black"
V(net)$size=3
pl<-plot(net)
plot(net,coord=pl)

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

###### temps d'inf??rences

temps<-readRDS("/home/momal/Git/these/pack1/R/Simu/erdos/d/temps.rds")



#############
# ggraph
set_graph_style()
type<-"erdos"
variable<-"d"
nbgraph<-2
valeur=20

build_net<-function(type, variable, nbgraph, valeur){
  omega<-readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLN.2.0/",type,"/",
                        variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds"))$omega
  net_from_matrix(omega,1e-16,FALSE)%>% 
  as_tbl_graph() %>% 
  ggraph(layout="kk")+
  geom_edge_link()+
  geom_node_point(size=3, color="blue")
}
build_net(type, variable, nbgraph=45,valeur)

#connexes : 7,23, 27,28










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


build_net<-function(type, variable, nbgraph, valeur){
  if(type=="scale-free"){
    omega<-readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLN.2.0/",type,"/",
                          variable,"/Sets_param/Graph",nbgraph,".rds"))$omega
  }else{
 
    omega<-readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLN_nonfav/",type,"/",
                          variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds"))$omega
  }

    net_from_matrix(omega,1e-16,FALSE)%>% 
  as_tbl_graph() %>% 
  ggraph(layout="kk")+
  geom_edge_link()+
  geom_node_point(size=3, color="blue")+labs(title=type)
}
type<-"cluster"
variable<-"r"
nbgraph<-2
valeur=16

g1<-build_net("erdos", "prob", nbgraph=45,0.225)
g2<-build_net("cluster", "r", nbgraph=45,16)
g3<-build_net("scale-free", "n", nbgraph=45)

#connexes : 7,23, 27,28

grid.arrange(g1,g2,g3,nrow=1,ncol=3)







######################
# I. idée pour calculer la covariance partielle (et avoir son signe) :
# une fois un réseau obtenu,
# 1.identifier les voisins directs de deux noeuds dont on veux connaître la covariance partielle
# 2. réduire la matrice Sigma à la matrice S11-2 de dimension 2+#(voisins directs de la paire)
# 3. inverser le bloc correspondant aux voisins de la paire
# 4. calculer la covariance partielle

setwd("/Users/raphaellemomal/these/R/codes")

# Réseau jouet:

Graph1_30 <- readRDS("/Users/raphaellemomal/simulations/Simu/PLN.2.0/erdos/d/Sets_param/Graph1_30.rds")
graph=as_tbl_graph(Graph1_30$omega)

g=graph%>%
  ggraph()+
  geom_edge_link()+
  geom_node_point()+
  geom_node_text(label = 1:30, color="red", size=5)+theme_graph()
# on vise l'arête 13-2

# simu données et ajustement modèle
source("F_inference.R")
library(tidyverse)
n=200
Sigma=Graph1_30$sigma
Y=generator_PLN(Sigma,covariates=NULL,n=n)[[1]]

infEM=EMtree(PLNobject = PLN(Y~1))
resampl_inf=ResampleEMtree(Y,cores = 3,S = 50, maxIter=40)
adjmat_freq=F_Vec2Sym(colMeans(1*(resampl_inf$Pmat>0.5)))

graph2=as_tbl_graph(infEM$ProbaCond) %>%  activate(edges)

graph3=adjmat%>%as_tbl_graph() %>%  activate(edges) %>% filter(weight>0.1)

graph3 %>%  ggraph()+
  geom_edge_link(aes(edge_width=weight))+
  geom_node_point()+theme_graph()+
  scale_edge_width_continuous(range=c(0,2))+
  geom_node_text(label = 1:30, color="red", size=5)+theme_graph()

######################
#  1.
# voisins directs de 13 et 20
neighb= graph3 %>% activate(nodes) %>%
  mutate(name=1:30, neighbor=  ( (node_is_adjacent(13) | node_is_adjacent(20)) & !(name %in%c(13,20)) ) ) %>%
  select(neighbor) %>% as.tibble()

indexNeighb= which(neighb==TRUE)
######################
#  2., 3. and 4.
bloc1=1:2
bloc2 = 3:8
Sigma_red=Sigma[c(13,20,indexNeighb),c(13,20,indexNeighb)]
Sigma_partiel = Sigma_red[bloc1,bloc1] - Sigma_red[bloc1,bloc2]%*%solve(Sigma_red[bloc2,bloc2])%*%Sigma_red[bloc2,bloc1]
cov_partielle=Sigma_partiel[1,2]

# toutes les covariances partielles doivent être négatives ici car leur signe est l'opposé des signes
# dans la matric de précision. L'opposition de signe itervient lors de l'inversion

######################
#  Automatisation
# the following functions allow to get the partial covariances between two variables of a given dependency
# graph. This graph is described with an adjacency matrix ; here is tested the matrix of conditional probabilities
# computed by EMtree and the matrix of edges selection frequency computed with ResampleEMtree.

# @get_partialCov
# input : a particular edge in a given network, and a covariance matrix
# output : the partial covariance of this edge conditionnally on its neighborhood
get_partialCov<-function(edge,graph,Sigma){
  from=edge[1]
  to = edge[2]
  p=dim(graph%N>%as_tibble())[1]
  neighb= graph %>% activate(nodes) %>%
    mutate(name=1:p, neighbor=  ( (node_is_adjacent(from) | node_is_adjacent(to)) & !(name %in%c(from,to)) ) ) %>%
    select(neighbor) %>% as.tibble()

  indexNeighb= which(neighb==TRUE)
  v=length(indexNeighb)+2
  if(length(indexNeighb)!=0){
    bloc1=1:2
    bloc2 = 3:v
    reduc=c(from,to,indexNeighb)
    Sigma_red=Sigma[reduc,reduc]
    Sigma_partiel = Sigma_red[bloc1,bloc1] - Sigma_red[bloc1,bloc2]%*%solve(Sigma_red[bloc2,bloc2])%*%Sigma_red[bloc2,bloc1]
    partialCov=Sigma_partiel[1,2]
  }else{
    partialCov=Sigma[from,to]
  }

  return(partialCov)
}

# @edges_partialCov
# input : an adjacency describing a graph matrix that can be weighted, a covariance matrix
#         and a bound to filter the initial weights
# output : a tibble with all edges from the graph, their weight and additional partial covariance
edges_partialCov=function(adjmat,Sigma,lowbound=0.1){
  graph=adjmat%>%as_tbl_graph() %>%  activate(edges) %>% filter(weight>lowbound)
  edges = graph %E>% as_tibble()
  list_edges=edges %>% select(from,to) %>% t() %>% as_tibble() %>% as.list()
  vec_partialCov=unlist(lapply(list_edges,function(x) get_partialCov(x,graph,Sigma)))
  vec_partialCov[which(vec_partialCov>0 & vec_partialCov<1e-15)]=0
  edges$partialCov = vec_partialCov
  print(summary(vec_partialCov))

  return(edges)
}

# @visu_partialCov
# input : an adjacency matrix for a graph, a covariance matrix and a bound to filter edges according to weights
# output : a figure with original graph on the left and infered one on the right, with colors corresponding to the
#       sign of the partial covariances that are computed.
# the filter boolean allows to drop edges with 0 partial covariance. This happen when the true Sigma is
# given.
visu_partialCov<-function(bound,Sigma,adjmat, filter=FALSE, title=""){
  set.seed(1)
  edges=edges_partialCov(adjmat,Sigma,bound)
  edges %>% ggplot(aes(weight,partialCov))+geom_point()+theme_bw()

  original_Layout<-create_layout(graph,"nicely")
  graph_hat=adjmat%>%as_tbl_graph() %>%  activate(edges) %>% filter(weight>bound)

  g=graph%N>%
    mutate(x=original_Layout$x,y=original_Layout$y) %>%
    ggraph(layout="auto")+
    geom_edge_link()+
    geom_node_point()+
    geom_node_text(label = 1:30, color="red", size=5)+theme_graph()

  graph_hat=graph_hat %E>%  mutate(partialCov=edges$partialCov)
  if(filter) graph_hat=graph_hat %>%  filter(partialCov!=0)

  g1= graph_hat %N>%  mutate(x=original_Layout$x,y=original_Layout$y) %>%
    ggraph(layout="auto")+
    geom_edge_link(aes(edge_width=weight, color=as.factor(sign(partialCov))))+
    geom_node_point()+theme_graph()+
    scale_edge_width_continuous(range=c(0.1,2))+
    geom_node_text(label = 1:30, color="blue", size=5)+theme_graph()

  grid.arrange(g,g1,ncol=2,nrow=1, top=title)
}

adjmat_freq1=F_Vec2Sym(colMeans(1*(resampl_inf$Pmat>2/p)))
adjmat_freq2=F_Vec2Sym(colMeans(1*(resampl_inf$Pmat>0.5)))

plnY=PLN(Y~1)
visu_partialCov(2/30,plnY$model_par$Sigma, infEM$ProbaCond,FALSE, title="Conditionel probabilities" ) # probabilities are threshold at 10%
visu_partialCov(0.8, plnY$model_par$Sigma,adjmat_freq1,FALSE, title="Edge selection frequencies") # frequencies are threshold at 80%
visu_partialCov(0.05, plnY$model_par$Sigma,adjmat_freq2,FALSE, title="Edge selection frequencies") # frequencies are threshold at 80%

#####################
# study of partial covariances as function of scores on the edges
partialCov_freq=edges_partialCov(adjmat_freq,plnY$model_par$Sigma,0, TRUE)
partialCov_prob=edges_partialCov( infEM$ProbaCond,plnY$model_par$Sigma,0, TRUE)

partialCov_freq %>% ggplot(aes(weight,partialCov))+geom_point()+geom_hline(yintercept=0,color="red")+
  geom_vline(xintercept=0.5, col="blue")+labs(title="edges selection frequencies")+theme_bw()
partialCov_prob %>% ggplot(aes(weight,partialCov))+geom_point()+geom_hline(yintercept=0,color="red")+
  geom_vline(xintercept=2/30, col="blue")+labs(title="edges conditional probabilities")+theme_bw()
# we are happier with probabilities
# 2/p threshold illustration :
list_prop=lapply(seq(0,1,0.01), function(thresh){
  n_lower=length(which(partialCov_prob$weight<thresh))
  n_upper=length(which(partialCov_prob$weight>thresh))
  prop_low=length(which(partialCov_prob$partialCov>0 & partialCov_prob$weight<thresh))/n_lower
  prop_up=length(which(partialCov_prob$partialCov>0 & partialCov_prob$weight>thresh))/n_upper
  return(c(prop_low,prop_up))
})
props=as_tibble(do.call(rbind,list_prop)[-c(1,101),]);colnames(props)=c("prop_low","prop_up")
props$thresh=seq(0,1,0.01)[-c(1,101)]
props %>% gather(key,value,-thresh) %>%
  ggplot(aes(thresh,value,color=key))+geom_point()+geom_line()+theme_bw()+
  coord_cartesian(xlim = c(0,0.2))+geom_vline(xintercept = 2/30, linetype="dashed")
# we want less positive partial covariance in proportion after the threshold than before it. So we want
# prop_up to be lower than prop_low, which starts to be the case at about 2/p.

#####################
# study of the goodness of prediction as a function of the threshold
original_edges=graph%E>%as_tibble()
join=full_join(partialCov_prob,original_edges, by=c("from","to")) %>%
  mutate(original = !is.na(weight.y))
join %>% ggplot(aes(original,weight.x))+geom_violin()+geom_hline(yintercept=2/30)+theme_bw()

list_prop=lapply(seq(0,1,0.01), function(thresh){
  n_lower=length(which(join$weight.x<thresh))
  n_upper=length(which(join$weight.x>thresh))
  low_true=length(which(join$original==TRUE & join$weight.x<thresh))/n_lower
  high_false=length(which(join$original==FALSE & join$weight.x>thresh))/n_upper
  return(c(low_true,high_false))
})
props=as_tibble(do.call(rbind,list_prop)[-c(1,101),]);colnames(props)=c("FN","FP")
props$thresh=seq(0,1,0.01)[-c(1,101)]
props %>% gather(key,value,-thresh) %>%
  ggplot(aes(thresh,value,color=key))+geom_point()+geom_line()+theme_bw()+
  coord_cartesian(xlim = c(0,0.8))+geom_vline(xintercept = 2/30, linetype="dashed")+geom_hline(yintercept = 0.1)
# for a 50% prob threshold, FP and FN are under 10%

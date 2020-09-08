library(nestor)
library(useful)
library(mclust)
library(gridExtra)
i=1
set.seed(i)
data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=FALSE)
g=EMtree::draw_network(data$G,layout="kk",curv=0,btw_rank=1,groupes=c(rep(1,10),2),
 nodes_size=c(5,5),nodes_label=c(1:10,"H"), pal_nodes= c("#adc9e0","#e7bd42"),
 pal_edges = "#31374f")$G
g
PLNfit<-norm_PLN(data$Y)
Sigma_hat=PLNfit$sigma_obs
#-- original true clique
data$TC
#-- clique found by mclust
clique_mclust=init.mclust(Sigma_hat, r=1)
clique_mclust$init
plotInitMclust(clique_mclust, title="")
i=i+1
pg=grid.arrange(g,p,ncol=2)
ggsave("SFmclust.png", plot=g, width=4, height=4,path= "/Users/raphaellemomal/these/R/images")
ggsave("plotmclust.png", plot=p, width=4, height=4,path= "/Users/raphaellemomal/these/R/images")
ggsave("facetmclust.png", plot=pg, width=6, height=3,path= "/Users/raphaellemomal/these/R/images")

i=i+1

#### blockmodels
i=1
set.seed(7)
data=generate_missing_data(n=100,p=10,r=1,type="cluster",dens=0.3, plot=TRUE)
g=EMtree::draw_network(data$G,layout="kk",curv=0,btw_rank=1,groupes=c(rep(1,10),2),
                       nodes_size=c(5,5),nodes_label=c(1:10,"H"), pal_nodes= c("#adc9e0","#e7bd42"),
                       pal_edges = "#31374f")$G
g
PLNfit<-norm_PLN(data$Y)
Sigma_hat=PLNfit$sigma_O
#-- original true clique
data$TC

init=init.blockmodels(data$Y,Sigma_hat, MO, SO, k=2)

g2<-EMtree::draw_network(init$net,layout="kk",curv=0,btw_rank=1,
                     nodes_size=3, pal_nodes= c("#adc9e0"),
                     pal_edges = "#31374f")$G

init2=init
init2$net[3,4]=init2$net[4,3]=1
init2$net[5,2]=init2$net[2,5]=1
graphorigin=EMtree::draw_network(init$net,layout="kk",curv=0,btw_rank=1,groupes=c(1*(1:10)%in%init$clique[[1]][[1]][[1]]),
                                 nodes_size=c(5,5),nodes_label=1:10, pal_nodes= c("#e7bd42","#adc9e0"),
                                 pal_edges = "#31374f")
lay=create_layout(graphorigin$graph_data,layout="kk")
otherdraw<-function(adj_matrix,title="", size=4, curv=0,width=1, shade=FALSE, filter_deg=FALSE,btw_rank=2,
         layout=NULL,nodes_label=NULL,nodes_size=c(2,5),pal_edges=NULL, pal_nodes=NULL, groupes=NULL){
  adj_matrix=as.matrix(adj_matrix)
  p=nrow(adj_matrix) ; binary=FALSE
  if(is.null(nodes_label)){ nodes_label=1:p ; nonames=TRUE}else{nonames=FALSE}
  if(sum(unique(adj_matrix))==1) binary=TRUE
  # edges width
  min.width=ifelse(binary,0,0.1)
  #edges colour
  pal_edges <-  ifelse(is.null(pal_edges), viridisLite::viridis(5, option = "C")[c(3,2,4,1)], pal_edges)
  #betweenness computations
  res<- as_tbl_graph(adj_matrix, directed=FALSE) %>% activate(edges) %>%
    mutate(btw.weights=ifelse(weight==0,0,log(1+1/weight))) %>%
    activate(nodes) %>%
    mutate( btw=centrality_betweenness(weights=btw.weights),
            bool_btw=(btw>sort(btw, decreasing = TRUE)[btw_rank]),
            bool_deg=(centrality_degree()>0),
            deg=centrality_degree(), title=title, name=nodes_label )
  if(!is.null(groupes)) res<-res %>% mutate(groupes=as.factor(groupes))
  if(nonames){
    res<-res %>% mutate(label=ifelse(bool_btw,name,""))
  }else{
    res<-res %>% mutate(label=ifelse(bool_deg,name,""))
  }
  if(filter_deg) res <- res %>% activate(nodes) %>% filter(deg!=0)
  res<-res %>%
    activate(edges)  %>%
    filter(weight !=0) %>%
    mutate(neibs=edge_is_incident(which(.N()$bool_btw)), title=title)
  
  # define nodes color
  if(!is.null(groupes)){
    if(is.null(pal_nodes)) pal_nodes<-c("#31374f","#adc9e0","#e7bd42")
  }else{
    groupes=unlist(res %>%
                     activate(nodes)  %>% dplyr::select(bool_btw) %>% tibble::as_tibble())
    if(is.null(pal_nodes)) pal_nodes<-c("#31374f","#e7bd42")
  }
  res<-res %>%
    activate(nodes)  %>%
    mutate(finalcolor=groupes)
  
  #draw graph
  set_graph_style(family="sans")
 # layout = ifelse(is.null(layout), "circle", layout)
  g=res %>% ggraph(layout = layout)
  if(shade){ #shading
    g<-g+
      geom_edge_arc(aes(edge_width=weight, alpha=neibs,color = title), strength = curv, show.legend = FALSE) +
      scale_edge_alpha_manual(values=c(0.2,1))
  }else{ g<-g+
    geom_edge_arc(aes(edge_width=weight,color = title), strength = curv, show.legend = FALSE)
  }
  g<-g+
    geom_node_point(aes(color = finalcolor, size = groupes), show.legend = FALSE) +
    scale_edge_colour_manual(values = pal_edges) +
    scale_color_manual(values = pal_nodes)+
    scale_size_manual(values = nodes_size)+
    geom_node_text(aes(label = label), color = "black", size = size) +#,nudge_x = 0.3
    labs(title = title) + theme(plot.title = element_text(hjust = 0.5))+
    scale_edge_width_continuous(range=c(min.width,width))
  
  return(list(G=g,graph_data=res))
}
g3<-otherdraw(init2$net,layout=lay[,1:2],curv=0,btw_rank=1,
                         nodes_size=3, pal_nodes= c("#adc9e0"),
                         pal_edges = "#31374f")$G
g3
p=grid.arrange(g2,g3, ncol=2)
ggsave("chordal.png", plot=p, width=5, height=3.5,path= "/Users/raphaellemomal/these/R/images")

ggimage(init$net[unlist(init$clique),unlist(init$clique)])
        


init.blockmodels<-function(Y, sigma_obs, MO, SO, k=3,poisson=FALSE, alpha=0.1, cores=1){
  init=initVEM(Y = Y,cliqueList=NULL, cov2cor(sigma_obs),MO,r = 0)
  #--- fit nestor with 0 missing actor
  resVEM0<- tryCatch(nestor(Y,MO,SO,initList=init, eps=1e-3, alpha=alpha, maxIter=100,verbatim = 0),
                     error=function(e){e}, finally={})
  # if(length(resVEM0)>3){
  #   sbm.0 <- BM_bernoulli("SBM_sym",1*(resVEM0$Pg>0.5), plotting="", verbosity=0,ncores=cores)
  #   net=1*(resVEM0$Pg>0.5)
  # }else{
    p=ncol(Y)
    resEM0 = EMtree(sigma_obs)
    sbm.0 <- BM_bernoulli("SBM_sym",1*(resEM0$edges_prob>2/p), plotting="", verbosity=0,ncores=1)
    net=1*(resEM0$edges_prob>2/p)
 # }
  sbm.0$estimate()
  paramEstimSBMPoisson <- extractParamBM(sbm.0,k)
  #--- extract k groups from inferred probabilities
  clique=list()
  clique$cliqueList= lapply(1:k, function(z){
    list(which(paramEstimSBMPoisson$Z==z))
  })
  return(list(clique=clique, net=net))
}


##########
# sPCA

set.seed(11)
data=missing_from_scratch(n=100,p=10,r=1,type="scale-free",dens=0.2, plot=TRUE)

data$TC #1236

FitSparsePCA(data$Y,r=3, min.size = 1)$cliques
g=EMtree::draw_network(data$G,layout="kk",curv=0,btw_rank=1,groupes=c(rep(1,10),2),
                     nodes_size=c(5,5),nodes_label=c(1:10,"H"), pal_nodes= c("#adc9e0","#e7bd42"),
                     pal_edges = "#31374f")$G

ggsave("SF_spca.png", plot=g, width=3, height=3,path= "/Users/raphaellemomal/these/R/images")

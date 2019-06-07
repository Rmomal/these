# code pour obtenir les matrices d'adjacence

# 1. matrice S*#(possible edges) des probas de chaque edge pour chaque rééchantillon (S=100)
#     plusieurs modèles -> liste de taille nb models : StabselOak et StabselBarans
#
# 2. Créer les matrices d'adjacence pour chaque modèle des jeux de données, en ajustant le seuil f



# data
StabselOak
StabselBarans

# functions
freq_selec_list<-function(list_Pmat,x,p,f){
  return(F_Vec2Sym( 1*(colMeans( 1*(list_Pmat[[x]]$Pmat>2/p))>f)))
}

adjMatrices<-function(Stabdata,f, models){
  x=ncol(Stabdata[[1]]$Pmat)
  p=(1+sqrt(8*x+1))/2
  
  matrices<-lapply(seq_along(Stabdata), function(x){
    df<-freq_selec_list(Stabdata,x,p,f)
    df
  })
  names(matrices)=models
  return(matrices)
}
#run
f=0.9
OakMatrices=adjMatrices(StabselOak,f,c("Null","Tree","Tree+D1+D2+D3"))
BaransMatrices=adjMatrices(StabselBarans,f,c("Null","Site","Date","Site+Date"))



#####
# package influenceR
library(influenceR)

#exemple du package

ig.ex <- igraph::erdos.renyi.game(100, p.or.m=0.3) # generate an undirected 'igraph' object
keyplayer(ig.ex, k=10) # key-player set consisting of 10 actors


# influenceR est utilisé par tidygraph dans la fonction node_is_keyplayer.
# ça rend un vecteur de booléens, utliisé pour dessiner le réseau.

library(tidyverse)
library(ggraph)
library(tidygraph)

create_Net_Data<-function(Stabdata,f,p,models){ 
  # obtiens les matrices d'adjacence + mise en forme pour tidygraph
  mat<-data.frame(freq_selec_list(Stabdata,1,p,f))
  allNets<-tibble(P = list(mat), models =models )  %>%
    mutate(P=map( seq_along(P), function(x) {
      df<-freq_selec_list(Stabdata,x,p,f)
      df[lower.tri(df, diag = TRUE)]<-0
      df<-data.frame(df)
      colnames(df)<-1:ncol(df)
      df
    })) %>%
    mutate(P = map(P,~rownames_to_column(.) %>%
                     gather(key, value , -rowname))) %>%
    unnest()
  allNets<-allNets[,c(3,2,1,4)]
  return(allNets)
}

models=c("Null","Tree","Tree+D1+D2+D3")

oak_plot<-function(data,seuil){
  allNets=create_Net_Data(data,seuil,p,models)
  groupes<-as.factor(substr(colnames(Y),1,1)) #Y est le tableau des comptages Oak
  set_graph_style()
  mods<-unique(allNets$models)
  
  spliT<-data.frame(allNets) %>% #spliT permet de traiter chaque modèle séparément. On identifie les 
    split(allNets$models) %>%    # key-players, les voisins de EA, plus des manipulations de labels
    tibble(P=map(.,function(x){
      model<-x$models[1]
      res<- as_tbl_graph(x, directed=FALSE) %>%
        activate(edges) %>%
        filter(value !=0) %>% 
        activate(nodes) %>% 
        mutate( groupe=groupes, lab=colnames(Y), keyplayer = node_is_keyplayer(k=6),
                model=model,EA=(name%in%c("48"))) 
      res %>% activate(edges) %>% 
        mutate(neibs=edge_is_incident(which(.N()$EA)), model=model) %>% 
        activate(nodes) %>% 
        mutate(neibEA=name%in%unique(c(.E()$from[which(.E()$neibs)],.E()$to[which(.E()$neibs)]))) %>% 
        mutate(label=ifelse(keyplayer|EA|neibEA,lab,""))
    }))
  # partie graphique : les palettes, ajustement des coordonnées des points,
  #                    jointure des données splitées, et facettage des réseaux.
  pal <- viridisLite::viridis(5, option = "C")
  pal2<-c("gray15","#9de1e1")
  coords<-create_layout(spliT$P[3][[1]],layout="star",center=48)[,c("x","y")]
  spliT$P[1][[1]] %>%
    bind_graphs(spliT$P[2][[1]] )%>%
    bind_graphs(spliT$P[3][[1]] )%>%
    activate(nodes) %>% 
    mutate(model=factor(model,levels=mods), x=rep(coords$x,3),y=rep(coords$y,3)) %>% 
    ggraph(layout="auto")+
    geom_edge_arc(aes(color=model,alpha=neibs),curvature=0.1,show.legend=FALSE)+ 
    scale_edge_alpha_manual(values=c(0.2,1))+
    geom_node_point(aes(color=groupe, size=keyplayer), show.legend=FALSE)+
    scale_edge_colour_manual("model",values=pal[c(2,1,4,3)], labels=mods)+
    scale_color_manual(values=c("indianred1","steelblue4","orange2"))+
    scale_size_manual(values=c(1.5,6))+
    geom_node_text(aes(label = label),color="black")+
    facet_nodes(~model, scales="free")+
    th_foreground(border=FALSE)+
    theme(strip.background = element_rect(fill="white",color="white"),
          strip.text = element_text(color="black",size=12))
}
oak_plot(StabselOak,0.90)
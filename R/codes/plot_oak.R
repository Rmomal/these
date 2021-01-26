library(tidyverse)
library(tidygraph)
library(ggraph)
countEdges<-function(Stabdata,Pt,m_names){
  mat<-data.frame(freq_selec(Stabdata[[1]]$Pmat,Pt)) # the first element is initialized
  
  allNets<-tibble(P = list(mat), mods =m_names )  %>%
    mutate(P=map( seq_along(P), function(x) {
      df<-freq_selec(Stabdata[[x]]$Pmat,Pt=Pt)
      df[lower.tri(df, diag = TRUE)]<-NA
      df<-data.frame(df)
      colnames(df)<-1:ncol(df)
      df
    })) %>%
    mutate(P = map(P,~rownames_to_column(.) %>%
                     gather(key, value , -rowname) %>%filter(!is.na(value))
    ),
    mods=unlist(mods)
    ) %>%
    unnest(cols = c(P))
  allNets<-allNets[,c(1,2,4,3)]
  colnames(allNets) = c("node1","node2","model","weight")
  return(allNets)
}


data.dir = '/Users/raphaellemomal/these/Data_Oak_remmoa/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))
theme_set(theme_bw())
# Parms
Y = as.matrix(Data$count)
colnames(Y)<-unlist(lapply(strsplit(colnames(Y),split="_"),function(x){paste0(toupper(x[1]),x[2])}))
colnames(Y)[48]<-"EA"
models=c("Null","Tree","Tree + D1 + D2 + D3")

########################################################
########################################################
compar_oak<-function(allNets, curv=0.2, width=1, alpha=FALSE,Ft=0,
                     nodes_label=NULL,seed=123, nb=3,groupes=NULL, layout="circle", 
                     base_model=NULL){
  
  mods=unique(allNets$model)
  if(is.null(base_model)) base_model = allNets$model[1]
  nbmod<-length(mods)
  binary=(sum(unique(allNets$weight))==1)
  
  if(is.null(nodes_label)){
    nb_sp = max(as.numeric(allNets$node2))
    nodes_label=1:nb_sp
  }
  
  spliT<-data.frame(allNets) %>%
    base::split(allNets$model) %>%
    tibble(P=map(.,function(x){
      mod<-x$model[1]
      
      res<- as_tbl_graph(x, directed=FALSE) %>%
        activate(edges) %>% filter(weight!=0) %>%
        mutate(btw.weights=log(1+1/weight)) %>%
        activate(nodes) %>%
        mutate( groupe=groupes,importance=centrality_degree(),btw=centrality_betweenness(weights = btw.weights),
                boolbtw=(btw>sort(btw, decreasing = TRUE)[nb]),
                mod=mod,EA=(name%in%c("48")), names=nodes_label) %>%
        activate(edges)  %>% filter(weight>Ft) %>%
        mutate(neibs=edge_is_incident(which(.N()$EA)), mod=mod) %>% 
        activate(nodes) %>% 
        mutate(neibEA=name%in%unique(c(.E()$from[which(.E()$neibs)],.E()$to[which(.E()$neibs)]))) %>% 
        mutate(label=ifelse(boolbtw|EA|neibEA,names,""))
      return(res)
    }))
  
  
  pal_edges <- viridisLite::viridis(5, option = "C")
  pal_nodes<-c("gray15","goldenrod1")
  # lay<-create_layout(spliT$P[[base_model]],layout=layout)
  
  set_graph_style(family="sans")
  
  lay<-create_layout(spliT$P[[base_model]],layout="star",center=48)[,c("x","y")]
  
  
  set.seed(seed)
  plot<- spliT$P %>%
    reduce(bind_graphs) %>%
    activate(nodes) %>%
    mutate(mod=factor(mod,levels=mods),x=rep(lay$x,nbmod),y=rep(lay$y,nbmod)) %>%
    ggraph(layout="nicely")
  if(alpha){
    plot<-plot+
      geom_edge_arc(aes(edge_width=weight,color=mod, alpha=neibs),strength=curv,show.legend=FALSE)+
      scale_edge_alpha_manual(values=c(0.2,1))
    
  }else{plot<-plot+
    geom_edge_arc(aes(edge_width=weight,color=mod),strength=curv,show.legend=FALSE)
  }
  
  g= plot+
    geom_node_point(aes(color=groupe, size=boolbtw), show.legend=FALSE)+
    scale_edge_colour_manual(values=pal_edges[c(2,1,4)], labels=mods)+
    scale_color_manual(values=c("indianred1","orange2"))+
    scale_size_manual(values=c(1.5,6))+
    geom_node_text(aes(label = label),color="black", repel = FALSE,size=4)
  if(length(mods)>1){
    g=g+facet_nodes(~mod, scales="free",ncol=nbmod)+
      th_foreground(border=FALSE)+
      theme(strip.background = element_rect(fill="white",color="white"),
            strip.text = element_text(color="black",size=14))
  }
  
  if(!binary){ g=g+
    scale_edge_width_continuous(range=c(0.1,width))
  }else{
    g=g+
      scale_edge_width_continuous(range=c(0,width))
  }
  return(list(g,spliT))
}
p=ncol(Y)
groupes<-(substr(colnames(Y),1,1))
groupes[48]="F"
groupes=as.factor(groupes)
allNets=countEdges(StabselOak,2/p,models)
graph=compar_oak(allNets =allNets,nb = 9,groupes=groupes, nodes_label =colnames(Y) ,#ifelse(colnames(Y)=="EA","EA",substr(colnames(Y),2,4))
           Ft=0.9,  alpha=TRUE)
graph[[1]]
spliT=graph[[2]]$P

spliT=lapply(spliT,function(graph){
  graph %>% activate(nodes) %>% as_tibble() %>% select(btw,names,mod)
})
spliT=do.call(rbind,spliT)
spliT %>% group_by(mod) %>% summarize(quantEA=pemp(btw[names=="B49"],btw))
#Null                  0.163
#Tree                  0.933
#Tree + D1 + D2 + D3   0.960

ggsave(plot=graph[[1]], filename = "OakProbNets.png",path = "/Users/raphaellemomal/these/R/images",
       height=5.1, width=13)

infectedNet=countEdges(list(finalinftree_covs),2/p,c("infected tree"))
compar_oak(allNets =infectedNet,nb = 7,groupes=groupes, nodes_label =colnames(Y) ,#ifelse(colnames(Y)=="EA","EA",substr(colnames(Y),2,4))
           Ft=0.8,  alpha=TRUE)



########################################################
########################################################
jakuschlabels=c("B1045", "B109", "B1093", "B11", "B112", "B1191", "B123", "B13", "B17", 
                "B171", "B18", "B20", "B22", "B235", "B24", "B27", "B29", "B304", "B31",
                "B33", "B37", "B42", "B443", "B444", "B447", "B47", "B49", "B58", "B59",
                "B625", "B63", "B662", "B69", "B73", "B8", "B87", "B90", "F1", "F10", "F1278", 
                "F15", "F1567", "F19", "F2", "F20", "F23", "F25", "F26", "F28", "F3", "F8", "F9")
jakuschneibs=c("F19","F10","F1","F1278","F28","F15","F2","F1567","F9","F25","F23","F20","F26","B625",
               "B27","B49","B63","B90","B112","B87","B31","B444","B1191","B42","B33","B20")
compare_jak_oak<-function(allNets, nb=5, nodes_label=NULL, Ft=0.8){
  
  mods=unique(allNets$model)
  spliT<-data.frame(allNets) %>%
    base::split(allNets$model) %>%
    tibble(P=map(.,function(x){
      mod<-x$model[1]
      
      res<- as_tbl_graph(x, directed=FALSE) %>%
        activate(edges) %>% filter(weight!=0) %>%
        mutate(btw.weights=log(1+1/weight)) %>%
        activate(nodes) %>%
        mutate( groupe=groupes,importance=centrality_degree(),btw=centrality_betweenness(weights = btw.weights),boolbtw=(btw>sort(btw, decreasing = TRUE)[nb]),
                mod=mod,EA=(name%in%c("48")), names=nodes_label) %>%
        activate(edges)  %>% filter(weight>Ft) %>%
        mutate(neibs=edge_is_incident(which(.N()$EA)), mod=mod) %>% 
        activate(nodes) %>% 
        mutate(neibEA=name%in%unique(c(.E()$from[which(.E()$neibs)],.E()$to[which(.E()$neibs)]))) %>% 
        mutate(label=ifelse(boolbtw|EA|neibEA,names,""))
      return(res)
    }))
  
  # datatest=spliT %>% 
  #   mutate( neibors=map(P,function(x){
  #     nodes =x %N>% as.tibble() %>% filter(names!="EA")
  #     res=nodes[nodes$neibEA,c("names","mod")]
  #     return(res) }))
  # 
  # datajak=data.frame(names=jakuschneibs, jak="jak")
  # final_join= full_join(datatest$neibors[[1]],datajak, by="names")
  # final_join= final_join %>% mutate(EMtree=ifelse(is.na(mod),0,1), jak=ifelse(is.na(jak),0,1))
  # 
  # res=data.frame(`00`= length(which((final_join$EMtree+final_join$jak)==0)),
  #                    `11` = length(which((final_join$EMtree+final_join$jak)==2)),
  #                    `01`=length(which(final_join$EMtree==0 & final_join$jak==1)),
  #                    `10`=length(which(final_join$EMtree==1 & final_join$jak==0)))
  # 
  return(spliT$P)
}

########################################################
#infected with covariates
finalinftree_covs=readRDS(paste0(data.dir,"infectedTree.rds"))
prob=data.frame(x=finalinftree_covs$Pmat[1,])
prob %>% ggplot(aes(x))+ geom_histogram(binwidth = 0.005,alpha=0.4, fill="steelblue", color="steelblue")+theme_light()+
  coord_cartesian(xlim=c(-0.001,0.1))+
  geom_vline(xintercept = 2/114, color="red", linetype="dashed")+
  labs(x="Probabilities", y="")
freqs=freq_selec(finalinftree_covs$Pmat,Pt=2/114)
hist(F_Sym2Vec(freqs), breaks=30,prob=TRUE)
hist(finalinftree_covs$Pmat[1,], breaks=30, prob=TRUE)
ggimage(F_Vec2Sym(finalinftree_covs$Pmat[1,]))
ggimage(F_Vec2Sym(finalinftree_covs$Pmat[1,]>2/114))
ggimage(freqs)
ggimage(1*(freqs>0.9))
infectedNet=countEdges(list(finalinftree_covs),2/p,c("infected tree"))
reseaux=compare_jak_oak(infectedNet, nodes_label = colnames(Y), Ft=0)

proba_neigh_infected=reseaux[["infected tree"]] %>% activate(edges) %>% as_tibble() %>% 
  filter(neibs) %>% select(from,to,weight) %>% 
  mutate(neighb_node=ifelse(from==48,to,from)) %>% mutate(label=colnames(Y)[neighb_node]) %>% 
  mutate(injak=as.factor(1*label%in%jakuschneibs),injaklarge=as.factor(1*label%in%jakuschlabels)) %>%
  select(weight,injak,injaklarge,label)
proba_neigh_infected %>% ggplot(aes(as.factor(injak),weight))+geom_beeswarm()+theme_bw()

#infected without covariates
infectedNet_null=countEdges(list(inftree_null),2/p,c("infected tree no covar"))
reseau_null=compare_jak_oak(infectedNet_null, nodes_label = colnames(Y), Ft=0)

proba_neigh_infected_null=reseau_null[["infected tree no covar"]] %>% activate(edges) %>% as_tibble() %>% 
  filter(neibs) %>% select(from,to,weight) %>% 
  mutate(neighb_node=ifelse(from==48,to,from)) %>% mutate(label=colnames(Y)[neighb_node]) %>% 
  mutate(injak=as.factor(1*label%in%jakuschneibs),injaklarge=as.factor(1*label%in%jakuschlabels) )%>% 
  select(weight,injak,injaklarge,label)

proba_neigh_infected_null %>% ggplot(aes(as.factor(injak),weight))+geom_beeswarm()+theme_bw()





proba_neigh=rbind(cbind(proba_neigh_infected,Covariates="~ leaf \nposition"),
                  cbind(proba_neigh_infected_null,Covariates="~1")) %>% 
  mutate(injak=fct_recode(injak,"EA neighbors"="1","Others"="0"))


proba_neigh_infected=cbind(proba_neigh_infected,Covariates="~ leaf \nposition") %>% 
  mutate(injak=fct_recode(injak,"Ea neighbors"="1","Others"="0"))
EAneighbors=proba_neigh_infected %>% ggplot(aes(fct_reorder(injak,-weight),weight,color=fct_reorder2(Covariates, injak,-weight)))+
  geom_boxplot(width=0.4)+labs(x="Jakuschkin et al.",y="Ea neighbors selection frequency")+
  scale_color_brewer("",palette="Dark2")+
  theme_light()+ coord_flip()

proba_neigh %>% ggplot(aes(injaklarge,weight,color=Covariates))+
  geom_boxplot(width=0.3)+labs(x="Jakuschkin et al.",y="EA neighbors probability")+
  scale_color_brewer(palette="Dark2")+
  theme_light()


proba_neigh %>% filter(Covariates!="~1" & weight > 0.15 & injak=="Other")
proba_neigh %>% group_by(Covariates,injak) %>% summarise(min=min(weight),med=median(weight),
                                                                          max=max(weight))
ggsave(plot=EAneighbors,filename = "EAneighbors.png",path="/Users/raphaellemomal/these/R/images",
       height=2.5,width=5.5)
########################################################
########################################################

datapct=lapply(seq(0,0.99,0.02),function(x){
  return(cbind(compare_jak_oak(infectedNet, nodes_label = colnames(Y),Ft=x),"seuil"=x))
})
datapct=do.call(rbind,datapct)
dataplot=datapct %>% as_tibble() %>% 
  mutate(sumEMtree=X10+X11, sumjak=X01+X11, FDR=X10/sumEMtree, PPVEMtree=X11/sumEMtree,
         PPVjak=X11/sumjak) 



plot=dataplot %>%  select(-PPVEMtree,-PPVjak,-FDR,-sumjak,-sumEMtree) %>% gather(type, value, -seuil) %>% 
  ggplot( aes(seuil, value, col=type))+geom_point()+geom_line()+theme_minimal()+
  labs(x="Selection threshold",y="Number of OTUs",color = "EMtree vs.\nJakuschkin \net al.:")+
  geom_vline(xintercept = 0.9, linetype="dashed")+scale_color_brewer(palette="Dark2")
ggsave("curves_OTUS.png",plot=plot,height=5,width=7.5 )


brewer.pal(8,"Dark2")

plot2=dataplot %>%  select(PPVEMtree,seuil) %>% gather(type, value, -seuil) %>% 
  ggplot( aes(seuil, value, col=type))+geom_point()+geom_line()+theme_minimal()+
  scale_color_manual(name="",labels="(OTUs in common)/\n(EMtree discoveries)", values="#66A61E")+
  geom_vline(xintercept = 0.9, linetype="dashed")+
  labs(x="Selection threshold",y="%")

plot3=grid.arrange(plot,plot2, nrow=1, ncol=2)
ggsave("propOTUS.png",plot=plot2,height=4,width=6)



ggsave("Oaknet_90_btw.png", width=11, height=4.5)
################
# compare tables accross models
load("/Users/raphaellemomal/these/Data/Oaks-CVacher/oaks-StabSel.Rdata")
null1<-Stab.sel[[1]][[1]][[1]] # comparaison des proba non seuillées pour 1 sous échantillon
null2<-Stab.sel[[1]][[2]][[1]]
null3<-Stab.sel[[1]][[3]][[1]]
null4<-Stab.sel[[1]][[4]][[1]]
datafreq<-tibble(it1=colMeans(null1),it2=colMeans(null2),it3=colMeans(null3),it4=colMeans(null4))
datafreq %>% 
  ggplot(aes(x=it4))+
  geom_point(aes(y=it2),color="orange", size=0.4)+
  geom_point(aes(y=it1),color="deepskyblue", size=0.4)+
  geom_point(aes(y=it3),color="red", size=0.6)+
  geom_abline(linetype="dashed")+
  geom_hline(yintercept = 2/p)+
  geom_vline(xintercept = 2/p)+
  coord_cartesian(xlim=c(0,0.08), ylim=c(0,0.16))




tree<-1*(colMeans( 1*(Stab.sel[[2]]$Pmat>2/p))>freq.sel.thres)
all<-1*(colMeans( 1*(Stab.sel[[3]]$Pmat>2/p))>freq.sel.thres)

table(null,tree)
table(null,all)
table(tree,all)
table(tree,all,null) # 8 en commun aux trois
table(site, date_Site)



summary(Stab.sel[[1]]$iter)
summary(Stab.sel[[2]]$iter)
summary(Stab.sel[[3]]$iter)

deg1<-colSums(F_Vec2Sym(null))
deg2<-colSums(F_Vec2Sym(tree))
deg3<-colSums(F_Vec2Sym(all))

summary(deg3-deg1)

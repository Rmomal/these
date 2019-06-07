
# nombre d'arêtes par modèle en fonction du seuil de sélection
countEdges<-function(Stabdata,f,models){
  x=ncol(Stabdata[[1]]$Pmat)
  p=(1+sqrt(8*x+1))/2
  
  mat<-data.frame(freq_selec_list(Stabdata,1,p,f))
  #models=names(covar)
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
  res<-allNets %>% group_by(models) %>% summarise(sum=sum(value))
  res$seuil=f
  return(res)
}
countEdges_notEMtree<-function(Stabdata,f){
  #browser()
  res<-c(sum=sum((Stabdata[upper.tri(Stabdata, diag = FALSE)]>f)*1),seuil=f)
  
  return(res)
}
load(file = paste0(data.dir, data.name, '-StabselOak_allSpecies.Rdata'))

# StabselOak_shuffle<-(c(StabselOak,shuffleOakTree))
# StabselBarans_shuffleDate<-c(Faatfreq[[1]]$`100`,shuffleBaransDate)
# StabselBarans_shuffleSite<-c(Faatfreq[[1]]$`100`,shuffleBaransSite)

# plot de l'évolution de la densité d'arête en fonction du seuil, et accolage d'une légend agrandie
evolPlot<-function(Stabdata, mod, labs=NULL, orderpal,ordershapes){
# browser() 
  if(is.null(labs)) labs=mod
  nmods=length(mod)-1
  pal=c(pal,"#274BF3")
  shapes=c(1,8,15,16,17)
  countsEdges<-do.call(rbind,lapply(seq(0,1,0.02), function(x) countEdges(Stabdata, x, models=mod)))
 setseuil=unique(countsEdges$seuil)
 # Oak_shuffle_selec=Oak_shuffle_tree %>% filter(seuil %in% setseuil, type=="shuffled") %>% select(-"type")
 #  countsEdges=rbind(countsEdges,Oak_shuffle_selec)
  ylimit=max(countsEdges$sum[countsEdges$seuil>0.75])

  plot<-countsEdges %>% as_tibble() %>%# mutate( models=rename_factor(as.factor(models),  "null"="~ 1","tree"="~ tree","treeDist"="~ tree + D1 + D2 + D3")) %>% 
    ggplot(aes(seuil, sum, color=models, shape=models))+ theme_minimal()+geom_line(alpha=0.4)+
    labs(x="Selection threshold", y="Quantity of edges")+ theme(axis.text=element_text(size=12),
                                                                axis.title=element_text(size=12), legend.position="bottom",
                                                                legend.text = element_text(size=12))+
    scale_colour_manual("",values=pal[orderpal], breaks=unique(countsEdges$models),labels=labs)+
    scale_shape_manual("",values=shapes[ordershapes], breaks=unique(countsEdges$models),labels=labs)

  legend=g_legend(plot+geom_point(size=2.2,alpha=0.8))
  plot=plot+geom_point(alpha=0.8)+guides(color=FALSE, shape=FALSE)#+coord_cartesian(xlim=c(0.8,1), ylim=c(0,ylimit))
  plot<-grid.arrange(plot,legend, nrow=2, ncol=1,heights=c(9,1))
}


EvolOak<-evolPlot(StabselOak, mod=c("null","tree","treeDist"),
                  labs=c("Null","Tree", "Tree+D1+D2+D3"),orderpal=c(2,1,4),
                  ordershapes=c(4,3,5))
ggsave("QET_Oak.png",plot=EvolOak, width=5, height=5.5)


EvolBaransSite<-evolPlot(StabselBarans_shuffleSite,  c("Null","Site","Date","Site+Date","shuffleSite"),
                         labs=c("Null","Site","Date","Site+Date","Shuffle Site"),c(1,2,6,4,3),c(4,3,1,5,2))
EvolBaransDate<-evolPlot(StabselBarans_shuffleDate,  c("Null","Site","Date","Site+Date","shuffleDate"),
                         labs=c("Null","Site","Date","Site+Date","Shuffle Date"),c(1,2,6,4,3),c(4,3,1,5,2))

#plotShuffle<-evolPlot(shuffleEMtree,  c("Null","Site","Date","Site+Date"))
assembl=grid.arrange(EvolBaransSite,EvolBaransDate,EvolOak,nrow=1,ncol=3)
ggsave("QET_shuffletwoDataSets.png",plot=assembl, width=14, height=5)
ggsave("QET_zoomShuffle.png",plot=assembl, width=14, height=5)


##################

models=c("Null","Site","Date","Site+Date")


countsEdgesBarans<-do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges(Faatfreq[[1]]$`100`, x, models=models)))
countsEdgesShuffle<-do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges(shuffleEMtree, x, models=models)))
countsEdgesBarans$type="original"
countsEdgesShuffle$type="shuffled"
countsEdgesShuffle_full<-do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges(shuffleNull, x, models="Null")))
countsEdgesShuffle_full$type="shuffled"
countsBind=rbind(countsEdgesBarans[which(countsEdgesBarans$models=="Null"),],countsEdgesShuffle_full)

countsBind %>% as_tibble() %>%# mutate( models=rename_factor(as.factor(models),  "null"="~ 1","tree"="~ tree","treeDist"="~ tree + D1 + D2 + D3")) %>% 
  ggplot(aes(seuil, sum, color=models, shape=type))+ theme_minimal()+geom_line(alpha=0.4)+
  labs(x="Selection threshold", y="Quantity of edges")+ theme(axis.text=element_text(size=12),
                                                              axis.title=element_text(size=12), legend.position="bottom",
                                                              legend.text = element_text(size=12))+
  geom_point()#+coord_cartesian(xlim=c(0.8,1),ylim=c(0,60))#+
scale_colour_manual("",values=pal[c(1,2,4,3)])+ scale_shape_manual("",values=c(15,16,17,8))

##########################################
##########################################
#mélanger et inférer

fillFreq<-function(counts, models){
  # models<- list("1")
  res<-lapply(models, function(model){
    ResampleEMtree(counts,model,covariate=NULL, S=100, maxIter=300,cond.tol=1e-12, cores=1)
  })
  return(res)
}
shuffle<-function(counts,var){
  levels<-levels(var)
  counts_shuffle=0*counts
  for (bloc in levels){
    indices<-which(var==bloc)
    counts_shuffle[indices,]<-apply(counts[indices,], 2, function(x){
      sample(x)
    })
  }
  return(counts_shuffle)
}
# shuffled data
data_shuffleTree<-shuffle(Y,X$tree)

# inferences on shuffled data
# EMtree
shuffled_EMtree<-function(data=Y,var=X$tree){
  textvar<-deparse(substitute(var))
  model=strsplit(textvar,split="\\$")[[1]][2]
  data_shuffle<-shuffle(data,var)
  print("a")
  shuffleOak<-fillFreq(data_shuffle,textvar) 
  print("b")
  OakcountsEdgesShuffle<-do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges(shuffleOak, x,models=model)))
  print("c")
  OakcountsEdgesShuffle$type="shuffled"
  OakcountsEdges<-do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges(StabselOak, x, models=model)))
  OakcountsEdges$type="original"
  OakcountsBind=rbind(OakcountsEdges,OakcountsEdgesShuffle)
  return(OakcountsBind)
}
Oak_shuffle_tree<-shuffled_EMtree()
Oak_shuffle_null<-shuffled_EMtree(var=1)
Oak_shuffle_tree$models="Tree shuffled"
#le tableau des fréquences de sélection de chacun des arêtes pour les données oaks shufflées, 
#ainsi que les données originales,



#gCoda
shuffle_gCoda<-function(data1,data2,shuffleRes1,shuffleRes2, covar="X$tree", shufflevar){
  originalOak_gCoda<-gcoda(data1, counts=T, covar=covar, shuffleRes=shuffleRes1)
  original_gCoda_score <- Reduce("+",originalOak_gCoda$path)
  original_gCoda_score<- original_gCoda_score / max(original_gCoda_score)
  OakcountsEdgesoriginal_tree_gcoda<-as_tibble(do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges_notEMtree(original_gCoda_score, x))))
  OakcountsEdgesoriginal_tree_gcoda$type="original"
  
  shuffleOak_gCoda<-gcoda(data2, counts=T, covar=covar, shuffleRes=shuffleRes2, shufflevar =shufflevar)
  shuffle_gCoda_score <- Reduce("+",shuffleOak_gCoda$path)
  shuffle_gCoda_score<- shuffle_gCoda_score / max(shuffle_gCoda_score)
  OakcountsEdgesShuffle_tree_gcoda<-as_tibble(do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges_notEMtree(shuffle_gCoda_score, x))))
  OakcountsEdgesShuffle_tree_gcoda$type="shuffled"
  
  OakcountsBind_gCoda=rbind(OakcountsEdgesoriginal_tree_gcoda,OakcountsEdgesShuffle_tree_gcoda)
  return(OakcountsBind_gCoda)
}
gCoda_shuffle<-shuffle_gCoda(data1=Y, data2=data_shuffleTree, shuffleRes=FALSE )
#spiecEasi
shuffle_spiecEasi<-function(data1,data2,shuffleRes, covar=X, shufflevar,reg=TRUE){
  U<-t(clr.matrix(data1,mar=1))
  m<- model.matrix(~tree,covar)
  resid<-lm(U~m)$residuals
  if(!reg) resid=U
  originalOak_spiec<- spiec.easi(resid, method="glasso",icov.select = FALSE, nlambda = 50, verbose = FALSE)
  original_spiec_score <- Reduce("+",originalOak_spiec$est$path)
  original_spiec_score<- original_spiec_score / max(original_spiec_score)
  OakcountsEdgesoriginal_tree_spiec<-as_tibble(do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges_notEMtree(original_spiec_score, x))))
  OakcountsEdgesoriginal_tree_spiec$type="original"
  
  U<-t(clr.matrix(data2,mar=1))
  resid<-lm(U~m)$residuals
  if(shuffleRes) resid=shuffle(resid,shufflevar)
  if(!reg) resid=U
  shuffleOak_spiec<- spiec.easi(resid, method="glasso",icov.select = FALSE, nlambda = 50, verbose = FALSE)
  shuffle_spiec_score <- Reduce("+",originalOak_spiec$est$path)
  shuffle_spiec_score<- shuffle_spiec_score / max(shuffle_spiec_score)
  OakcountsEdgesShuffle_tree_spiec<-as_tibble(do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges_notEMtree(shuffle_spiec_score, x))))
  OakcountsEdgesShuffle_tree_spiec$type="shuffled"
  OakcountsBind_spiec=rbind(OakcountsEdgesoriginal_tree_spiec,OakcountsEdgesShuffle_tree_spiec)
}
spiecEasi_shuffle<-shuffle_spiecEasi(data1=Y, data2=data_shuffleTree,shuffleRes=FALSE,shufflevar=X$tree )

#####################
# tests de vérification sur gCoda et SpicEasi avec les données Oak
# T1 : mélanger les résidus et au lieu des données (on doit retrouver le même résultat)
# T2 : donner les données non régressées (on veut voir moins d'arêtes inférées avec les données mélangées, même 
#      sans prise en compte de tree)

# gCoda
gCoda_T1<-shuffle_gCoda(data1=Y, data2=Y, covar="X$tree",shuffleRes1=FALSE,shuffleRes2=TRUE, shufflevar=X$tree)
gCoda_T2<-shuffle_gCoda(data1=Y, data2=data_shuffleTree, covar=NULL,shuffleRes1=FALSE,shuffleRes2=FALSE)

#spiec
spiecEasi_T1<-shuffle_spiecEasi(data1=Y, data2=Y,shuffleRes=TRUE,shufflevar=X$tree )
spiecEasi_T2<-shuffle_spiecEasi(data1=Y, data2=data_shuffleTree,shuffleRes=FALSE,shufflevar=X$tree, reg=FALSE )


##########################################
##########################################
# graphiques QET
QETplot<-function(data,coordx=NULL, main=NULL){
  pal <- viridisLite::viridis(5, option = "C")
  g<-data %>% as_tibble() %>%# mutate( models=rename_factor(as.factor(models),  "null"="~ 1","tree"="~ tree","treeDist"="~ tree + D1 + D2 + D3")) %>% 
    ggplot(aes(seuil, sum, color=type, shape=type))+ theme_minimal()+geom_line(alpha=0.4)+
    labs(x="Selection threshold", y="Quantity of edges")+ theme(axis.text=element_text(size=12),
                                                                axis.title=element_text(size=12),
                                                                legend.position="bottom",
                                                                legend.text = element_text(size=12))+
    geom_point()+scale_colour_manual("",values=pal[c(1,2,4,3)])+ scale_shape_manual("",values=c(15,16,17,8))+labs(title=main)
  if(!is.null(coordx)){
    yl=max(data$sum[data$seuil==coordx])
    xl=max(data$seuil)
    g<-g+coord_cartesian(xlim=c(coordx,xl),ylim=c(0,yl))
  }
  return(g)
}

x=NULL
q1<-QETplot(OakcountsBind,coordx=x)
q2<-QETplot(gCoda_T2,coordx=x)
q3<-QETplot(OakcountsBind_spiec,coordx=x)
grid_arrange_shared_legend(q1,q2,q3, nrow=1, ncol=3)

q1<-QETplot(gCoda_shuffle,coordx=x, main="gCoda sur résidus originaux \nvs. sur résidus des données mélangées")
q2<-QETplot(gCoda_T1,coordx=x, main="gCoda sur résidus originaux \nvs.sur résidus mélangés")
q3<-QETplot(gCoda_T2,coordx=x, main="gCoda sur données brutes \nvs. sur données mélangés (pas d'ajustement)")
plot=grid_arrange_shared_legend(q1,q2,q3, nrow=1, ncol=3)
ggsave("QET_gCoda.png",plot=plot,width=14, height=6)


q1<-QETplot(spiecEasi_shuffle,coordx=x, main="spiecEasi sur résidus originaux \nvs. sur résidus des données mélangées")
q2<-QETplot(spiecEasi_T1,coordx=x, main="spiecEasi sur résidus originaux \nvs.sur résidus mélangés")
q3<-QETplot(spiecEasi_T2,coordx=x, main="spiecEasi sur données brutes \nvs. sur données mélangés (pas d'ajustement)")
plot=grid_arrange_shared_legend(q1,q2,q3, nrow=1, ncol=3)
##########################################
##########################################
#gCOda, SpiecEasi et MInt optimaux sur les données Oak
# à faire : modifier la grille des lambdas
#gCoda
covarexp<-X[,1:2]
colnames(covarexp)<-c("X1","X2")
infShufflegCoda <- F_Sym2Vec(gcoda(Oakcounts_fullshuffle, counts=T)$opt.icov>1e-16)*1
sum(infShufflegCoda)
#MInt
covarexp<-covarexp %>% mutate(X1=as.numeric(X1), X2=as.numeric(X2))
inShuffleMInt<- F_Sym2Vec(eval_store_mint(Oakcounts_fullshuffle,covarexp,path)>1e-16)*1
sum(inShuffleMInt)
#SpiecEasi
m<- model.matrix(~X1+X2,covarexp)
U<-t(clr.matrix(counts_shuffle,mar=1))
model<-lm(U~m)
inShuffleSpiec<-spiec.easi(Oakcounts_fullshuffle, method='mb', lambda.min.ratio=1e-2, nlambda=30,icov.select.params=list(rep.num=100 ))
spiecgraph<-F_Sym2Vec(as.matrix(inShuffleSpiec$refit[[1]]))
sum(spiecgraph)


##########################################
##########################################
# mélanger des données simulées:
# a simulated cluster exmaple
Graph1_30 <- readRDS("/Users/raphaellemomal/simulations/Simu/PLN_nonfav/cluster/d/Sets_param/Graph1_30.rds")
fixeCovar <- readRDS("/Users/raphaellemomal/simulations/Simu/PLN_nonfav/fixeCovar.rds")
colnames(fixeCovar)[3]="factor"#la 3e est un facteur
simulDataClust<-generator_PLN(Graph1_30$sigma,fixeCovar) #rend Y et cor(Z)

# shuffling on simulated data
shuffleDataClust<-shuffle(simulDataClust[[1]],fixeCovar$factor)

# inferences on shuffled data
infshuffleDataClust<-fillFreq(shuffleDataClust,"fixeCovar$factor")




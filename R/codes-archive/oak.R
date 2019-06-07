rm(list=ls()); par(pch=20);
#devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
library(PLNmodels);
library(igraph)
library(RColorBrewer)
library(tidyverse)
library(parallel)
source('/Users/raphaellemomal/simulations/codes/FunctionsMatVec.R')
source('/Users/raphaellemomal/simulations/codes/FunctionsTree.R')
source('/Users/raphaellemomal/simulations/codes/FunctionsInference.R')
#source('TreeMixture-mac.R')
#source('/Users/raphaellemomal/simulations/codes/TreeMixture-RML.R')
source('/Users/raphaellemomal/simulations/codes/fonctions.R')
# Data
data.dir = '/Users/raphaellemomal/these/Data/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))
theme_set(theme_bw())
# Parms
Y = as.matrix(Data$count);
O = Data$offset; X = Data$covariates
# Selection des especes
YO = Y/O
Rank = rank(colSums(YO))
#plot(cumsum(sort(colSums(YO))))
Seuil = 20
Ycum = colSums(Y); Order = order(Ycum)
#plot(cumsum(Ycum[Order]), col = 1+(Rank[Order]<Seuil))

Y = Y[, Rank > Seuil];
O = O[, Rank >Seuil];
n = nrow(Y); p = ncol(Y)
colnames(Y)<-unlist(lapply(strsplit(colnames(Y),split="_"),function(x){paste0(toupper(x[1]),x[2])}))
colnames(Y)[44]<-"EA"
#################################################################################
#################################################################################

# Inférences
PLN.vide<-PLN(Y ~ 1 + offset(log(O)))
PLN.tree<-PLN(Y ~ 1 + X$tree+  offset(log(O)))
PLN.orient<-PLN(Y ~ 1 + X$orientation+  offset(log(O)))
PLN.dist1<-PLN(Y ~ 1 + X$distTOground+  offset(log(O)))
PLN.dist2<-PLN(Y ~ 1 + X$distTObase+  offset(log(O)))
PLN.dist3<-PLN(Y ~ 1 + X$distTOtrunk+  offset(log(O)))
PLN.distOrient<-PLN(Y ~ 1 + X$distTOground+ X$orientation+  offset(log(O)))
PLN.infect<-PLN(Y ~ 1 +  X$pmInfection+  offset(log(O)))
PLN.infectDist<-PLN(Y ~ 1 +  X$pmInfection + X$distTOground+  offset(log(O)))
PLN.infectAllDist<-PLN(Y ~ 1 +  X$pmInfection +X$distTObase+X$distTOtrunk+X$distTOground+   offset(log(O)))
PLN.infectTree<-PLN(Y ~ 1 +  X$pmInfection + X$tree+  offset(log(O)))
PLN.infectOrient<-PLN(Y ~ 1 +  X$pmInfection + X$orientation+  offset(log(O)))
PLN.infectTreeOrient<-PLN(Y ~ 1 +  X$pmInfection + X$tree+X$orientation+  offset(log(O)))
PLN.Alldist<-PLN(Y ~ 1 +  X$distTObase+X$distTOtrunk+X$distTOground+ offset(log(O)))
PLN.treeOrient<-PLN(Y ~ 1 + X$tree+ X$orientation+  offset(log(O)))
PLN.treeAlldist<-PLN(Y ~ 1 +  X$tree+ X$distTObase+X$distTOtrunk+X$distTOground+ offset(log(O)))

PLN.treeDist<-PLN(Y ~ 1 +
                    X$tree+
                    X$distTOground+
                    offset(log(O)))
PLN.treeallDist<-PLN(Y ~ 1 +
                       X$tree+
                       X$distTObase+X$distTOtrunk+X$distTOground+
                       offset(log(O)))
PLN.allcovariates<-PLN(Y ~ 1 +
                         X$tree+
                         X$distTObase+X$distTOtrunk+X$distTOground+
                         X$pmInfection+X$orientation+ 
                         offset(log(O)))

datamodels<-data.frame(cbind(rbind(PLN.vide$criteria,
                                   PLN.tree$criteria,
                                   PLN.orient$criteria,
                                   PLN.dist1$criteria,
                                   PLN.dist2$criteria,
                                   PLN.dist3$criteria,
                                   PLN.treeDist$criteria,
                                   PLN.treeOrient$criteria,
                                   PLN.Alldist$criteria,
                                   PLN.treeDist$criteria,
                                   PLN.allcovariates$criteria,
                                   PLN.infect$criteria,
                                   PLN.infectDist$criteria,
                                   PLN.infectAllDist$criteria,
                                   PLN.infectTree$criteria,
                                   PLN.infectOrient$criteria,
                                   PLN.infectTreeOrient$criteria,
                                   PLN.treeAlldist$criteria,
                                   PLN.distOrient$criteria),
                             model=c("vide","tree","orient","dist1","dist2","dist3","treeDist",
                                     "treeOrient","alldist","treeDist","allcovariates",
                                     "infect","infectDist","infectAllDist", "infectTree", "infectOrient",
                                     "infectTreeOrient","TreeAllDist","distOrient")))

datamodels[,1:5]<-apply(datamodels[,1:5],2,function(x) as.numeric(as.character(x)))
p1<-ggplot(datamodels, aes(x = reorder(model, loglik), y = loglik)) +
  geom_point()+
  # geom_point(aes(y=ICL),color="red")+
  # geom_point(aes(y=BIC),color="blue")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2<-ggplot(datamodels, aes(x = reorder(model, -ICL), y = ICL)) +
  geom_point(color="red")+
  # geom_point(aes(y=loglik))+
  # geom_point(aes(y=BIC),color="blue")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p3<-ggplot(datamodels, aes(x = reorder(model, -BIC), y = BIC)) +
  geom_point(color="blue")+
  # geom_point(aes(y=ICL),color="red")+
  # geom_point(aes(y=loglik))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(p1,p2,p3,nrow=3,ncol=1)
png(filename="~/these/fig.png", width=4, height=7, units="in",res=300)
grid.arrange(p1,p2,p3,nrow=3,ncol=1)
dev.off()
# datamodels<-datamodels %>% 
#   gather(key, value, -"model")
# ggplot(datamodels,aes(model))+
#   geom_col(aes(y=value,fill=key))

Z.vide = PLN.vide$model_par$Sigma
Z.tree=PLN.tree$model_par$Sigma
Z.treeDist=PLN.treeAlldist$model_par$Sigma



# trace de la diff, de la loglik. affichages des histogrammes de beta et de leur summary (save des summary)
# 1. trajectoires vraisemblances EM
track_convergence<-function(inf.vide,title){
  loglik<-data.frame(logL=inf.vide$L,diff=inf.vide$diff_beta)
  
  p1<-loglik %>% 
    mutate(index=1:nrow(loglik)) %>% 
    ggplot(aes(x=index,y=logL))+
    geom_point(color="deepskyblue")+
    geom_line()+
    theme_bw()
  p2<-loglik %>% 
    mutate(index=1:nrow(loglik)) %>% 
    ggplot(aes(x=index,y=diff))+
    geom_point(color="deepskyblue")+
    geom_line()+
    theme_bw()
  p3<-loglik %>% 
    mutate(index=1:nrow(loglik)) %>% 
    ggplot(aes(x=index,y=diff))+
    geom_point(color="deepskyblue")+
    geom_line()+
    coord_cartesian(ylim = c(0,2e-4))+
    geom_abline(intercept = 1e-4,slope=0,color="deepskyblue",size=1.3)+
    geom_abline(intercept = 1e-6,slope=0,color="deepskyblue3")+
    theme_bw()+
    labs(y="diff zoom1")
  p4<-loglik %>% 
    mutate(index=1:nrow(loglik)) %>% 
    ggplot(aes(x=index,y=diff))+
    geom_point(color="deepskyblue")+
    geom_line()+
    coord_cartesian(ylim = c(0,2e-6))+
    geom_abline(intercept = 1e-6,slope=0,color="deepskyblue3", size=1.3)+
    theme_bw()+
    labs(y="diff zoom2")
  grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2,top=title)
}

T1<-Sys.time()
inf.vide<-TreeGGM(cov2cor(Z.vide),print=FALSE,step="FALSE",maxIter = 150, n=nrow(Y) )
# itermax=154, 13.3 min // alpha=0.1034, 6.13 mins itermax=22 mais NAN // 1.34e_17 // NAN et Error in while ((diff > eps) & (iter < maxIter)) { : 
# missing value where TRUE/FALSE needed
# tolérance trop faible ?
#// 5.7 min, itermax=24
T2<-Sys.time()
difftime(T2,T1)
T1<-Sys.time()
inf.sansnorm<-TreeGGM(cov2cor(Z.tree),print=FALSE,step="FALSE",maxIter = 150,n=nrow(Y),
                      cond.tol=1e-12)
# itermax = 243 20min // alpha = 0.224, 2.56 min, itermax=20, et again NaN 1.497385e-17 NaN 
# Error in while ((diff > eps) & (iter < maxIter)) { : 
#     missing value where TRUE/FALSE needed
#// 2.6 mins, itermax=21
T2<-Sys.time()
difftime(T2,T1)

#inf.betamin1, 40, 500.26
# inf.betamin2, 150, 500.26
# inf.norm1, 37, 488.6
# inf.norm2, 42, 489

# impact sur les valeurs de proba ?
# attendu : différence entre betalin1 et betamin2, mais pas entre norm1 et norm2
save(inf.betamin1,inf.betamin2,inf.norm1,inf.norm2,inf.norm3,file="testsBetas")
summary(diff(F_Sym2Vec(inf.betamin1$P)-F_Sym2Vec(inf.betamin2$P)))
summary(diff(F_Sym2Vec(inf.norm1$P)-F_Sym2Vec(inf.norm2$P)))
summary(diff(F_Sym2Vec(inf.norm3$P)-F_Sym2Vec(inf.norm2$P)))
#############################
data.frame(inf.tree$L )%>% 
  mutate(index=1:length(inf.tree$L)) %>% 
  ggplot(aes(index,inf.tree.L))+
  geom_point(color="deepskyblue")+
  geom_line() + 
  facet_zoom(y= inf.tree.L >488.6)+
  theme_bw()


T1<-Sys.time()
inf.treeDist<-TreeGGM(cov2cor(Z.treeDist),print=FALSE,step="FALSE",maxIter = 150,n=nrow(Y))
# itermax=183, 13 min // alpha=0.25, 6.52min, itermax=29, pas d'erreur.
T2<-Sys.time()
difftime(T2,T1)

track_convergence(inf.vide,"Vide")
track_convergence(inf.tree,"Tree")
track_convergence(inf.treeDist,"Tree and distances")

save(inf.tree,inf.vide, inf.treeDist,file="Res_inf_maxlikelihood.RData")
############################################################
#compare beta1 beta2 inf finale
vide.beta1<-readRDS("seuil1vide.rds")
vide.beta2<-readRDS("seuil2vide.rds")
vide.finale<-inf.vide$P
plot(F_Sym2Vec(vide.beta1),F_Sym2Vec(vide.finale))
# 4. Coudes selon les modèles
compar_coudes<-data.frame(vide=sort(F_Sym2Vec(inf.vide$probaCond)), tree=sort(F_Sym2Vec(inf.tree$probaCond)),
                          treeDist=sort(F_Sym2Vec(inf.treeDist$probaCond)))
compar_coudes %>% 
  mutate(index=1:nrow(.)) %>% 
  gather(model,Edgeproba,-index) %>% 
  ggplot(aes(x=index,y=Edgeproba,color=model))+
  geom_point(size=0.3)+
  coord_cartesian(xlim=c(3800,nrow(compar_coudes)),ylim=c(0,0.5)) # !! zoom !!

# 5. jitter des degrés selon les modèles
maxDegrees<-1*(colSums(inf.vide$probaCond) %in% tail(sort(colSums(inf.vide$probaCond))))+
  2*(colSums(inf.tree$probaCond) %in% tail(sort(colSums(inf.tree$probaCond))))+
  3*(colSums(inf.treeDist$probaCond) %in% tail(sort(colSums(inf.treeDist$probaCond))))

compar_degres<-data.frame(vide=colSums(inf.vide$probaCond), tree=colSums(inf.tree$probaCond),
                          treeDist=colSums(inf.treeDist$probaCond),color=rep(c(1,2,1,3,1),c(13,1,29,1,50)),
                          color2= maxDegrees)
# color according to EA/F19/else
compar_degres %>% 
  mutate(color=as.factor(color),color2=as.factor(color2),index=1:nrow(compar_degres)) %>% 
  gather(model,degree,-color,-color2,-index) %>% 
  ggplot(aes(x=reorder(model,degree),y=degree,color=color, size=color))+
  scale_color_manual(values=c("lightslategrey","royalblue1","chartreuse3"),breaks=c(2,3,1),
                     labels=c("F19","EA","else"))+
  scale_size_manual(values=c(0.5,2.5,2.5),breaks=c(2,3,1),
                    labels=c("F19","EA","else"))+
  geom_jitter(width=0.2)+
  labs(x="",y="nodes degree estimates")

# color according to max degree
compar_degres %>% 
  mutate(color=as.factor(color),color2=as.factor(color2),index=1:nrow(compar_degres)) %>% 
  gather(model,degree,-color,-color2,-index) %>% 
  ggplot(aes(x=reorder(model,degree),y=degree,color=color2, size=color2))+
  scale_color_manual(name="max in:",values=c("lightslategrey","hotpink2","deepskyblue2","orange2"),breaks=c(6,5,1,0),
                     labels=c("all models","models with \ncovariates","empty","else"))+
  scale_size_manual(name="max in:",values=c(0.5,2.5,2.5,2.5),breaks=c(6,5,1,0),
                    labels=c("all models","models with \ncovariates","empty","else"))+
  geom_jitter(width=0.2)+
  labs(x="",y="nodes degree estimates")

# 6. ridges
compar_degres %>% 
  mutate(color=as.factor(color),index=1:nrow(compar_degres)) %>% 
  gather(model,degree,-color,-index) %>% 
  ggplot(aes(y=reorder(model,degree),x=degree,fill=factor(..quantile..))) +
  stat_density_ridges(scale=1.5,geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      quantiles = c(0.025, 0.975)) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
  )





############################################################
T1<-Sys.time()
inf.allcovariates<-TreeGGM(cov2cor(Z.allcovariates),print=FALSE,
                           step="FALSE",maxIter = 150)# 8 min avec eps=1e-6, maxier=79
T2<-Sys.time()
difftime(T2,T1)

saveRDS(inf.vide,"inf_vide.rds")
saveRDS(inf.tree,"inf_tree.rds")
saveRDS(inf.treeDist,"inf_treeDist.rds")
saveRDS(inf.allcovariates,"inf_allcovariates.rds")
inf.vide<-readRDS("inf_vide.rds")
inf.tree<-readRDS("inf_tree.rds")
inf.allcovariates<-readRDS("inf_allcovariates.rds")
scores<-data.frame(vide=c(inf.vide$P[upper.tri(inf.vide$P,diag = FALSE)]),
                   tree=c(inf.tree$P[upper.tri(inf.tree$P,diag = FALSE)]),
                   treeDist=c(inf.treeDist$P[upper.tri(inf.treeDist$P,diag = FALSE)]),
                   allcov=c(inf.allcovariates$P[upper.tri(inf.allcovariates$P,diag = FALSE)]),
                   index=1:length(c(inf.vide$P[upper.tri(inf.vide$P,diag = FALSE)])))
scores<-scores %>% mutate(color=vide) %>% 
  gather(model,proba,-"index",-"color") 
ggplot(scores,aes(model,proba,group=index,color=color))+
  geom_line()+ 
  scale_x_discrete(limits=c("vide","tree","treeDist","allcov"))

#################################################################################
#################################################################################
# Stability selection
# F_ResampleTreePLN <- function(Y, X, O, v = 0.8, B = 1e2){
#   n = nrow(Y)
#   p = ncol(Y)
#   P = p * (p - 1) / 2
#   V = round(v * n)
#   # Pmat = matrix(0, B, P)
#   
#   # browser()
#   # if(is.null(vecVar)){
#   #    fmla <- as.formula(paste("Y.sample ~ -1 + offset(log(O.sample))")) 
#   #  }else{
#   #    fmla <- as.formula(paste("Y.sample ~ -1 + offset(log(O.sample))+", 
#   #                             paste(paste0("X.sample$",vecVar), collapse= "+")))
#   #  }
#   res<-mclapply(1:B,function(x){
#     sample = sample(1:n, V, replace = F)
#     Y.sample = Y[sample,]
#     X.sample = X[sample,]
#     O.sample = O[sample,]
#     # PLN.sample = switch(model,"1"=PLN(Y.sample ~ 1 + offset(log(O.sample))),
#     #                     "2"=PLN(Y.sample ~ 1 + X.sample$tree+  offset(log(O.sample))),
#     #                     "3"=PLN(Y.sample ~ 1 + X.sample$tree+ X.sample$orientation+  offset(log(O.sample))),
#     #                     "4"=PLN(Y.sample ~ 1 + X.sample$tree+ X.sample$orientation+ X.sample$distTObase+  offset(log(O.sample))))
#     try({PLN.sample=PLN(Y.sample ~ 1 + X$tree+ offset(log(O.sample)))},silent=TRUE)
#     if(exists("PLN.sample")){
#       Sigma.sample = PLN.sample$model_par$Sigma
#     inf<-TreeGGM(cov2cor(Sigma.sample), "FALSE", FALSE, maxIter = 150)
#     res<-list(F_Sym2Vec(inf$P), F_Sym2Vec(inf$probaCond),inf$L)
#     }else{
#       res<-list(NA,NA,NA)
#     }
#     
#     return(res)
#   }, mc.cores=3)
#   return(res)
# }
if(Tree.res){
  Stab.sel = list()
  covar = list(); 
  covar[[1]] = matrix(1, n, 1); 
  covar[[2]] = as.matrix(lm(Y ~ X$tree, x=T)$x)
  covar[[3]] = as.matrix(lm(Y ~ X$tree+ X$distTObase+ X$distTOtrunk+X$distTOground, x=T)$x)
  T1<-Sys.time()
  sapply(1:3, function(m){
    Stab.sel[[m]] <<- F_ResampleTreePLN(Y, covar[[m]],log(O), B=B.resample, maxIter=3,
                                        cond.tol=1e-12)
  })
  T2<-Sys.time()
  difftime(T2,T1)
  save(Stab.sel, file = paste0(data.dir, data.name, '-StabSel.Rdata'))
}

load(paste0(data.dir, data.name, '-StabSel.Rdata'))
Pmats<-lapply(Stab.sel, function(x) {x[["Pmat"]]})

build_freqoak<-function(list){
  df<-data.frame(freq=double(), model=character(), index_edge=integer())
  names(list)<-c('~ 1', '~ tree', '~ tree + D1 +D2 + D3')
  lapply(seq_along(list), function(i){
    model<-names(list)[i]
    mat<-list[[i]][["Pmat"]] # mat B*nb(edges)
    freq<-colMeans(1 * (mat> 2 / p)) # freq nb(edges)
    df<<-rbind(df,data.frame(freq=freq, model=model, index_edge=seq_along(freq)))
  }) 
  return(df)
}
Fatfreq<-as_tibble(build_freqoak(Stab.sel))
library(ggraph)
library(tidygraph)
M<-4
p<-ncol(Y)
freq.sel.thres<-0.8
mat<-data.frame(F_Vec2Sym( 1*(colMeans( 1*(Stab.sel[[1]]$Pmat>2/p))>freq.sel.thres)))

allNets<-tibble(P = list(mat), seuil =c("1","2","3") )  %>% 
  mutate(P=map( seq_along(P), function(x) {
    df<-F_Vec2Sym( 1*(colMeans( 1*(Stab.sel[[x]]$Pmat>2/p))>freq.sel.thres))
    df[lower.tri(df, diag = TRUE)]<-0
    df<-data.frame(df)
    colnames(df)<-1:ncol(df)
    df
  })) %>% 
  mutate(P = map(P,~rownames_to_column(.) %>% 
                   gather(key, value , -rowname))) %>% 
  unnest() 
allNets<-allNets[,c(3,2,1,4)]
groupes<-as.factor(substr(colnames(Y),1,1))
set_graph_style()
mods<-c('~ 1', '~ tree', '~ tree + D1 +D2 + D3')
spliT<-data.frame(allNets) %>% 
  split(allNets$seuil) %>% 
  tibble(P=map(.,function(x){
    model<-switch(x$seuil[1],"1"="~ 1","2"="~ tree","3"="~ tree + D1 +D2 + D3")
#  browser()
    res<- as_tbl_graph(x, directed=FALSE) %>%
      activate(edges) %>%
      filter(value !=0) %>% 
      activate(nodes) %>% 
      mutate( groupe=groupes, importance=centrality_degree(),lab=colnames(Y),
              keyplayer = node_is_keyplayer(k=6), model=model,EA=(name%in%c("44"))) 
    res %>% activate(edges) %>% 
      mutate(neibs=edge_is_incident(which(.N()$EA)), model=model) %>% 
      activate(nodes) %>% 
      mutate(neibEA=name%in%unique(c(.E()$from[which(.E()$neibs)],.E()$to[which(.E()$neibs)]))) %>% 
      mutate(label=ifelse(keyplayer|EA|neibEA,lab,""))# %>% 
      #filter(importance!=0)
  }))
pal <- viridisLite::viridis(5, option = "C")
pal2<-c("gray15","#9de1e1")
coords<-create_layout(spliT$P[3][[1]],layout="star",center=44)[,c("x","y")]
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





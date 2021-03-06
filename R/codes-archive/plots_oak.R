library(RColorBrewer)
library(ggplot2)
###########################
# DATA
data.dir = '/home/momal/Git/these/Data/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))
theme_set(theme_bw())
Y = as.matrix(Data$count);
O = Data$offset; X = Data$covariates
YO = Y/O
Rank = rank(colSums(YO))
Seuil = 20
Ycum = colSums(Y); Order = order(Ycum)
Y = Y[, Rank > Seuil];
O = O[, Rank >Seuil];
n = nrow(Y); p = ncol(Y)

###########################
# MODELS & INFERENCES
PLN.vide<-PLN(Y ~ 1 + offset(log(O)))
PLN.tree<-PLN(Y ~ 1 + X$tree+  offset(log(O)))
PLN.treeDist<-PLN(Y ~ 1 +
                    X$tree+
                    X$distTObase+X$distTOtrunk+X$distTOground+
                    offset(log(O)))
inf.vide<-readRDS("inf_vide.rds")
inf.tree<-readRDS("inf_tree.rds")
inf.treeDist<-readRDS("inf_treeDist.rds")

###########################
# PLOTS
# 1. trajectoires vraisemblances EM
loglik<-data.frame(logL=c(inf.vide$L,inf.tree$L,inf.treeDist$L),
                   model=rep(c("vide","tree","treeDist"),c(length(inf.vide$L),
                                                           length(inf.tree$L),length(inf.treeDist$L))),
                   index=c(1:length(inf.vide$L),
                           1:length(inf.tree$L),1:length(inf.treeDist$L)))

ggplot(loglik,aes(x=index,y=logL,color=model))+
  geom_point()+
  geom_line(aes(group=model))


# 2. coudes entre proba  marg et cond
inf=inf.treeDist
coudes<-data.frame(marg=sort(F_Sym2Vec(inf$P)), cond=sort(F_Sym2Vec(inf$probaCond)))
coudes %>% 
  mutate(index=1:nrow(.)) %>% 
  gather(group,value,-index) %>% 
  ggplot(aes(x=index,y=value,color=group))+
  geom_point(size=0.3)

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
  ggplot(aes(y=reorder(model,degree),x=degree,fill = ..x..)) +
  geom_density_ridges_gradient(scale = 1.2, rel_min_height = 0.0001) +
  scale_fill_viridis(name = "Temp. [F]", option = "C")+
  labs(y="",x="nodes degree estimates")
  
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

# 7. degree variation

compar_degres %>% 
  mutate(btwvideandtree=vide-tree,btwtreeandtreeDist=tree-treeDist,
         btwvideandtreeDist=vide-treeDist,index=1:nrow(compar_degres)) %>% 
  gather(model,degree,-color,-color2,-index,-vide,-tree,-treeDist) %>% 
  ggplot(aes(y=reorder(model,degree),x=degree,fill=factor(..quantile..))) +
  stat_density_ridges(scale=3,geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      quantiles = c(0.025, 0.975)) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
  )
  
# 8. Spaghetti plot for fun
compar_coudes<-data.frame(vide=sort(F_Sym2Vec(inf.vide$P)), tree=sort(F_Sym2Vec(inf.tree$P)),
                          treeDist=sort(F_Sym2Vec(inf.treeDist$P)))
compar_coudes<-data.frame(vide=sort(F_Sym2Vec(inf.vide$probaCond)), tree=sort(F_Sym2Vec(inf.tree$probaCond)),
                          treeDist=sort(F_Sym2Vec(inf.treeDist$probaCond)))

compar_coudes %>% 
  mutate(index=1:nrow(.),color=vide) %>% 
  gather(model,proba,-index,-color) %>% 
  ggplot(aes(model,proba,group=index,color=color))+
  geom_line()+ 
  scale_x_discrete(limits=c("vide","tree","treeDist"))

  
  

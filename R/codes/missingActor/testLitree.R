library(saturnin)
library(reshape2)
library(patchwork)
#############################
# données avec structure de corrélation par bloc

p=40
k=3
dens=15/p
r=50

gr<-SimCluster2(p,k,dens,r)
omega=Laplacian(gr$G)+diag(0.01,ncol(gr$G))
sigma=solve(omega)
ordre=order(gr$groupes)
p1<-omega[ordre,ordre] %>% melt() %>% ggplot(aes(Var1,Var2,fill=value))+geom_tile()+
  labs(title="Omega",x="",y="")
p2<-sigma[ordre,ordre] %>% melt() %>% ggplot(aes(Var1,Var2,fill=value))+geom_tile()+
  labs(title="Sigma",x="",y="")
p3<-draw_network(gr$G,pal="black", layout="fr", groupes=gr$groupes, curv=0.1)$G+
  labs(title=paste0("p=",p,", density=",dens))
p3/(p1+p2)+plot_layout(ncol=1)

#############################
# test initEM avec structure de bloc qui vient du fait qu'on observe pas une variable

cliques=findCliques(sigma,3)
#verif trouve les groupes

groupestrouves=lapply(seq_along(cliques),function(x){
 L= length(cliques[[x]])
 rep(x,L)
})

found=data.frame(noeuds=unlist(cliques),find=unlist(groupestrouves))
init=data.frame(initgroupes=gr$groupes,noeuds=1:p)

found=left_join(found,init,by="noeuds")

table(found$find, found$initgroupes)

p1<-sigma[found$noeuds,found$noeuds] %>% melt() %>% ggplot(aes(Var1,Var2,fill=value))+geom_tile()+
  labs(title="cliques found",x="",y="")
p2<-sigma[ordre,ordre] %>% melt() %>% ggplot(aes(Var1,Var2,fill=value))+geom_tile()+
  labs(title="original cliques",x="",y="")
p1+p2+plot_layout(ncol=2)

#on est heureux.
#############################
# initEM

library(MASS)
initial.param<-initEM(sigma,n=500,cliquelist = cliques,pca=TRUE) # ajout de trois variables manquantes

# créer des nouvelles var très corrélées avec les cliques
# on retrouve que les cliques sont les variables les plus corrélées
# avec les nouvelles variables
sigma0=initial.param$Sigma0
K0=initial.param$K0
data.frame(cov2cor(sigma0)[41:43,-c(41,42,43)],missingvar=41:43) %>% 
  gather(key=noeuds,value,-missingvar) %>% 
  mutate(noeuds=as.numeric(substr(noeuds,2,3))) %>% 
  left_join(.,found,by="noeuds") %>% 
  ggplot(aes(y=abs(value),x=as.factor(find),color=as.factor(missingvar)))+
  geom_boxplot()+theme_light()+
  labs(x="Cliques",y="Absolute correlation")+
  scale_color_brewer("covariate",palette="Dark2")

#Chacune des variables ajoutées est préférentiellement corrélée avec une des cliques, 
# mais le niveau moyen de corrélation diffère
#globalement content. Peut être moyenne aussi bien que acp à tester


# ordre=rank(cov2cor(sigma0)[41,-c(41,42,43)])
# which(ordre%in%sort(ordre,decreasing=TRUE)[1:length(cliques$`1`)])
# cliques$`1`


#############################
# treeAgr
res<-treeAgr.EM(S = sigma,k=3, K0 = initial.param$K0,
                Sigma0 = initial.param$Sigma0,
                pii=0.5, n=500, max.iter = 100)
data.frame(res$alpha[41:43,-c(41,42,43)],missingvar=41:43) %>% 
  gather(key=noeuds,value,-missingvar) %>% 
  mutate(noeuds=as.numeric(substr(noeuds,2,3))) %>% 
  left_join(.,found,by="noeuds") %>% 
  ggplot(aes(y=abs(value),x=as.factor(find),color=as.factor(missingvar)))+
  geom_boxplot()+theme_light()+
  labs(x="Cliques")+
  scale_color_brewer("covariate",palette="Dark2")

#

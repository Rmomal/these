
#Sim15=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15.rds")
# object de taille 6*200 , chaque ligne est une liste de taille 200 contenant le paramètre en question
library(ggbeeswarm)
library(ggridges)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")
Sim15=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15_r1_200SF.rds")
p=14 ; r=1
H=p+r
#-------------------
# general processing
#times
timeboots<-mean(do.call(rbind,lapply(Sim15, function(seed){seed$time_boots})))
time1<-mean(do.call(rbind, lapply(Sim15,function(seed){seed$VEM_1$time})))

# links of missing actor
linksH<-data.frame(do.call(rbind,lapply(Sim15, function(seed){seed$omega[H,-H]})))
cpH<-data.frame(do.call(rbind,lapply(Sim15, function(seed){vecpartcor=-cov2cor(seed$omega)[H,-H]})))
colnames(cpH)<-1:ncol(cpH)
voisvois<-do.call(rbind,lapply(Sim15, function(seed){
  vois=1*(seed$omega[H,-H]!=0) # vect size 14 1 vois 0 pas vois
  di=round(diag(seed$omega)[-H],0) # degré des 14 premiers noeuds
  di[vois==0]=0 # mise à zéro des degré des non voisins de H
  return(di)
}))

voisMaxDeg<-do.call(rbind,lapply(Sim15, function(seed){
  #faudrait mise à zero de ceux qui ne sont pas voisin d'un voisin de H
  index.vois=which(seed$omega[H,-H]!=0)
  vecmax<-rep(0,H-1)
  sapply(index.vois, function(vois){
    if(round(diag(seed$omega)[vois],0)>1){
      vv=setdiff(which(seed$omega[vois,-H]!=0),vois) # vecteur des noeuds voisins de vois à part H
      maxDeg=round(max(diag(seed$omega)[vv]),0)
      vecmax[vois]<<-maxDeg
    }
  })
  return(vecmax)
}))
voisBetween<-do.call(rbind,lapply(Sim15, function(seed){
  between=draw_network(seed$omega)$graph_data %>% activate(nodes) %>% as_tibble() %>% select(btw)
  between=between$btw/max(between$btw)
  between[seed$omega[H,-H]==0] = 0
  between=between[-H]
  return(between)
}))

Hbtw<-do.call(rbind,lapply(Sim15, function(seed){
  between=draw_network(seed$omega)$graph_data %>% activate(nodes) %>% as_tibble() %>% select(btw)
  betweenH=between$btw[H]
  return(betweenH)
}))

misAct<-data.frame(nH=rowSums(1*(cpH!=0)), 
                   mean.cpH=apply(cpH, 1,function(vec){mean(vec[vec!=0])}),
                   mean.vv=apply(voisvois, 1,function(vec){mean(vec[vec!=0])}),
                   max.maxDeg=apply(voisMaxDeg, 1,function(vec){max(vec[vec!=0])}), 
                   mean.maxDeg=apply(voisMaxDeg, 1,function(vec){mean(vec[vec!=0])}), 
                   mean.btw=apply(voisBetween, 1,function(vec){mean(vec[vec!=0])}),
                   btwH=Hbtw)


# omega conditioning
cond.omgraph<-do.call(rbind, lapply(Sim15, function(seed){
  lambda=svd(seed$omega)$d
  cond=min(abs(lambda)/max(abs(lambda)))
  return(cond)}))
cond.siginf<-do.call(rbind, lapply(Sim15, function(seed){
  M=seed$VEM_1$M
  S=seed$VEM_1$S
  n=nrow(M)
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  lambda=svd(solve(SigmaTilde))$d
  cond=min(abs(lambda)/max(abs(lambda)))
  return(cond)}))
cond.sigobs<-do.call(rbind, lapply(Sim15, function(seed){
  M=seed$VEM_1$M[,-H]
  S=seed$VEM_1$S[,-H]; n=nrow(M)
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  lambda=svd(solve(SigmaTilde))$d
  cond=min(abs(lambda)/max(abs(lambda)))
  return(cond)}))

cond.ominf<-do.call(rbind, lapply(Sim15, function(seed){
  omega=seed$VEM_1$omega
  lambda=svd(solve(omega))$d
  cond=min(abs(lambda)/max(abs(lambda)))
  return(cond)}))


# inference results
PgH<-data.frame(Pg=do.call(rbind,lapply(Sim15, function(seed){seed$VEM_1$Pg[H,-H]})))
colnames(PgH)<-1:ncol(PgH)
qual<-data.frame(do.call(rbind, lapply(Sim15, function(seed){
  Pg=seed$VEM_1$Pg
  ome=seed$omega
  diag(ome)=0
  auc=round(auc(pred = Pg, label = ome),4) 
  vec=accppvtpr(Pg,ome,h=H,seuil=0.5)
  return(c(auc=auc, vec))
}))) %>% as_tibble() %>% mutate(seed=1:200)
colnames(qual)=c("auc","Acc","AccH","AccO","PPV","PPVH","PPVO","TPR","TPRH","TPRO","seed")


#ZH MH
rebuild<-do.call(rbind,lapply(Sim15, function(seed){
  ZH<-seed$ZH
  MH<-seed$VEM_1$M[,H]
  abs(cor(ZH, MH))
}))
qual$rebuild=rebuild[,1]

# Convergence

convergence<-data.frame(iter=do.call(rbind,lapply(Sim15,function(seed){
  iter=seed$VEM_1$finalIter })), 
  nbconv=do.call(rbind,lapply(Sim15,function(seed){
    nbconv=seed$nbconv })))


################
# Detection
#------------
# l'acteur manquant est détecté 100% du temps
# vBICs %>% ggplot(aes(vBIC0,vBIC1))+geom_point()+mytheme+
#   geom_abline()+coord_cartesian(ylim=c(-17000,-7000))+
#   labs(title="Missing actor was always detected")
# 


################
# Quality
#------------
# Graph total

qual %>% gather(key, value,-seed, -rebuild,-auc)  %>%
  mutate(key = fct_reorder(key, value, .fun='median')) %>%
  mutate(type=substr(key,1,3)) %>% 
  ggplot(aes( value,key, color=key, fill=key))+
  geom_density_ridges(alpha=0.3)+theme_light()+guides(color=FALSE, fill=FALSE)+
  facet_wrap(~type)


# PPVH a deux modes : pourquoi ?
data=cbind(qual, misAct,convergence, cond.omgraph,cond.ominf, cond.siginf)  %>%
  mutate(bad=as.factor(ifelse(PPVH<=0.5, "<0.5", ">0.5"))) %>% as_tibble()

data %>%as_tibble() %>%  filter(!is.na(bad))   %>% 
  dplyr::select(mean.cpH,mean.vv,nH,iter,nbconv,bad) %>% 
  gather(key, value,-bad) %>% 
  ggplot(aes(bad,value, color=key, fill=key))+
  geom_boxplot(alpha=0.3, width=0.5)+theme_light()+facet_wrap(~key, scales = "free", nrow=2)+
  guides(color=FALSE, fill=FALSE)+labs(x="PPVH")

data %>%as_tibble() %>%  filter(!is.na(bad))   %>% 
  dplyr::select(max.maxDeg, mean.btw,btwH,nH,bad) %>% 
  gather(key, value,-bad) %>% 
  ggplot(aes(bad,value, color=key, fill=key))+
  geom_boxplot(alpha=0.3, width=0.5)+mytheme.dark("")+facet_wrap(~key, scales = "free", nrow=1)+
  guides(color=FALSE, fill=FALSE)+labs(x="PPVH")

data %>%as_tibble() %>%  filter(!is.na(bad))   %>% 
  dplyr::select(mean.maxDeg,nH,bad) %>%  
  ggplot(aes(bad,mean.maxDeg, color=as.factor(nH), fill=as.factor(nH)))+
  geom_boxplot(alpha=0.3, width=0.5)+theme_light()+facet_wrap(~nH, nrow=1)+
  guides(color=FALSE, fill=FALSE)+labs(x="PPVH")
# ça c'est intéressant
data %>%as_tibble() %>%  filter(!is.na(bad))   %>% 
  dplyr::select(btwH,mean.btw,nH,bad) %>%  
  ggplot(aes(bad,mean.btw, color=as.factor(nH), fill=as.factor(nH)))+
  geom_boxplot(alpha=0.3, width=0.5)+theme_light()+facet_wrap(~nH, nrow=1)+
  guides(color=FALSE, fill=FALSE)+labs(x="PPVH")


data %>%as_tibble() %>%  
  dplyr::select(nH,btwH) %>%
  ggplot(aes(nH,btwH/nH))+geom_point()+mytheme.dark("") 

data %>% ggplot(aes(TPRH, PPVH,color=(mean.maxDeg)))+geom_point()+theme_light()
data %>% ggplot(aes(TPRH, (nH)))+geom_point()+theme_light()
data %>% ggplot(aes(auc, mean.cpH, color=bad))+geom_point()+mytheme.dark("")
data %>% ggplot(aes(auc, mean.cpH, color=btwH))+geom_point()+theme_light()
data %>% filter(!is.na(bad)) %>%  ggplot(aes(as.factor(nH), mean.cpH, color=bad))+
  geom_point()+mytheme.dark("")+facet_wrap(~bad)


data %>%  
  ggplot(aes(as.factor(nH), PPVH, color=as.factor(nH), fill=as.factor(nH)))+
  geom_boxplot(alpha=0.3)+ theme_light()+guides(color=FALSE, fill=FALSE)

#effet du conditionnement du graph de départ sur PPVH est un biais du nombre de nH
# data %>% as_tibble() %>% 
#   ggplot(aes(cond.omgraph, PPVH, color=as.factor(nH), fill=as.factor(nH)))+
#   geom_point(alpha=0.3)+ theme_light()+guides(color=FALSE, fill=FALSE)
# 
# data %>% as_tibble() %>% 
#   ggplot(aes(bad,cond.omgraph))+
#   geom_boxplot(alpha=0.3)+ theme_light()+guides(color=FALSE, fill=FALSE)+facet_wrap(~nH)

# effet de nH sur PPVH
data %>% as_tibble() %>% 
  ggplot(aes(nH, PPVH, color=as.factor(nH), fill=as.factor(nH)))+
  geom_point(alpha=0.3)+ theme_light()+guides(color=FALSE, fill=FALSE)

data %>% as_tibble() %>% 
  ggplot(aes(bad,nH  ))+
  geom_boxplot(alpha=0.3)+ theme_light()+guides(color=FALSE, fill=FALSE)
# effet mean.vv sur PPVH
data %>% 
  ggplot(aes(bad,mean.vv))+geom_boxplot()+mytheme.dark("")+facet_wrap(~nH)

# valeurs de mean.cpH difficilement interprétable
data %>% as_tibble() %>% 
  ggplot(aes(mean.cpH, PPVH, color=as.factor(nH), fill=as.factor(nH)))+
  geom_point(alpha=0.3)+ theme_light()+guides(color=FALSE, fill=FALSE)
data %>% as_tibble() %>% 
  ggplot(aes(bad, abs(mean.cpH), color=as.factor(nH), fill=as.factor(nH)))+
  geom_boxplot(alpha=0.3)+ theme_light()+guides(color=FALSE, fill=FALSE)

# PPVH et rebuild
data %>% as_tibble() %>% 
  ggplot(aes( PPVH,rebuild, color=as.factor(nH), fill=as.factor(nH)))+
  geom_abline(linetype="dashed", color="gray")+
  geom_hline(yintercept = 0.6, linetype="dashed",color="gray")+
  geom_point(alpha=0.4)+ theme_light()+guides(color=FALSE, fill=FALSE)+facet_wrap(~nH)

# PPVH et mean.btw
data %>% as_tibble() %>% 
  ggplot(aes( PPVH,mean.btw, color=as.factor(nH), fill=as.factor(nH)))+
  geom_point(alpha=0.3)+ theme_light()+guides(color=FALSE, fill=FALSE)+facet_wrap(~nH)

# iter, une vraie info
data %>% as_tibble() %>% 
  ggplot(aes( PPVH,iter, color=as.factor(nH), fill=as.factor(nH)))+
  geom_point(alpha=0.3)+ theme_light()+guides(color=FALSE, fill=FALSE)+facet_wrap(~nH)

# certains auc <0.5, pourquoi
data=cbind(qual, misAct,convergence, cond.omgraph,cond.ominf, cond.siginf)  %>%
  mutate(bad=as.factor(ifelse(auc<=0.5, "auc<0.5", "auc>0.5")))

#------------
# par arêtes
tmp1<-t(PgH)%>% as_tibble()%>% gather(key, prob) %>% mutate(seed=as.numeric(as.factor(key))) 
tmp2<-t(cpH) %>% as_tibble()%>% gather(key2, cpH)%>% mutate(seed2=as.numeric(as.factor(key2)))
tmp3<-t(linksH) %>% as_tibble() %>% gather(key3, neigb)%>% mutate(seed3=as.numeric(as.factor(key3)))
tmp4<-t(voisvois) %>% as_tibble() %>% gather(key4, vv)%>% mutate(seed4=as.numeric(as.factor(key4)))

probInf<-cbind(tmp1, tmp2, tmp3, tmp4) %>% dplyr::select(seed, prob, cpH,vv, neigb)%>% 
  as_tibble() 
datanH<-data.frame(nH=misAct$nH,seed=1:200)
join<-left_join(probInf,datanH, by="seed")
join %>% 
  ggplot(aes(prob,as.factor(neigb), color=as.factor(neigb), fill=as.factor(neigb)))+
  geom_density_ridges(alpha=0.3)+mytheme.dark("")+guides(color=FALSE, fill=FALSE)+
  facet_wrap(~nH)
join %>% 
  ggplot(aes(prob,as.factor(vv), color=as.factor(vv), fill=as.factor(vv)))+
  geom_density_ridges(alpha=0.3)+mytheme.dark("")+guides(color=FALSE, fill=FALSE)+
  facet_wrap(~nH)


# beaucoup de voisins non trouvés
join %>% filter(neigb==1) %>% 
  ggplot(aes(prob,(nH),  color=as.factor(vv)))+geom_point()+theme_light()

join %>% filter(neigb==1) %>% 
  ggplot(aes(vv,prob))+geom_point()+theme_light()

join %>% filter(neigb==1) %>% 
  ggplot(aes(prob, abs(cpH), color=as.factor(vv)))+geom_point()+theme_light()


#################@
# etude d'un cas de non convergence
converged %>% group_by(converged) %>% summarise(nb=n())
# 7 sur 200 n'ont pas convergé malgré les 40 essais d'initialisation
badseeds=which(converged$converged==0)
do.call(cbind,lapply(Sim15[badseeds], function(seed){
  seed$VEM_1$feature$diffW
}))
do.call(grid.arrange,lapply(Sim15[badseeds[1]], function(seed){
  seed$VEM_1$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>% 
    filter(key=="J") %>% 
    ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(parameter)), size=3)+geom_line()+
    facet_wrap(~key, scales="free")+
    labs(x="iteration",y="", title="Lower bound")+mytheme+
    scale_color_discrete("")+coord_cartesian(xlim=c(0,30))
}))


do.call(grid.arrange,lapply(Sim15[badseeds], function(seed){
  seed$VEM_1$features %>% rowid_to_column() %>%  gather(key,value,-rowid) %>% 
    filter(key=="diffW") %>% 
    ggplot(aes(rowid,value, group=key))+geom_point()+geom_line()+
    facet_wrap(~key, scales="free")+
    labs(x="iteration",y="", title="")+mytheme+
    scale_color_discrete("")+coord_cartesian(xlim=c(99,100))
}))

#################
#vBIC 0 quand pas de manquants

Sim15_r0=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15_r0.rds")

vBICs<-data.frame(vBIC0=do.call(rbind,lapply(Sim15_r0,function(seed){seed$vBIC_0})),
                  vBIC1=do.call(rbind,lapply(Sim15_r0,function(seed){seed$vBIC_1})))%>%
  as_tibble() %>% 
  mutate(best=ifelse(vBIC1>vBIC0,"vBIC1","vBIC0"))
vBICs %>% ggplot(aes(vBIC0,vBIC1))+geom_point()+mytheme+
  geom_abline()#+  coord_cartesian(ylim=c(-17000,-7000))+
labs(title="Missing actor was always detected")

alphas<-data.frame(alpha1=do.call(rbind,  lapply(Sim15_r0,function(seed){seed$VEM_1$alpha})),
                   alpha0=do.call(rbind,  lapply(Sim15_r0,function(seed){seed$VEM_0$alpha})))
alphas %>%  
  ggplot(aes(alpha0,alpha1 ))+geom_point()+theme_light()+
  geom_abline()

vem=Sim15_r0[[1]]
True_lowBound<-function(Y, M,S,theta,X, W, Wg, Pg, omega)
  
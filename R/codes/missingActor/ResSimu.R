
#Sim15=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15.rds")
# object de taille 6*200 , chaque ligne est une liste de taille 200 contenant le paramètre en question
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")
Sim15=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15conv.rds")

#-------------------
# general processing
#times
timeboots<-mean(do.call(rbind,lapply(Sim15, function(seed){seed$time_boots})))
time1<-mean(do.call(rbind, lapply(Sim15,function(seed){seed$VEM_1$time})))
time0<-mean(do.call(rbind,  lapply(Sim15,function(seed){seed$VEM_0$time})))
#alphas
alphas<-data.frame(alpha1=do.call(rbind,  lapply(Sim15,function(seed){seed$VEM_1$alpha})),
                   alpha0=do.call(rbind,  lapply(Sim15,function(seed){seed$VEM_0$alpha})))
alphas %>% mutate(mean.cpH=misAct$mean.cpH) %>% 
  ggplot(aes(alpha0,alpha1, color=mean.cpH))+geom_point()+theme_light()+
  geom_abline()
# links of missing actor
linksH<-data.frame(do.call(rbind,lapply(Sim15, function(seed){seed$omega[15,-15]})))
cpH<-data.frame(do.call(rbind,lapply(Sim15, function(seed){vecpartcor=-cov2cor(seed$omega)[15,-15]})))
colnames(cpH)<-1:ncol(cpH)
voisvois<-do.call(rbind,lapply(Sim15, function(seed){
  vois=1*(seed$omega[15,-15]!=0)
  di=round(diag(seed$omega)[-15],0)
  di[vois==0]=0
  return(di)
}))
misAct<-data.frame(nH=rowSums(1*(cpH!=0)), 
                   mean.cpH=apply(cpH, 1,function(vec){mean(vec[vec!=0])}),
                   mean.vv=apply(voisvois, 1,function(vec){mean(vec[vec!=0])}))


# omega conditioning
cond.omgraph<-do.call(rbind, lapply(Sim15, function(seed){
  lambda=svd(seed$omega)$d
  cond=min(abs(lambda)/max(abs(lambda)))
  return(cond)}))
cond.siginf<-do.call(rbind, lapply(Sim15, function(seed){
  M=seed$VEM_1$M
  S=seed$VEM_1$S
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  lambda=svd(solve(SigmaTilde))$d
  cond=min(abs(lambda)/max(abs(lambda)))
  return(cond)}))
cond.sigobs<-do.call(rbind, lapply(Sim15, function(seed){
  M=seed$VEM_1$M[,-15]
  S=seed$VEM_1$S[,-15]
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  lambda=svd(solve(SigmaTilde))$d
  cond=min(abs(lambda)/max(abs(lambda)))
  return(cond)}))

cond.ominf<-do.call(rbind, lapply(Sim15, function(seed){
  omega=seed$VEM_1$omega
  lambda=svd(solve(omega))$d
  cond=min(abs(lambda)/max(abs(lambda)))
  return(cond)}))

# vBIC
vBICs<-data.frame(vBIC0=do.call(rbind,lapply(Sim15,function(seed){seed$vBIC_0})),
                  vBIC1=do.call(rbind,lapply(Sim15,function(seed){seed$vBIC_1})))%>%
  as_tibble() %>% 
  mutate(best=ifelse(vBIC1>vBIC0,"vBIC1","vBIC0"))

# inference results
PgH<-data.frame(Pg=do.call(rbind,lapply(Sim15, function(seed){seed$VEM_1$Pg[15,-15]})))
colnames(PgH)<-1:ncol(PgH)
qual<-data.frame(do.call(rbind, lapply(Sim15, function(seed){
  Pg=seed$VEM_1$Pg
  ome=seed$omega
  diag(ome)=0
  auc=round(auc(pred = Pg, label = ome),4) 
  vec=accppvtpr(Pg,ome,h=15,seuil=0.5)
  return(c(auc=auc, vec))
}))) %>% as_tibble() %>% mutate(seed=1:200)
colnames(qual)=c("auc","Acc","AccH","AccO","PPV","PPVH","PPVO","TPR","TPRH","TPRO","seed")
################
# Detection
#------------
# l'acteur manquant est détecté 100% du temps
vBICs %>% ggplot(aes(vBIC0,vBIC1))+geom_point()+mytheme+
  geom_abline()+coord_cartesian(ylim=c(-17000,-7000))+
  labs(title="Missing actor was always detected")


vBICs %>%mutate(diff=vBIC1-vBIC0) %>% 
  ggplot(aes(y=log(diff), x=" "))+geom_beeswarm()+mytheme+
  labs(title="Missing actor was always detected", x="",y="log(vBIC1-vBIC0)")+
  coord_cartesian(y=c(5,log(8000)))

#------------
# Convergence

converged<-data.frame(diffW=do.call(rbind,lapply(Sim15,function(seed){
  diffW=seed$VEM_1$features$diffW
  return(diffW[length(diffW)])})), 
  diffOm=do.call(rbind,lapply(Sim15,function(seed){
    diffOm=seed$VEM_1$features$diffOm
    return(diffOm[length(diffOm)])})))

converged=converged %>% as_tibble() %>% 
  mutate(converged=ifelse(diffW<1e-3 & diffOm< 1e-3, 1, 0))
################
# Quality
#------------
# Graph total
qual %>% filter(auc>0.5)
# qual %>% gather(key, value,-seed) %>%
#   mutate(key = fct_reorder(key, value, .fun='median')) %>%
#   ggplot(aes(key, value, color=key))+geom_boxplot()+theme_light()+guides(color=FALSE)
qual %>% gather(key, value,-seed)  %>%
  mutate(key = fct_reorder(key, value, .fun='median')) %>%
  ggplot(aes( value,key, color=key, fill=key))+geom_density_ridges(alpha=0.3)+theme_light()+guides(color=FALSE, fill=FALSE)
# qual %>%
#   ggplot(aes( PPVH,TPRH))+geom_point()+theme_light()+guides(color=FALSE, fill=FALSE)

# PPVH a deux modes : pourquoi ?
data=cbind(qual, misAct, alphas,converged, cond.omgraph,cond.ominf,cond.sigobs, cond.siginf)  %>%
  mutate(bad=as.factor(ifelse(PPVH<=0.5, 1, 0)))

data %>% filter(!is.na(bad))  %>% mutate(ldiffW=log(diffW)) %>% 
  dplyr::select(mean.cpH,mean.vv,alpha1,nH,cond.omgraph, cond.ominf,cond.sigobs,
                cond.siginf,ldiffW,bad) %>% 
  gather(key, value,-bad) %>% 
  ggplot(aes(bad,value, color=key, fill=key))+
  geom_boxplot(alpha=0.3, width=0.5)+theme_light()+facet_wrap(~key, scales = "free", nrow=4)+
  guides(color=FALSE, fill=FALSE)+labs(x="PPVH < 0.5")

data %>% ggplot(aes(auc, mean.cpH, color=bad))+geom_point()+mytheme.dark
data %>% filter(!is.na(bad)) %>%  ggplot(aes(as.factor(nH), mean.cpH, color=bad))+
  geom_point()+mytheme.dark+facet_wrap(~bad)

# étudier les qualités en fonction de converged
# effet de la non convergence sur les performances
data %>% 
  dplyr::select(converged, auc, Acc, AccH, AccO, PPV, PPVH, PPVO, TPR, TPRO,TPRH) %>% 
  gather(key, value,-converged) %>% 
  ggplot(aes(as.factor(converged),value, color=key, fill=key))+
  geom_boxplot(alpha=0.3, width=0.5)+theme_light()+facet_wrap(~key, scales = "free", nrow=4)+
  guides(color=FALSE, fill=FALSE)+labs(x="converged")

# pourquoi la non-convergence
data  %>% 
  dplyr::select(mean.cpH,mean.vv,alpha1,nH,cond.omgraph,cond.sigobs,
                converged) %>% 
  gather(key, value,-converged) %>% 
  ggplot(aes(as.factor(converged),value, color=key, fill=key))+
  geom_boxplot(alpha=0.3, width=0.5)+theme_light()+facet_wrap(~key, scales = "free", nrow=4)+
  guides(color=FALSE, fill=FALSE)+labs(x="PPVH < 0.5")

qual %>%  
  ggplot(aes(PPVH, TPRO, color=TPRH))+geom_point()+theme_light()+guides(color=FALSE)

qual %>%  
  ggplot(aes((PPVH), auc, color=as.factor(PPVH)))+geom_boxplot()+theme_light()+guides(color=FALSE)

#------------
# par arêtes
tmp1<-t(PgH)%>% as_tibble()%>% gather(key, prob) %>% mutate(seed=as.numeric(as.factor(key))) 
tmp2<-t(cpH) %>% as_tibble()%>% gather(key2, cpH)%>% mutate(seed2=as.numeric(as.factor(key2)))
tmp3<-t(linksH) %>% as_tibble() %>% gather(key3, neigb)%>% mutate(seed3=as.numeric(as.factor(key3)))
tmp4<-t(voisvois) %>% as_tibble() %>% gather(key4, vv)%>% mutate(seed4=as.numeric(as.factor(key4)))

probInf<-cbind(tmp1, tmp2, tmp3, tmp4) %>% dplyr::select(seed, prob, cpH,vv, neigb)%>% 
  as_tibble() 
alpha2<-cbind(alphas, misAct) %>% as_tibble() %>%  mutate(seed=1:nrow(alphas))
join=left_join(probInf, alpha2, by="seed")
join %>% 
  ggplot(aes(prob,as.factor(neigb), color=as.factor(neigb), fill=as.factor(neigb)))+
  geom_density_ridges(alpha=0.3)+mytheme.dark+guides(color=FALSE, fill=FALSE)+
  facet_wrap(~nH)
join %>% 
  ggplot(aes(prob,as.factor(vv), color=as.factor(vv), fill=as.factor(vv)))+
  geom_density_ridges(alpha=0.3)+mytheme.dark+guides(color=FALSE, fill=FALSE)+
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
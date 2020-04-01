
#Sim15=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Sim15.rds")
# object de taille 6*200 , chaque ligne est une liste de taille 200 contenant le paramètre en question
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")
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
misAct<-data.frame(nH=rowSums(1*(cpH!=0)), 
                   mean.cpH=apply(cpH, 1,function(vec){mean(vec[vec!=0])}))

misAct %>% ggplot(aes(as.factor(nH), mean.cpH))+geom_beeswarm()+mytheme

# omega conditioning
cond<-do.call(rbind, lapply(Sim15, function(seed){
  lambda=svd(seed$omega)$d
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
#-------------------
# Detection
# l'acteur manquant est détecté 100% du temps
vBICs %>% ggplot(aes(vBIC0,vBIC1))+geom_point()+mytheme+
  geom_abline()+coord_cartesian(ylim=c(-17000,-7000))+
  labs(title="Missing actor was always detected")


vBICs %>%mutate(diff=vBIC1-vBIC0) %>% 
  ggplot(aes(y=log(diff), x=" "))+geom_beeswarm()+mytheme+
  labs(title="Missing actor was always detected", x="",y="log(vBIC1-vBIC0)")+
  coord_cartesian(y=c(5,log(8000)))

#-------------------
# Convergence

converged<-data.frame(diffW=do.call(rbind,lapply(Sim15,function(seed){seed$VEM_1$features$diffW})),
                   diffOm=do.call(rbind,lapply(Sim15,function(seed){seed$VEM_1$features$diffOm})))


#-------------------
# Quality
qual %>% filter(auc>0.5)
qual %>% gather(key, value,-seed) %>%
  mutate(key = fct_reorder(key, value, .fun='median')) %>%
  ggplot(aes(key, value, color=key))+geom_boxplot()+theme_light()+guides(color=FALSE)
qual %>% gather(key, value,-seed)  %>%
  mutate(key = fct_reorder(key, value, .fun='median')) %>%
  ggplot(aes( value,key, color=key, fill=key))+geom_density_ridges(alpha=0.3)+theme_light()+guides(color=FALSE, fill=FALSE)
qual %>%
  ggplot(aes( PPVH,TPRH))+geom_point()+theme_light()+guides(color=FALSE, fill=FALSE)

# PPVH a deux modes : pourquoi ?
data=cbind(qual, misAct, alphas, cond)  %>% mutate(bad=as.factor(ifelse(PPVH<=0.5, 1, 0)))
data %>% group_by(bad) %>% summarise(pct=n()*100/nrow(data))
data %>% filter(!is.na(bad)) %>% dplyr::select(mean.cpH,alpha1,alpha0,nH,cond, bad) %>% 
  gather(key, value,-bad) %>% 
  ggplot(aes(bad,value, color=key, fill=key))+
  geom_boxplot(alpha=0.3)+mytheme.dark+facet_wrap(~key, scales = "free")+
  guides(color=FALSE, fill=FALSE)+labs(x="PPVH < 0.5")

data %>% ggplot(aes(auc, mean.cpH, color=as.factor(bad)))+geom_point()+mytheme.dark
data %>% filter(!is.na(bad)) %>%  ggplot(aes(as.factor(nH), mean.cpH, color=bad))+
  geom_point()+mytheme.dark+facet_wrap(~bad)

qual %>%  
  ggplot(aes((PPVH), auc, color=as.factor(PPVH)))+geom_boxplot()+theme_light()+guides(color=FALSE)

qual %>%  
  ggplot(aes((PPVH), auc, color=as.factor(PPVH)))+geom_boxplot()+theme_light()+guides(color=FALSE)


tmp1<-t(PgH)%>% as_tibble()%>% gather(key, prob) %>% mutate(seed=as.numeric(as.factor(key))) 
tmp2<-t(cpH) %>% as_tibble()%>% gather(key2, cpH)%>% mutate(seed2=as.numeric(as.factor(key2)))
tmp3<-t(linksH) %>% as_tibble() %>% gather(key3, neigb)%>% mutate(seed3=as.numeric(as.factor(key3)))
probInf<-cbind(tmp1, tmp2, tmp3) %>% dplyr::select(seed, prob, cpH, neigb)%>% 
  as_tibble() 
alpha2<-cbind(alphas, misAct) %>% as_tibble() %>%  mutate(seed=1:nrow(alphas))
join=left_join(probInf, alpha2, by="seed")
join %>% 
  ggplot(aes(prob,as.factor(neigb), color=as.factor(neigb), fill=as.factor(neigb)))+
  geom_density_ridges(alpha=0.3)+mytheme.dark+guides(color=FALSE, fill=FALSE)+
  facet_wrap(~nH)
# beaucoup de voisins non trouvés
join %>% filter(neigb==1) %>% 
  ggplot(aes(as.factor(nH), prob, color=mean.cpH))+geom_point()+theme_light()
 
join %>% filter(neigb==1) %>% 
  ggplot(aes(prob, abs(cpH), color=alpha1))+geom_point()+theme_light()


cbind(qual, misAct, alphas) %>% 
  ggplot(aes(auc, mean.cpH, color=alphas_1))+geom_point()+mytheme

 
cbind(qual, misAct, alphas) %>% 
  ggplot(aes(PPVH, mean.cpH))+geom_point()+mytheme

cbind(qual, misAct, alphas) %>% 
  ggplot(aes((abs(mean.cpH)),(alpha1)))+geom_point()+mytheme
cbind(qual, misAct, alphas) %>% 
  ggplot(aes(PPVH,(alpha1)))+geom_point()+mytheme

cbind(qual, misAct) %>% group_by(PPVH) %>% summarise(mcpH=mean(mean.cpH)) %>% 
  ggplot(aes( PPVH, mcpH))+geom_point()+mytheme


cbind(qual, misAct) %>% 
  ggplot(aes(auc, mean.cpH, color=as.factor(auc)))+geom_boxplot()+mytheme+
  guides(color=FALSE)

cbind(qual, misAct) %>% 
  ggplot(aes(PPVH,(auc), color=nH))+geom_point()+mytheme

cbind(qual, misAct) %>% 
  ggplot(aes(as.factor(nH), PPVH))+geom_boxplot()+mytheme

cbind(qual, misAct, 1:200) %>%  filter(PPVH<0.3)

 
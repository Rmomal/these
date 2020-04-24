source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")
# data
simus=list()

simus=lapply(1:200, function(seed){
  readRDS(paste0(
    "/Users/raphaellemomal/these/R/codes/missingActor/SimResults/15nodes_1r/SF_seed",seed,".rds"))
  })
#--- Structure de l'objet simus:
#   simus est une list de 200 éléments contenant les résultats de simulation pour chacune des 200 graines.
#   simus[[graine]]$ omega || ZH || VEM_1 || time_boots || nbinit
#   simus[[graine]]$VEM_1$ M || S || Pg || Wg || W || omega || lowbound || 
#                         features || finalIter || time || clique || nbocc
#-------------------
#-- run time & cie
# toutes les simus ont foncitonné :
nb.bug=sum(do.call(rbind, lapply(simus, function(graine){
  length(graine$VEM_1)!=12
}))) #nb.bug=0

# nombre moyen de cliques initiales testées pour chaque graine :
nbinit=do.call(rbind, lapply(simus, function(graine){
  graine$nbinit
})) # moy.nbinit = 31.94, sd = 13.41215
boxplot(nbinit, horizontal=TRUE)

#temps moyen pour trouver les cliques d'initialisation par spca
time_boots=do.call(rbind, lapply(simus, function(graine){
  graine$time_boots
})) # temps faible : entre 4 et 14 secondes, mean= 6.5 sd = 1.2
boxplot(time_boots, horizontal=TRUE)

# run time de VEMtree
time_VEMtree=do.call(rbind, lapply(simus, function(graine){
  graine$VEM_1$time
})) # temps moyen : entre 1 et 60s, mean=20 sd = 11
boxplot(time_VEMtree, horizontal=TRUE)

#20*30*200/(60*60) : environ 33h pour 200 graines sans calcul parallèle, environ 13 avec.

#-------------------
# graphs générés
misAct=build_misAct(simus, H=15)

#-------------------
#-- VEM_1 convergence
# nb iterations
nb.iter=do.call(rbind, lapply(simus, function(graine){
  graine$VEM_1$finalIter
})) # longue queue : entre 12 et 180, mean=36 sd = 22
boxplot(nb.iter, horizontal=TRUE)
# maxIter = 200 donc tous les modèles ont convergés. Une seule valeure au-dessus de 150

# nb occurence de la meilleure clique trouvée :
nb.occ=do.call(rbind, lapply(simus, function(graine){
  graine$VEM_1$nbocc
})) # très longue queue : entre 1 et 40, med= 1, mean= 3.5, sd = 5
boxplot(nb.occ, horizontal=TRUE)

# nb de stops par J au lieu des paramètres : 
stopJ=do.call(rbind, lapply(simus, function(graine){
  diffW = tail(graine$VEM_1$features$diffW,1)
  diffOm = tail(graine$VEM_1$features$diffOm,1)
  return(diffW>1e-3 || diffOm>1e-3)
}))
nb.stopJ = sum(stopJ)# nb.stopJ = 145
# est-ce qu'on était loin de eps pour W et Om quand J a stoppé ?
diffOm=do.call(rbind, lapply(simus, function(graineJ){
  diffOm = tail(graineJ$VEM_1$features$diffOm,1)
}))
sum(diffOm[which(stopJ)]<1e-3) # en fait tous les Omega ont convergé
diffW=do.call(rbind, lapply(simus, function(graineJ){
  diffW = tail(graineJ$VEM_1$features$diffW,1)
}))
hist(log(diffW[which(stopJ)]))
sum(diffW[which(stopJ)]>0.1) # 79
weird.conv = which(do.call(rbind, lapply(simus, function(graine){
  diffW = tail(graine$VEM_1$features$diffW,1)
  return(diffW > 0.1)
})))
# donc divergence des poids W pour 40% des graines
data.frame(diffW= diffW, diffOm= diffOm, nH=misAct$nH) %>% as_tibble() %>% 
  group_by(nH) %>% summarise(count=n(), noWconv = sum(diffW>1), prop=noWconv/count)
# iter des convergences bizarres :
weird.iter=do.call(rbind, lapply(simus[weird.conv], function(graine){
  graine$VEM_1$finalIter
})) # longue queue : entre 13 et 135, mean=41.5 sd = 25
good.iter=do.call(rbind, lapply(simus[-weird.conv], function(graine){
  graine$VEM_1$finalIter
})) # entre 12 et 180, mean = 32.5, sd = 18, (31 et 12.5 si on enlève la valeur extrême)
data.frame(iter = rbind(weird.iter, good.iter), 
           Wconv = as.factor(rep(c(0,1),c(length(weird.iter),length(good.iter))) ) )%>% 
  ggplot(aes(Wconv, iter, color=Wconv))+geom_boxplot()+theme_light()


#-------------------
# Qualité d'inférence
#-- réseau 
AUC<-data.frame(auc=do.call(rbind, lapply(simus, function(seed){
  Pg=seed$VEM_1$Pg
  ome=seed$omega
  diag(ome)=0
  auc=round(auc(pred = Pg, label = ome),4)
  return(auc)
})), seed = 1:200) %>% as_tibble() %>% mutate(Wconv = as.factor(1*!(seed %in% weird.conv)))
AUC %>% ggplot(aes(y=Wconv, x=auc, color=Wconv, fill=Wconv))+geom_density_ridges(alpha=0.3)+
  geom_vline(xintercept = 0.5, linetype="dashed")+
  theme_light()
# partie non liée au noeud H
AUC_O<-data.frame(auc0=do.call(rbind, lapply(simus, function(seed){
  Pg=seed$VEM_1$Pg[-H,-H]
  ome=seed$omega[-H,-H]
  diag(ome)=0
  auc=round(auc(pred = Pg, label = ome),4)
  return(auc)
})), seed = 1:200) %>% as_tibble() 

left_join(AUC, AUC_O, by="seed")%>% gather(key, value, -seed, -Wconv, - nH) %>% 
  ggplot(aes(y=key, x=value, color=key, fill=key))+geom_density_ridges(alpha=0.3)+
  theme_light() # la bosse vient de la partie OH
AUC$nH = misAct$nH
AUC %>% group_by(nH) %>% summarise( count=n(),prop80 = sum(auc>0.8)/count,
                                    prop60 = sum(auc>0.6)/count)
#-- position de l'acteur manquant dans le réseau
PPVH<-data.frame(ppvh=do.call(rbind, lapply(simus, function(seed){
  Pg=seed$VEM_1$Pg
  ome=seed$omega
  diag(ome)=0
  H=15
  ppvh=accppvtpr(Pg,ome,h=H,seuil=0.5)[5]
  return(ppvh)
})),seed=1:200) %>% as_tibble()  %>% mutate(Wconv = as.factor(1*!(seed %in% weird.conv)))
PPVH %>% ggplot(aes(y=Wconv, x=ppvh, color=Wconv, fill=Wconv))+geom_density_ridges(alpha=0.3)+
  geom_vline(xintercept = 0.5, linetype="dashed")+
  theme_light()
PPVH$nH = misAct$nH
PPVH %>% group_by(nH) %>% summarise( count=n(),prop80 = sum(ppvh>0.8)/count,
                                    prop60 = sum(ppvh>0.6)/count)
PPVH %>% mutate(major=as.factor(nH>6)) %>% 
  ggplot(aes(y=major, x=ppvh, color=major, fill=major))+geom_density_ridges(alpha=0.3)+
  geom_vline(xintercept = 0.5, linetype="dashed")+
  theme_light()
AUC %>% mutate(major=as.factor(nH>6)) %>% 
  ggplot(aes(y=major, x=auc, color=major, fill=major))+geom_density_ridges(alpha=0.3)+
  geom_vline(xintercept = 0.5, linetype="dashed")+
  theme_light()

data=data.frame(ppvh=PPVH$ppvh, auc=AUC$auc, nH=misAct$nH) %>%
  gather(key, value, -nH)%>%
  mutate(influence=unlist(purrr::map(nH, function(x){
    if(x<5) res="minor"
    if(x>=5 & x<7) res="medium"
    if(x>=7) res="major"
    return(res)}))) %>% as_tibble()
data %>% group_by(influence) %>% summarise(count=n())
ggplot(data,aes(key, value, color=key, fille=key))+geom_boxplot(alpha=0.3)+
  facet_grid(~as.factor(influence))+theme_light()

weird.auc=which(AUC$auc < 0.5 & AUC$nH >6)
bad=intersect(weird.conv, weird.auc)
sort(diffW[bad],decreasing = TRUE )
sort(diffW[which(misAct$nH>6)],decreasing = TRUE )[1:20]
mis1=misAct[weird.auc,]
mis2=misAct[-weird.auc,]
hist(mis1$btwH)
hist(mis2$btwH)
hist(misAct$btwH/misAct$nH)
AUC %>% mutate(major=as.factor(misAct$mean.cpH>(-0.3))) %>% 
  ggplot(aes(y=major, x=auc, color=major, fill=major))+geom_density_ridges(alpha=0.3)+
  geom_vline(xintercept = 0.5, linetype="dashed")+
  theme_light()
# Les mauvaises convergences de W n'empêchent pas de bonnes estimations, et n'expliquent pas les cas mal 
# inférés.
left_join(AUC, PPVH, by=c("seed","Wconv")) %>% ggplot(aes(auc, ppvh))+geom_point(color="blue", alpha=0.5)+
  theme_light()+geom_vline(xintercept = 0.5, linetype="dashed")+geom_hline(yintercept = 0.5, linetype="dashed")

#-- reconstruction de l'acteur manquant

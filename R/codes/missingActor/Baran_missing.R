source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")
library(ade4)
data(baran95)
counts=as.matrix(baran95$fau)

# Lancer VEM avec 2 manquants et représenter les sites sur le plan des 2 vecteurs de moyenne 
# initialization
p=ncol(counts) ; n=nrow(counts) ; B=100 ; O=1:p
PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M  
SO<-PLNfit$var_par$S  
sigma_obs=PLNfit$model_par$Sigma
theta=PLNfit$model_par$Theta 
matcovar=matrix(1, n,1)

#------ r=2
# r=2 ; H=(p+1):(p+r)
# cliques_spcaY2 <- boot_FitSparsePCA(scale(counts),B, r=r)
# tic()
# ListVEM2<-List.VEM(cliqueList=cliques_spcaY2, counts, sigma_obs, MO,SO,r=r, cores=3)
# toc() # 60.981s sur 3 coeurs avec B=100 (50 pour B=50)
# #------ r=1
# r=1 ;
# cliques_spcaY1 <- boot_FitSparsePCA(scale(counts),B, r=r)
# tic()
# ListVEM1<-List.VEM(cliqueList=cliques_spcaY1, counts, sigma_obs, MO,SO,r=r, cores=3)
# toc() # 20.346s with 3 cores et B=100

##########
#saveRDS(list(ListVEM1,ListVEM2,cliques_spcaY1,cliques_spcaY2),file="/Users/raphaellemomal/these/R/save_simu_baran95/FatalaVEM.rds")
##########

list=readRDS("/Users/raphaellemomal/these/R/save_simu_baran95/FatalaVEM.rds")
ListVEM1=list[[1]]
ListVEM2=list[[2]]
cliques_spcaY1=list[[3]]
cliques_spcaY2=list[[4]]
vBICs1<-as.numeric(vec.vBIC(ListVEM1,counts,theta, matcovar,r))
vBICs2<-as.numeric(vec.vBIC(ListVEM2,counts,theta, matcovar,r))
vBICdata=data.frame(values=c(vBICs1, vBICs2))
vBICdata$r=rep(c("1","2"),c(length(vBICs1),length(vBICs2)))
p1<-vBICdata %>%
  ggplot(aes(r,values, fill=r, color=r))+
  geom_boxplot(alpha=0.3, width=0.5)+mytheme.dark+guides(color=FALSE, fill=FALSE)+
  labs(x="number of missing actors",y="values",title="vBIC from boot.sPCA (B=100)")
best1=which.max(vBICs1)
best2=which.max(vBICs2)
vBICs1[best1]
vBICs2[best2]
#best cliques
cliques_spcaY2[[best2]]
#best missing actors
r=2
H=(p+1):(p+r)
BaranVEM<-ListVEM2[[best2]]
MH=BaranVEM$M[,H]
p2<-MH %>% as_tibble() %>% mutate(site=baran95$plan$site, date=baran95$plan$date) %>% 
  mutate(covar=paste(site, date), covarnum=as.numeric(as.factor(covar))) %>% 
  ggplot(aes(V1, V2, color=site))+geom_point()+mytheme.dark+
  labs(x="Actor 1",y="Actor 2", title="Means for best VEMtree with r=2")
plot=grid.arrange(p1, p2, ncol=2)
ggsave(filename = "Barans_missing_B100.png", plot = plot,
       path ="/Users/raphaellemomal/these/R/images/", width = 9, height = 4)

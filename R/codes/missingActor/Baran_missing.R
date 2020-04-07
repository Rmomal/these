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
# alpha<-tryCatch(expr={computeAlpha(solve(sigma_obs),default =1/n, MO, SO, plot=plot)},
#                 error=function(e){message("sythetic alpha")
#                   return(1/n)})
alpha=1
#------ r=5
r=5 ; H=(p+1):(p+r)
cliques_spcaY5 <- boot_FitSparsePCA(scale(counts),B, r=r)
tic()
ListVEM5<-List.VEM(cliqueList=cliques_spcaY5, counts, sigma_obs, MO,SO,alpha,r=r,maxIter = 100,
                   cores=3)
toc() # 60.981s sur 3 coeurs avec B=100 (50 pour B=50)
converged<-do.call(rbind,lapply(ListVEM5,function(vem){
  diffW=vem$features$diffW
  conv=(diffW[length(diffW)]<=1e-4)}))
if(sum(converged)!=0){
  ListVEM5=ListVEM5[converged]
}
# test
best=which.max(crit5$vBIC)
VEM_5<-ListVEM5[[best]]
#------ r=4
r=4 ; H=(p+1):(p+r)
cliques_spcaY4 <- boot_FitSparsePCA(scale(counts),B, r=r)
tic()
ListVEM4<-List.VEM(cliqueList=cliques_spcaY4, counts, sigma_obs, MO,SO,alpha,r=r,maxIter = 100,
                   cores=3)
toc() # 60.981s sur 3 coeurs avec B=100 (50 pour B=50)
converged<-do.call(rbind,lapply(ListVEM4,function(vem){
  diffW=vem$features$diffW
  conv=(diffW[length(diffW)]<=1e-4)}))
if(sum(converged)!=0){
  ListVEM4=ListVEM4[converged]
}
#test
best=which.max(crit4$vBIC)
VEM_4<-ListVEM4[[best]]
#------ r=3
r=3 ; H=(p+1):(p+r)
cliques_spcaY3 <- boot_FitSparsePCA(scale(counts),B, r=r)
tic()
ListVEM3<-List.VEM(cliqueList=cliques_spcaY3, counts, sigma_obs, MO,SO,alpha,r=r,maxIter = 100,
                   cores=3)
toc() # 228.9s sur 3 coeurs avec B=100 (50 pour B=50)
converged<-do.call(rbind,lapply(ListVEM3,function(vem){
  diffW=vem$features$diffW
  conv=(diffW[length(diffW)]<=1e-4)}))
if(sum(converged)!=0){
  ListVEM3=ListVEM3[converged]
}
#------ r=2
r=2 ; H=(p+1):(p+r)
cliques_spcaY2 <- boot_FitSparsePCA(scale(counts),B, r=r)
tic()
ListVEM2<-List.VEM(cliqueList=cliques_spcaY2, counts, sigma_obs, MO,SO,alpha,r=r,maxIter = 100, cores=3)
toc() # 183.9s sur 3 coeurs avec B=100 (50 pour B=50)
converged<-do.call(rbind,lapply(ListVEM2,function(vem){
  diffW=vem$features$diffW
  conv=(diffW[length(diffW)]<=1e-4)}))
if(sum(converged)!=0){
  ListVEM2=ListVEM2[converged]
}
#------ r=1
r=1 ;
cliques_spcaY1 <- boot_FitSparsePCA(scale(counts),B, r=r)
tic()
ListVEM1<-List.VEM(cliqueList=cliques_spcaY1, counts, sigma_obs, MO,SO,alpha,r=r,maxIter = 100, cores=3)
toc() # 62s with 3 cores et B=100
converged<-do.call(rbind,lapply(ListVEM1,function(vem){
  diffW=vem$features$diffW
  conv=(diffW[length(diffW)]<=1e-4)}))
if(sum(converged)!=0){
  ListVEM1=ListVEM1[converged]
} 
# test 
# cliqueinit=cliques_spcaY1[[1]]
# init1=initVEM(counts = counts, initviasigma = cliqueinit,  sigma_obs,r = 1)
# Wginit= init1$Wginit; Winit= init1$Winit; omegainit=init1$omegainit ; MHinit=init1$MHinit
# 
# VEM_1<-VEMtree(counts, MO, SO, MH=MHinit,ome_init = omegainit,W_init =Winit,eps=1e-3,
#                Wg_init =Wginit,plot = TRUE, maxIter = 100,print.hist = FALSE,
#                vraiOm = NULL, alpha=1, verbatim=FALSE , filterPg = TRUE,filterWg = TRUE)

#------ r=0
r=0;
init0=initVEM(counts = counts, initviasigma = NULL,  sigma_obs,r = 0)
Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 
VEM_0<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=1e-4,
               Wg_init =Wginit,plot = TRUE, maxIter = 100,print.hist = FALSE,
               vraiOm = NULL, alpha=alpha, verbatim=FALSE , filterPg = FALSE,filterWg = TRUE)
# lower bound check
VEM_5$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>%
  ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(parameter)), size=3)+geom_line()+
  facet_wrap(~key, scales="free")+
  labs(x="iteration",y="", title="Lower bound and components")+mytheme+
  scale_color_discrete("")+coord_cartesian(xlim=c(0,50))
VEM_5$features  %>%  rowid_to_column() %>%
  pivot_longer(-rowid,names_to="key",values_to = "values") %>%
  ggplot(aes(rowid,values, color=key))+
  geom_point()+geom_line() + facet_wrap(~key, scales="free")+
  labs(x="",y="", title="Parameters")+ mytheme.dark+guides(color=FALSE)
ggimage(VEM_5$Pg)
VEM_5$features$diffW
##########
saveRDS(list(VEM_0,ListVEM1,ListVEM2,ListVEM3,ListVEM4,ListVEM5,
             cliques_spcaY1,cliques_spcaY2, cliques_spcaY3, 
             cliques_spcaY4,cliques_spcaY5),file="/Users/raphaellemomal/these/R/save_simu_baran95/FatalaVEM.rds")
##########

list=readRDS("/Users/raphaellemomal/these/R/save_simu_baran95/FatalaVEM.rds")
VEM_0<-list[[1]]
ListVEM1=list[[2]]
ListVEM2=list[[3]]
ListVEM3=list[[4]]
ListVEM4=list[[5]]
ListVEM5=list[[6]]
cliques_spcaY1=list[[7]]
cliques_spcaY2=list[[8]]
cliques_spcaY3=list[[9]]
cliques_spcaY4=list[[10]]
cliques_spcaY5=list[[11]]
crit0<-criteria(list(VEM_0),counts=counts,theta = theta,matcovar = matcovar, r=0)
crit1<-criteria((ListVEM1),counts=counts,theta = theta,matcovar = matcovar, r=1)
crit2<-criteria(ListVEM2,counts=counts,theta = theta,matcovar = matcovar, r=2)
crit3<-criteria(ListVEM3,counts=counts,theta = theta,matcovar = matcovar, r=3)
crit4<-criteria(ListVEM4,counts=counts,theta = theta,matcovar = matcovar, r=4)
crit5<-criteria(ListVEM5,counts=counts,theta = theta,matcovar = matcovar, r=5)
# 36, 24, 20, 13, 16 : tailles des listes de 1 à 5
crit=rbind(crit0,crit1, crit2, crit3, crit4, crit5)
plot=crit %>% gather(key, value, -r) %>% 
  ggplot(aes(as.factor(r), value, color=key, fill=key))+geom_point(alpha=0.3)+
  facet_wrap(~key)+mytheme.dark+
  labs(x="number of missing actors",y="values",title="Crit from boot.sPCA (B=100, alpha=1, eps=1e-4, kept only converged)")
 
ggsave(filename = "Bad_Fatala_crit_B100.png", plot = plot,
       path ="/Users/raphaellemomal/these/R/images/", width = 6, height = 5)

best1=which.max(vBICs1)
best2=which.max(vBICs2)
bests<-c(which.max(vBICs1),which.max(vBICs2),which.max(vBICs3),
         which.max(vBICs4),which.max(vBICs5))
bestVEM<-list(ListVEM1[bests[1]],ListVEM2[bests[2]],
         ListVEM3[bests[3]],ListVEM4[bests[4]],
         ListVEM5[bests[5]])
bestVEM[[1]][[1]]$lowbound$J
bestVEM[[2]][[1]]$lowbound$J
bestVEM[[3]][[1]]$lowbound$J
bestVEM[[4]][[1]]$lowbound$J
bestVEM[[5]][[1]]$lowbound$J
r=5 ; d=1
nbparam<-2*p*r - 1
penvBIC=nbparam*log(n)/2

#low bounds of best cliques
do.call(grid.arrange,lapply(bestVEM, function(vem){
  r=nrow(vem[[1]]$omega)-ncol(counts)
    finalJ=round(tail(vem[[1]]$lowbound$J, 1),2)
  vem[[1]]$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>% 
    filter(key=="J") %>% 
    ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(parameter)), size=3)+geom_line()+
    facet_wrap(~key, scales="free")+
    labs(x="iteration",y="", title=paste0("final J=",finalJ," with r=",r))+mytheme+
    scale_color_discrete("")
}))
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

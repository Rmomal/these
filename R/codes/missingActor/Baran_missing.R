source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")
library(ade4)
data(baran95)
raw_counts=as.matrix(baran95$fau)
poor_species<-which(colSums(raw_counts==0)/nrow(raw_counts)>0.9)
counts=as.matrix(baran95$fau[,-poor_species])

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

D=.Machine$double.xmax
f<-function(u,q){u*exp(u)-D^(1/(q-1))/sqrt(q*(q-1))} # f est monotone croissante
fprim=function(u,q){
  exp(u)*(u-1)-u*D^(1/(q-1))/sqrt(q*(q-1))
}
minimum=optimize(fprim, c(0,100),maximum = FALSE, q=(p+3)) # max 3 missing actors
alpha_sup=minimum$minimum/n

#------ 

cliques_Fatala<-function(r,MO){
  set.seed(1); p=ncol(MO)
  cliques_spca <- boot_FitSparsePCA(scale(MO),B, r=r, cores = 3)
  saveRDS(cliques_spca,paste0("/Users/raphaellemomal/simulations/cliques_",r,".rds"))
  return(cliques_spca)
}
VEMr<-function(r, alpha, cores=3){
  set.seed(1)
  cliques_spca <- readRDS(paste0("/Users/raphaellemomal/simulations/cliques_",r,".rds"))
  tic()
  cat("VEM")
  ListVEM<-List.VEM(cliquesObj=cliques_spca, counts=counts, sigma_obs=sigma_obs,MO= MO,
                    SO=SO,alpha=alpha,r=r,
                    maxIter = 100,filterDiag = TRUE, filterWg = TRUE,nobeta = FALSE,
                    cores=cores, eps = 1e-3, save=TRUE)
  toc() # 60.981s sur 3 coeurs avec B=100 (50 pour B=50)
  return(ListVEM)
}

for(i in 1:3){
  tic()
  cliques_Fatala(i,MO)
  toc()
}
#----- RUN
#-- comparaison EM0 et VEM0
mean.sigT=mean((1/n)*(t(MO)%*%MO+diag(colSums(SO))))
EM0<-EMtree(PLNfit,maxIter = 20,cond.tol = 1e-6,plot = TRUE)
init0=initVEM(counts , initviasigma = NULL,  sigma_obs,r = 0)
Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 
VEM0<-VEMtree(counts,MO,SO,MH=NULL,omegainit,Winit,Wginit, eps=1e-3, 
              alpha=0.1,maxIter=20, plot=TRUE,print.hist=FALSE,filterWg=TRUE,
              verbatim = TRUE,nobeta = FALSE, filterDiag = TRUE)
g1<-ggimage(VEM0$Pg)
g2<-ggimage(EM0$edges_prob)
grid.arrange(g1, g2, ncol=2)
#-- 
alpha=round(alpha_sup,1)-0.1
VEM1<-VEMr(1, alpha=0.1)
VEM2<-VEMr(2, alpha=0.1)
VEM3<-VEMr(3, alpha=0.1)

#----- 
# bug coeur 3, relance
# mclapply(c(9,12,15,18,90), function(num){
# #  file=paste0("/Users/raphaellemomal/simulations/Fatala_missing/r3_num",num,".rds")
#  # if(!file.exists(file)){
#    # cat(paste0("num ", num," did not work in parallel\n"))
#   num=9
#     cliquesObj = readRDS(paste0("/Users/raphaellemomal/simulations/cliques_3.rds"))
#     c=cliquesObj$cliqueList[[num]]
#     init=initVEM(counts = counts, initviasigma=c, sigma_obs,r = 3)
#     Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
#     #run VEMtree
#     VEM<-VEMtree(counts=counts,MO=MO,SO=SO,MH=MHinit,ome_init=omegainit,W_init=Winit,
#                  Wg_init=Wginit, eps=1e-2, alpha=0.1,maxIter=100, 
#                  verbatim = TRUE, print.hist=FALSE, filterWg = TRUE, 
#                  filterDiag = TRUE, nobeta=FALSE)
#     VEM$clique=c
#     VEM$nbocc=cliquesObj$nb_occ[num]
#     VEM$nbocc_vec=cliquesObj$nb_occ
#     saveRDS(VEM, paste0("/Users/raphaellemomal/simulations/Fatala_missing/r3_num",num,".rds"))
#  # }
# }, mc.cores=3)
#######
# tout marche sauf la 9e clique pour r=3 : Error in d2q(as.matrix(A)) : argument not finite-valued
#######

#-----
# gather data
Fatala_1<-list()
Fatala_2<-list()
Fatala_3<-list()

Fatala_3<-lapply((1:90)[-c(9)], function(num){
  file =  paste0("/Users/raphaellemomal/simulations/Fatala_missing/r3_num",num,".rds")
  readRDS(file)
})
Fatala_selec=list(r1=Fatala_1,r2=Fatala_2,r3=Fatala_3)
saveRDS(Fatala_selec,"/Users/raphaellemomal/simulations/Fatala_missing/Fatala_selec.rds" )

# check all went right
mean(do.call(rbind,lapply(Fatala_selec[[1]], length)))
mean(do.call(rbind,lapply(Fatala_selec[[2]], length)))
mean(do.call(rbind,lapply(Fatala_selec[[3]], length))) #ok.

#----- initialization selection

JData123<-lapply(1:3, function(r){
  do.call(rbind,lapply(seq_along(Fatala_selec[[r]]),function(num){
    return(c(getJcor(Fatala_selec[[r]][[num]],p=ncol(MO)),num=num))
  }))}) # 2min
JData123=do.call(rbind,JData123) %>%  as_tibble()
JData0<-c(getJcor(VEM0,p=ncol(MO)),num=1)
Jdata = rbind(JData0, JData123)
Jdata %>% ggplot(aes(r, Jcor, color=as.factor(r)))+geom_point()+mytheme.dark("")

penT123<- lapply(1:3, function(r){ 
  do.call(rbind,lapply(seq_along(Fatala_selec[[r]]),function(num){
    vem=Fatala_selec[[r]][[num]]
    penT=-( sum( vem$Pg * log(vem$Wg+(vem$Wg==0)) ) - logSumTree(vem$Wg)$det) 
    return(c(penT=penT, r=r, num=num))}))})
penT123=do.call(rbind,penT123) %>% as_tibble()
penT0=c(-( sum( VEM0$Pg * log(VEM0$Wg+(VEM0$Wg==0)) ) - logSumTree(VEM0$Wg)$det) ,r=0,num=1)
penT=rbind(penT0, penT123)

Fatala_data=left_join(Jdata, penT, by=c("r","num"))
p=ncol(MO) ; n=nrow(MO)
Fatala_data=Fatala_data %>% mutate(penZH = 0.5*sum(log(rep(1,n))) + n*r*0.5*(1+log(2*pi)),
                       penvBIC=log(n)*0.5*(p+ (p*(p+1)/2 +r*p)+((p+r)*(p+r-1)/2 - 1))) %>% 
  mutate(ICL0 = Jcor - penvBIC ,
         ICL1 = Jcor - penvBIC - penZH,
         ICL2 = Jcor - penvBIC - penZH - penT )

Fatala_data %>%ggplot(aes(r, ICL, color=as.factor(r)))+geom_point()+mytheme.dark("")

Fatala_data %>%dplyr::select(r, Jcor, ICL0, ICL1, ICL2) %>% gather(key, value,-r,-Jcor) %>% 
  ggplot(aes(Jcor,value, color=as.factor(r)))+geom_abline()+geom_point()+ facet_wrap(~key)+
  mytheme.dark("")
Fatala_data %>%dplyr::select(r, ICL0, ICL1, ICL2) %>%   
  group_by(r) %>% mutate(maxICL0 = max(ICL0),maxICL1 = max(ICL1),maxICL2 = max(ICL2)) %>% ungroup() %>% 
  dplyr::select(r,maxICL0,maxICL1,maxICL2) %>% gather(key, value,-r) %>% 
  ggplot(aes(r,value, color=as.factor(r)))+geom_point(size=3)+ facet_wrap(~key)+
  mytheme.dark("")

#-----
p=ncol(MO); r=2 ; H=(p+1):(p+r)
bestnum=unlist(Fatala_data %>% filter(r==r) %>% filter(ICL0==max(ICL0)) %>% dplyr::select(num))
bestvem = Fatala_selec[[2]][[10]]
MH=bestvem$M[,H]
Fatala_means = data.frame(MH,baran95$plan$site)
colnames(Fatala_means) = c("M1","M2","Site")
ggplot(Fatala_means, aes(M1, M2, color=Site))+geom_point()+mytheme.dark("")
##########
# anciens résultats pour mémoire
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
crit=rbind(crit0,crit1, crit2, crit3, crit4)
plot=crit %>% gather(key, value, -r) %>% 
  ggplot(aes(as.factor(r), value, color=key, fill=key))+geom_boxplot(alpha=0.3)+
  facet_wrap(~key, scales="free")+mytheme.dark("")+
  labs(x="number of missing actors",y="values",title="Crit from boot.sPCA (B=100, alpha=1, eps=1e-4, kept only converged)")

ggsave(filename = "Bad_Fatala_ICL_B100.png", plot = plot,
       path ="/Users/raphaellemomal/these/R/images/", width = 7, height = 3)

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

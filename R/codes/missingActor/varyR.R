library(EMtree)
library(PLNmodels)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(mvtnorm)
library(useful)
library(mclust)
library(MASS)
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(kableExtra)
library(parallel)
library(sparsepca)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-exactDet.R")

#-- function
simu_vary_r<-function(p,n,r,B,rMax,maxIter,seed,type,alpha,eps=1e-3, cores=3, plot){
  set.seed(seed)
  # sim data
  missing_data<-missing_from_scratch(n,p,r=r,type,plot)
  counts=missing_data$Y 
  # Observed parameters
  PLNfit<-PLN(counts~1, control=list(trace=0))
  MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S ; 
  theta=PLNfit$model_par$Theta; sigma_obs=PLNfit$model_par$Sigma
  #-- normalize the PLN outputs
  D=diag(sigma_obs)
  matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
  MO=MO*matsig
  SO=SO*matsig^2
  #initialize
  init0=initVEM(counts , initviasigma = NULL,  cov2cor(sigma_obs),r = 0)
  Wginit= init0$Wginit; Winit= init0$Winit; upsinit=init0$upsinit 
  # vem with original counts
  VEM_r0<-VEMtree(counts, MO, SO, MH=NULL,upsinit,W_init =Winit,eps=eps,
                  Wg_init =Wginit,plot = FALSE, maxIter = maxIter,print.hist = FALSE,
                  alpha=alpha, verbatim=FALSE, trackJ=FALSE )
  
  #  vary r
  VEMr<-lapply(1:rMax, function(r){
    cat(paste0(r," missing actors: "))
    t1<-Sys.time()
    cliques_spca <- boot_FitSparsePCA(scale(MO),B,r=r,cores=3)
    ListVEM<-List.VEM(cliquesObj =cliques_spca, counts, cov2cor(sigma_obs), MO,SO,r=r,alpha=alpha, cores=cores,
                      eps=eps,maxIter,trackJ=FALSE)
    t2<-Sys.time()
    runtime=difftime(t2,t1)
    cat(paste0(round(runtime,3), attr(runtime, "units"),"\n"))
    return(ListVEM)
  })
  saveRDS(list(VEM_r0=VEM_r0, VEMr=VEMr, counts=counts, theta=theta),
          file=paste0("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/",type,"_seed",
                      seed,"_r",r,"_1-",rMax ,".rds") )
}

############
#-- run
seed=1; p=14 ; B=100; type="scale-free" ; n=200 
t1<-Sys.time()
simu_vary_r(seed=seed,B=B, n=n, p=p,r=0,maxIter=200,rMax=2, 
            type = type,nobeta=FALSE, plot=FALSE , cores=3, alpha=0.1)
t2<-Sys.time()
difftime(t2, t1)
t1<-Sys.time()
simu_vary_r(seed=seed,B=B, n=n, p=p,r=1,maxIter=200,
            rMax=2, type = type, nobeta=FALSE, plot=FALSE, cores=3, alpha=0.1)
t2<-Sys.time()
difftime(t2, t1)

#--

get_allData<-function(seed, rMax,eps1=1e-14,eps2=1e-14){
  simr1<-readRDS(paste0("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/scale-free_seed",seed,"_r1_1-",rMax,".rds") )       
  simr0<-readRDS(paste0("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/scale-free_seed",seed,"_r0_1-",rMax,".rds") )      
  p=14 ;n=200
  Data123<-lapply(1:rMax, function(r){
    do.call(rbind,lapply(seq_along(simr0$VEMr[[r]]),function(num){
      vem=simr0$VEMr[[r]][[num]]
      if(length(vem)==14){
        penT=-( 0.5*sum( vem$Pg * log(vem$Wg+(vem$Wg==0)) ) - logSumTree(vem$Wg)$det) 
        
        penZH=0.5*sum(apply(matrix(vem$S[,(p+1):(p+r)],n,r),2,
                            function(x){ log(sum(x))}))+0.5*r*n*(1+log(2*pi))
        sumPg = sum(vem$Pg)-2*(p+r-1)
        values=c(J=tail(vem$lowbound$J,1),getJcor(vem,p=p,eps1,eps2),penT=penT,penZH=penZH,sumPg=sumPg)
      }else{
        values=c(J=NaN,Jcor=NaN, diff=NaN,detEg=NaN, delta=NaN,penT=NaN,penZH=NaN,sumPg=NaN)
      }
      return(c(values,r=r, trueR=0, num=num))
    }))}) 
  Data123=do.call(rbind,Data123) %>%  as_tibble()
  Data0<-c(J=tail(simr0$VEM_r0$lowbound$J,1),
           getJcor(simr0$VEM_r0,p=14,eps1,eps2),
           penT=-(0.5*sum( simr0$VEM_r0$Pg * log(simr0$VEM_r0$Wg+(simr0$VEM_r0$Wg==0)) ) - 
                    logSumTree(simr0$VEM_r0$Wg)$det) ,
           penZH=0,sumPg=sum(simr0$VEM_r0$Pg)-26,
           r=0, trueR=0, num=1)
  Data_r0 = rbind(Data0, Data123)
  
  Data123<-lapply(1:rMax, function(r){
    do.call(rbind,lapply(seq_along(simr1$VEMr[[r]]),function(num){
      vem=simr1$VEMr[[r]][[num]]
      if(length(vem)==14){
        penT=-( 0.5*sum( vem$Pg * log(vem$Wg+(vem$Wg==0)) ) - logSumTree(vem$Wg)$det) 
        penZH=0.5*sum(apply(matrix(vem$S[,(p+1):(p+r)],n,r),2,
                            function(x){ log(sum(x))}))+0.5*r*n*(1+log(2*pi))
        sumPg = sum(vem$Pg)-2*(p+r-1)
        values=c(J=tail(vem$lowbound$J,1),getJcor(vem,p=p,eps1,eps2),penT=penT,penZH=penZH,sumPg=sumPg)
      }else{
        values=c(J=NaN,Jcor=NaN, diff=NaN,detEg=NaN, delta=NaN,penT=NaN,penZH=NaN,sumPg=NaN)
      }
      return(c(values,r=r, trueR=1, num=num))
    }))}) 
  Data123=do.call(rbind,Data123) %>%  as_tibble()
  Data0<-c(J=tail(simr1$VEM_r0$lowbound$J,1),getJcor(simr1$VEM_r0,p=14,eps1,eps2),
           penT=-(0.5*sum( simr1$VEM_r0$Pg * log(simr1$VEM_r0$Wg+(simr1$VEM_r0$Wg==0)) ) -logSumTree(simr1$VEM_r0$Wg)$det) ,
           penZH=0,sumPg=sum(simr1$VEM_r0$Pg)-26,
           r=0, trueR=1, num=1)
  Data_r1 = rbind(Data0, Data123)
  
  allData = rbind(Data_r0, Data_r1)
  allData=allData %>% mutate(goodsumP = abs(sumPg)<1e-3, 
                             penvBIC=0.5*log(n)*(p+ (p*(p+1)/2 +r*p)+((p+r)*(p+r-1)/2 - 1)),
                             ICL = Jcor-penT-penZH-penvBIC)
  return(allData)
}
allData1<-get_allData(1, rMax=2)
allData2<-get_allData(2, rMax=2)
allData7<-get_allData(7, rMax=3,eps1=1e-8,eps2=1e-7)
summarise=allData1 %>% filter(!is.na(ICL)) %>% #filter(detEg>-25) %>% #filter(goodsumP) %>% 
  group_by(r, trueR) %>%
  summarize( maxICL=max(ICL, na.rm=TRUE),
             maxJCor=max(Jcor, na.rm=TRUE),maxJ=max(J, na.rm=TRUE)) 

summarise %>%   gather(key, value, -r,-trueR) %>% 
  ggplot(aes(as.factor(r), value, color=as.factor(trueR)))+
  facet_wrap(~key, scales="free")+
  geom_point(shape=16, size=3)+mytheme.dark("True r:")+
  labs(x="r")

allData %>% ggplot(aes(Jcor, diff, color=as.factor(trueR)))+geom_point()+
  facet_wrap(~r)+mytheme.dark("")

allData %>% 
  #--- closer look
  vem = simr0$VEMr[[2]][[40]]
vem$lowbound$J

JfiltreDelta=Jdata %>% group_by(trueR, r) %>% 
  mutate(q10Delta=quantile(Delta, 0.1)) %>% 
  filter(Delta<q10Delta)

JfiltreDelta %>% mutate(maxJcor=max(Jcor)) %>% 
  ggplot(aes(as.factor(r), maxJcor, color=as.factor(trueR)))+
  geom_point(shape=16, size=3)+facet_wrap(~as.factor(trueR))+mytheme.dark("True r:")+
  labs(x="r", title="Parmi les petites valeurs de Delta :")
penTr0_0=data.frame(penT=-( sum( simr0$list.vem0$Pg * log(simr0$list.vem0$Wg+(simr0$list.vem0$Wg==0)) ) - 
                              logSumTree(simr0$list.vem0$Wg)$det) ,r=0, trueR = 0)
penT_r0<-data.frame(do.call(rbind, do.call(rbind, lapply(1:3, function(r){ 
  lapply(simr0$VEMr[[r]],function(vem){
    penT=-( sum( vem$Pg * log(vem$Wg+(vem$Wg==0)) ) - logSumTree(vem$Wg)$det) 
    return(c(penT=penT, r=r, trueR=0))})}))))
penTr1_0=data.frame(penT=-( sum( simr1$list.vem0$Pg * log(simr1$list.vem0$Wg+(simr1$list.vem0$Wg==0)) ) - 
                              logSumTree(simr1$list.vem0$Wg)$det) , r=0, trueR=1)
penT_r1<-do.call(rbind, do.call(rbind, lapply(1:3, function(r){ 
  lapply(simr1$VEMr[[r]],function(vem){
    penT=-( sum( vem$Pg * log(vem$Wg+(vem$Wg==0)) ) - logSumTree(vem$Wg)$det) 
    return(c(penT=penT, r=r, trueR=1))})
})))
penT=rbind(penTr0_0, penT_r0,penTr1_0,penT_r1) %>% as_tibble()
left_join(JfiltreDelta, penT, by=c("r","trueR")) %>% 
  mutate(nbparam=(p + (p*(p+1)/2 +r*p+r)+((p+r)*(p+r-1)/2 - 1)),
         ICL = Jcor - penT- nbparam*log(n)/2) %>% group_by(r, trueR) %>% 
  mutate(maxICL = max(ICL)) %>% 
  ggplot(aes(as.factor(r), maxICL, color=as.factor(trueR)))+
  geom_point(shape=16, size=3)+facet_wrap(~as.factor(trueR))+mytheme.dark("True r:")+
  labs(x="r", title="Parmi les petites valeurs de Delta :")


ggsave(plot, filename = "/Users/raphaellemomal/these/R/images/Selec_r.png", width=6, height=4)
# crit0<-do.call(rbind, lapply(simr0$list.vem0, function(vem0){
#  criteria(list(vem0),simr0$counts,simr0$theta, matcovar=matrix(1, nrow(simr0$counts),1),r=0)}))
crit0<- criteria(list(simr0$list.vem0),simr0$counts,simr0$theta,
                 matcovar=matrix(1, nrow(simr0$counts),1),r=0)
critr<-do.call(rbind, lapply(seq_along(simr0$VEMr), function(r){
  criteria(simr0$VEMr[[r]],simr0$counts,simr0$theta, matcovar=matrix(1, nrow(simr0$counts),1),r=r)
}))
critr0=rbind(crit0, critr)
crit0<- criteria(list(simr1$list.vem0),simr1$counts,simr1$theta,
                 matcovar=matrix(1, nrow(simr1$counts),1),r=0)
critr<-do.call(rbind, lapply(seq_along(simr1$VEMr), function(r){
  criteria(simr1$VEMr[[r]],simr1$counts,simr1$theta,  matcovar=matrix(1, nrow(simr1$counts),1),r=r)
}))
critr1=rbind(crit0, critr)
#----------------
# evolution du max de la borne inf avec r
critr0$trueR=0
critr1$trueR=1

critr0$max.prec=c(simr0$list.vem0$max.prec, do.call(rbind, lapply(simr0$VEMr, function(r){
  do.call(rbind, lapply(r, function(vem){
    vem$max.prec
  }))
})))
critr1$max.prec=c(simr1$list.vem0$max.prec, do.call(rbind, lapply(simr1$VEMr, function(r){
  do.call(rbind, lapply(r, function(vem){
    vem$max.prec
  }))
})))
critr0$nbocc=c(0, do.call(rbind, lapply(simr0$VEMr, function(r){
  do.call(rbind, lapply(r, function(vem){
    vem$nbocc
  }))
})))
critr1$nbocc=c(0, do.call(rbind, lapply(simr1$VEMr, function(r){
  do.call(rbind, lapply(r, function(vem){
    vem$nbocc
  }))
})))
crit=rbind(critr0,critr1) %>% as_tibble()
crit %>% ggplot(aes( J,ICL, color=(max.prec)))+geom_point()+geom_abline()+
  facet_grid(trueR~r)+mytheme.dark("")

testAUC<-crit %>% filter(trueR==1 & r==1)

set.seed(1)
missing_data<-missing_from_scratch(n,p,r=r,type,plot)
hidden=missing_data$H ; q=15
ome=omega[c(setdiff(1:q, hidden), hidden),c(setdiff(1:q, hidden), hidden)]
diag(ome)=0
AUC_critr<-do.call(rbind, lapply(simr1$VEMr[[1]], function(vem){
  Pg=vem$Pg
  auc<-round(auc(pred = Pg, label = ome),3)
}))
testAUC$AUC=AUC_critr
PPVH_critr<-do.call(rbind, lapply(simr1$VEMr[[1]], function(vem){
  Pg=vem$Pg
  ppvh=  accppvtpr(Pg,ome,h=15,seuil=0.5)[5]
}))
testAUC$ppvh=PPVH_critr
testAUC %>% ggplot(aes(AUC, vBIC, color=max.prec))+geom_point()+mytheme.dark("")


crit %>%as_tibble() %>% 
  group_by(trueR,r) %>%  summarise(max=max(ICL)) %>%
  ggplot(aes(r,max,color=as.factor(trueR)))+geom_point()+geom_line()+
  mytheme.dark("True r:")+  labs(x="r",y="max Lower Bound")

ggsave(filename = "EvolMaxJ_r.png", plot = plot,
       path ="/Users/raphaellemomal/these/R/images/", width = 4, height = 3)

#----------------
# nombre d'initialisations qui ont convergé pour chaque r
nbconv0<-data.frame(nbconv=do.call(rbind, lapply(simr0$VEMr, length)),
                    r=1:4,trueR=0)
nbconv1<-data.frame(nbconv=do.call(rbind, lapply(simr1$VEMr, length)),
                    r=1:4,trueR=1)
conv=rbind(nbconv0, nbconv1) 
conv %>% as_tibble() %>% 
  ggplot(aes(r,nbconv,color=as.factor(trueR),fill=as.factor(trueR)))+
  geom_col(alpha=0.3,position = "dodge")+mytheme.dark("True r:")
#----------------
# nombre de trim par r
nbtrim0<-data.frame(nbtrim=do.call(rbind, lapply(simr0$VEMr, function(r){
  sum(do.call(rbind,lapply(r, function(vem){vem$trim })))
})),  r=1:4,trueR=0)
nbtrim1<-data.frame(nbtrim=do.call(rbind, lapply(simr1$VEMr, function(r){
  sum(do.call(rbind,lapply(r, function(vem){vem$trim })))
})),  r=1:4,trueR=1)
trim=rbind(nbtrim0, nbtrim1) 
trim %>% as_tibble() %>% 
  ggplot(aes(r,nbtrim,color=as.factor(trueR),fill=as.factor(trueR)))+
  geom_col(alpha=0.3,position = "dodge")+mytheme.dark("True r:")
left_join(conv, trim, by=c("trueR","r")) %>% 
  ggplot(aes(nbconv, nbtrim, color=as.factor(trueR)))+geom_abline(linetype="dashed",color="gray")+
  geom_point(aes(shape=as.factor(r)))+geom_line()+
  mytheme.dark("True r:")

# meilleure J sans trim  
# r=3 pour trueR=1
data0=do.call(rbind,lapply(simr1$VEMr[[3]], function(vem){
  data.frame(trim=vem$trim,J=tail(vem$lowbound$J,1))})) 
data0$trueR=0

# r=2 pour trueR=0, exemple y a un trim meilleur
data1=do.call(rbind,lapply(simr0$VEMr[[2]], function(vem){
  data.frame(trim=vem$trim,J=tail(vem$lowbound$J,1))})) 
data1$trueR=1
rbind(data0,data1) %>% 
  ggplot(aes(as.factor(trueR), J, color=as.factor(trim)))+geom_boxplot()+mytheme.dark("")

#----------------
# nb_occ par r
nbocc0<-do.call(rbind, lapply(seq_along(simr0$VEMr), function(r){
  do.call(rbind,lapply(simr0$VEMr[[r]], function(vem){
    data.frame(nbocc=vem$nbocc ,trim=vem$trim,J=tail(vem$lowbound$J,1),r=r)}))
}))
nbocc0$trueR=0
nbocc1<-do.call(rbind, lapply(seq_along(simr1$VEMr), function(r){
  do.call(rbind,lapply(simr1$VEMr[[r]], function(vem){
    data.frame(nbocc=vem$nbocc ,trim=vem$trim,J=tail(vem$lowbound$J,1),r=r)}))
}))
nbocc1$trueR=1

nbocc0 %>% ggplot(aes(nbocc,J,color=as.factor(trim)))+geom_point()+
  facet_wrap(~r, scales="free")+mytheme.dark("")
nbocc1 %>% ggplot(aes(nbocc,J,color=as.factor(trim)))+geom_point()+
  facet_wrap(~r, scales="free")+mytheme.dark("")

#----------------
p1<-critr0 %>%as_tibble() %>% filter(r<4) %>%  gather(key, value, -r ) %>% 
  filter(key%in%c("vBIC","J")) %>% 
  ggplot(aes(as.factor(r), value, color=as.factor(r==0), fill=as.factor(r==0)))+
  geom_boxplot(alpha=0.3)+mytheme.dark("")+
  facet_wrap(~key, scale="free")+guides(color=FALSE,fill=FALSE)+
  labs(x="number of missing actors",y="",title="seed=1, r=0")
p2<-critr1 %>% as_tibble() %>% filter(r<4) %>%gather(key, value, -r ) %>% filter(key%in%c("vBIC","J")) %>% 
  ggplot(aes(as.factor(r), value, color=as.factor(r==1), fill=as.factor(r==1)))+
  geom_boxplot(alpha=0.3)+mytheme.dark("")+
  facet_wrap(~key, scale="free")+guides(color=FALSE,fill=FALSE)+
  labs(x="number of missing actors",y="",title="seed=1, r=1")
grid.arrange(p1, p2, ncol=2)
#----------------
stat=critr1 %>% as_tibble() %>% 
  mutate(penICL = J-ICL) %>% group_by(r) %>% 
  summarise(mvBIC=median(vBIC),mICL=median(ICL),mJ=median(J),
            mpenICL=median(penICL),maxICL=max(ICL), nr=n())
model=lm(stat$mJ~stat$r)
quant=summary(model)
quant$coefficients[2]/(p*n)

ggsave(filename = "selection_r0_r1.png", plot = plot,
       path ="/Users/raphaellemomal/these/R/images/", width = 11, height = 3)

vem=simr0$VEMr[[1]][[7]]
sum(F_Sym2Vec(vem$omega))/
  (ncol(vem$omega)*(ncol(vem$omega)-1)/2)
vem=simr0$list.vem0
sum(log(diag(vem$omega)))


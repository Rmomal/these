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

#-- function
# boot_VEM0<-function(Y,B=40, eps,alpha=1, maxIter){
#   list.res<-lapply(1:B, function(x){
#     n=nrow(Y); v=0.8; n.sample=round(0.8*n, 0)
#     ech=sample(1:n,n.sample,replace = FALSE)
#     Y.sample=Y[ech,]
#     PLNfit=PLN(Y.sample~1, control=list(trace=0))
#     sigma_obs=PLNfit$model_par$Sigma
#     MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; theta=PLNfit$model_par$Theta 
#     init0=initVEM(counts = Y.sample, initviasigma = NULL,  sigma_obs,r = 0)
#     Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 
#     VEM_r0<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
#                     Wg_init =Wginit,plot = FALSE, maxIter = maxIter,print.hist = FALSE,
#                     vraiOm = NULL, alpha=alpha, verbatim=FALSE, filterPg = TRUE, 
#                     filterWg = TRUE )
#     return(list(theta=theta, counts=Y.sample, VEM_r0=VEM_r0))
#   })
#   return(list.res)
# }
# boot_VEM0<-function(Y,MO,SO,B=40, eps,alpha=1, maxIter){
#   list.res<-lapply(1:B, function(x){
#     n=nrow(Y); v=0.8; n.sample=round(0.8*n, 0)
#     ech=sample(1:n,n.sample,replace = FALSE)
#     Y.sample=Y[ech,]
#     # bootstrap initialisation of omega, W and Wg
#     PLNfit=PLN(Y.sample~1, control=list(trace=0))
#     sigma_obs=PLNfit$model_par$Sigma
#     init0=initVEM(counts = Y.sample, initviasigma = NULL,  sigma_obs,r = 0)
#     Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 
#     # vem with original counts
#     VEM_r0<-VEMtree(Y, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
#                     Wg_init =Wginit,plot = FALSE, maxIter = maxIter,print.hist = FALSE,
#                     vraiOm = NULL, alpha=alpha, verbatim=FALSE, filterPg = TRUE, 
#                     filterWg = TRUE )
#     return(VEM_r0)
#   })
#   return(list.res)
# }
simu_vary_r<-function(p,n,r,B,rMax,maxIter,seed,type,alpha,eps=1e-3,nobeta=FALSE, cores=3, plot){
  # q=p+rMax
  # D=.Machine$double.xmax
  # alpha = (1/n)*((1/(q-1))*log(D) - 0.5*log(q*(q-1)))
  
  set.seed(seed)
  # sim data
  missing_data<-missing_from_scratch(n,p,r=r,type,plot)
  counts=missing_data$Y; omega=missing_data$Omega
  # Observed parameters
  PLNfit<-PLN(counts~1, control=list(trace=0))
  MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S ; 
  theta=PLNfit$model_par$Theta; sigma_obs=PLNfit$model_par$Sigma
  init0=initVEM(counts , initviasigma = NULL,  sigma_obs,r = 0)
  Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 
  # vem with original counts
  VEM_r0<-VEMtree(counts, MO, SO, MH=NULL,ome_init = omegainit,W_init =Winit,eps=eps,
                  Wg_init =Wginit,plot = FALSE, maxIter = maxIter,print.hist = FALSE,
                  alpha=alpha, verbatim=FALSE, filterWg = TRUE, filterDiag = FALSE, nobeta=nobeta )
  
  #  vary r
  VEMr<-lapply(1:rMax, function(r){
    cat(paste0(r," missing actors: "))
    t1<-Sys.time()
    cliques_spca <- boot_FitSparsePCA(scale(MO),B,r=r,cores=3)
    
    ListVEM<-List.VEM(cliquesObj =cliques_spca, counts, sigma_obs, MO,SO,r=r,alpha=alpha, cores=cores,
                      eps=eps,maxIter, nobeta=nobeta)
    t2<-Sys.time()
    runtime=difftime(t2,t1)
    cat(paste0(round(runtime,3), attr(runtime, "units"),"\n"))
    return(ListVEM)
  })
  saveRDS(list(list.vem0=VEM_r0, VEMr=VEMr, counts=counts, theta=theta),
          file=paste0("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/",type,"_seed",
                      seed,"_r",r,"_1-",rMax ,".rds") )       
}
# cliques_spca6 <- boot_FitSparsePCA(scale(counts),100,r=6)
# init=initVEM(counts = counts, initviasigma=cliques_spca6[[1]], sigma_obs,r = 6)
# Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
# VEM_6<-VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, 
#              eps=eps, alpha=1,verbatim = FALSE,
#              maxIter=maxIter, plot=TRUE,vraiOm=NULL, print.hist=FALSE, filterPg=FALSE,
#              filterWg = FALSE)
# lapply(cliques_spca6, function(clique){
#   length(unique(clique))==length(clique)
# })
# length((cliques_spca6))
# 
# `scale-free_seed1_r0_1-2`[[2]]$lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-V6) %>%
#   ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(V6)), size=3)+geom_line()+
#   facet_wrap(~key, scales="free")+
#   labs(x="iteration",y="", title="Lower bound and components")+mytheme+
#   scale_color_discrete("")

############
#-- run
seed=1; p=14 ; B=100; type="scale-free" ; n=200 
t1<-Sys.time()
simu_vary_r(seed=seed,B=B, n=n, p=p,r=0,maxIter=200,rMax=3, 
            type = type,nobeta=FALSE, plot=FALSE , cores=3, alpha=0.1)
t2<-Sys.time()
difftime(t2, t1)
t1<-Sys.time()
simu_vary_r(seed=seed,B=B, n=n, p=p,r=1,maxIter=200,
            rMax=3, type = type, nobeta=FALSE, plot=FALSE, cores=3, alpha=0.1)
t2<-Sys.time()
difftime(t2, t1)
#--
simr1<-readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/scale-free_seed1_r1_1-3.rds")        
simr0<-readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/scale-free_seed1_r0_1-3.rds")        

# - corriger la borne inf pour 0 acteurs manquants
Jcor_0r<-function(vem,p){
  O=1:p
  omega=vem$omega
  EhZZ=t(vem$M[,O])%*%vem$M[,O] + diag(colSums(vem$S[,O]))
  sigTilde =  sigTilde = (1/n)*EhZZ
  EsO=vem$Pg*vem$omega+diag(diag(vem$omega))
  EgO = nearPD(EsO, eig.tol=0.1)$mat
  JPLN_SigT = part_JPLN(sigTilde,EhZZ=EhZZ)
  JPLN_EgOm = part_JPLN(EgO,EhZZ=EhZZ, var=FALSE)
  diffJPLN = JPLN_SigT-JPLN_EgOm
  Jcor=tail(vem$lowbound$J,1)+diffJPLN
  return(Jcor)
}
Jcor_Delta<-function(vem,p){
  omega=vem$omega
  O=1:p ; H=(p+1):ncol(omega) ; r=length(H)
  EhZZ=t(vem$M[,O])%*%vem$M[,O] + diag(colSums(vem$S[,O]))
  sigTilde =  sigTilde = (1/n)*EhZZ
  EgO=vem$Pg*vem$omega+diag(diag(vem$omega))
  if(ncol(omega)==p){
    r=0
    EgOm = EgO[O,O]
  }else{
    if(r==1){
      EgOm = EgO[O,O] - matrix(EgO[O,H],p,r)%*%matrix(EgO[H,O],r,p)/EgO[H,H]
    }else{
      EgOm = EgO[O,O] - matrix(EgO[O,H],p,r)%*%solve(EgO[H,H])%*%matrix(EgO[H,O],r,p)
    }
  }
  EgOm = nearPD(EgOm, eig.tol=0.1)$mat
  JPLN_SigT = part_JPLN(sigTilde,EhZZ=EhZZ)
  JPLN_EgOm = part_JPLN(EgOm,EhZZ=EhZZ, var=FALSE)
  diffJPLN = JPLN_SigT-JPLN_EgOm
  Jcor=tail(vem$lowbound$J,1)+diffJPLN
  Delta = norm(solve(sigTilde) - EgOm,type = "F")
  return(c(Jcor=Jcor, Delta=Delta, r=r))
}

Jdata_r0<-do.call(rbind,c(list(Jcor_Delta(simr0$list.vem0,p=14)),
                          do.call(rbind, lapply(1:3, function(r){
                            lapply(simr0$VEMr[[r]],function(vem) Jcor_Delta(vem,p=14))}
                          ))))%>%
  as_tibble() %>% mutate(trueR = 0)
Jdata_r1<-do.call(rbind,c(list(Jcor_Delta(simr1$list.vem0,p=14)),
                          do.call(rbind, lapply(1:3, function(r){
                            lapply(simr1$VEMr[[r]],function(vem) Jcor_Delta(vem,p=14))}
                          )))) %>%
  as_tibble() %>% mutate(trueR = 1)

Jdata = rbind(Jdata_r0, Jdata_r1) %>% as_tibble()
plot=Jdata %>% group_by(trueR, r) %>% 
  mutate(q10Delta=quantile(Delta, 0.1)) %>% 
  filter(Delta<q10Delta) %>% mutate(maxJcor=max(Jcor)) %>% 
  ggplot(aes(as.factor(r), maxJcor, color=as.factor(trueR)))+
  geom_point(shape=16, size=3)+mytheme.dark("True r:")+
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


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
library(blockmodels)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-V2.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/modif_pkg.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/VEM_tools.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-exactDet.R")


mclust.init<-function(sigma_obs, trueClique, B,n){
  List.init=vector(mode = "list")
  p=ncol(sigma_obs)
  c=1
  for(i in 1:B){
    res<- init.mclust((cov2cor(sigma_obs)),title="Sigma",
                      trueClique = trueClique,n.noise=3*p,plot = FALSE)$init[[1]]
    if(length(res)>=1){
      List.init[[c]]<-list(res)
      c=c+1
    }
  }
  List.init<-unique(List.init)
  return(List.init)
}

varClust.init<-function(sigma_obs){
  vcPCA = F_VarClustPCA(sigma_obs)
  filter_content<-list()
  c=1
  lapply(seq_along(vcPCA$clustContent), function(x){
    if(length(vcPCA$clustContent[[x]])>1){
      filter_content[[c]]<<-list(vcPCA$clustContent[[x]])
      c<<-c+1
    } 
  })
  return(filter_content)
}
# vBIC computations
run.VEM<-function(clique,counts,sigma_obs,MO,SO,r){
  p=ncol(counts); O=1:p
  init=initVEM(counts = counts, initviasigma= (clique), sigma_obs=sigma_obs,MO=MO,r = r)
  omegainit=init$upsinit ; MHinit=init$MHinit; W_init=init$Winit ; Wg_init=init$Wginit
  #alpha<-computeAlpha(omegainit[O,O], MO, SO, plot=FALSE)
  vem<-tryCatch({VEMtree(counts,MO,SO,MH=MHinit,omegainit,W_init,Wg_init, eps=1e-3,
                         alpha=0.1, maxIter=100, plot=FALSE,verbatim=FALSE,trackJ=FALSE)}, 
                error=function(e){e},
                finally={}) 
  return(vem)
}

vec.vBIC<-function(List.vem,counts,theta, matcovar,r){
  
  n=nrow(counts);p= ncol(counts)
  unlist(lapply(List.vem,function(vem){
    J<-True_lowBound(counts,vem$M,vem$S, theta, matcovar,vem$W, vem$Wg, vem$Pg, vem$omega )
    vBIC<-VBIC(J,p,r=r, d=ncol(matcovar), n=n)
    return(vBIC)
  }))
}

run.method<-function(cliques,method,counts,sigma_obs,MO,SO, MHinit,ome,r, trueClique ){
  p=ncol(counts)
  vem=mclapply(cliques, function(init){
    run.VEM(clique = init,counts=counts, sigma_obs,MO=MO, SO=SO, r=r)
  }, mc.cores=3)
  vec.perf<-data.frame(t(sapply(vem, function(x){
    if(length(x)>5){
      auc=round(auc(pred = x$Pg, label = ome),4)
      h=(p+1):(p+r)
      vec=accppvtpr(x$Pg,ome,h=h,seuil=0.5)
      PPVH=vec[5]
      TPRH=vec[8]
      FPRH=vec[10]
      FNRH=vec[11]
      J=  tail(x$lowbound$J,1)
      res=c(auc=auc, PPVH=PPVH,TPRH=TPRH,FPRH=FPRH,FNRH=FNRH,J=J)
    }else{ res=c(auc=NaN, PPVH=NaN, TPRH=NaN,FPRH=NaN,FNRH=NaN, J=NaN) }
    return(res)
  })))
  data=computeFPN(cliques, trueClique[[1]], p) %>% 
    mutate(crit=FP+FN, crit.rank=rank(crit, ties.method = "min"), J=vec.perf$J, 
           sizes=unlist(lapply(cliques, function(x) length(x[[1]]))),
            auc=vec.perf$auc, PPVH=vec.perf$PPVH, TPRH=vec.perf$TPRH, 
           FPRH=vec.perf$FPRH, FNRH=vec.perf$FNRH,
           nbinit=length(cliques),method=c(method))
  # best=which.max(vec.perf$J)
  # choice=data[best,]
  
  return(data)
}
SelectInitClique<-function(seed, B, n=200, p=14, r=1){
  t1<-Sys.time()
  cat(paste0("\n seed ",seed, "\n"))
  set.seed(seed)
  type="scale-free" ;O=1:p;plot=FALSE
  #browser()
  missing_data<-missing_from_scratch(n,p,r,type,plot)
  counts=missing_data$Y;  omega=missing_data$Upsilon;
  trueClique=missing_data$TC; h=(p+1):(p+r)
  PLNfit<-PLN(counts~1, control=list(trace=0))
  MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S ; sigma_obs=PLNfit$model_par$Sigma
  ome=missing_data$G
  #-- normalize the PLN outputs
  D=diag(sigma_obs)
  matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
  MO=MO*matsig
  SO=SO*matsig^2
  
  #true perfs
  vem=run.VEM(clique = trueClique ,counts=counts, sigma_obs, MO=MO, SO=SO,r=r)
  if(length(vem)>5){
    vec=accppvtpr(vem$Pg,ome,h=h,seuil=0.5)
    oracle.perf<-data.frame(auc=round(auc(pred =vem$Pg, label = ome),4),
                            PPVH=vec[5],  TPRH=vec[8], FPRH=vec[10],FNRH=vec[11],
                            J= tail(vem$lowbound$J,1))
  }else{oracle.perf<-data.frame(auc=NaN,
                                PPVH=NaN,TPRH=NaN,FPRH=NaN,FNRH=NaN,
                                J= NaN)}


  #=== MCLUST
  # cat(paste0("\n ------ MCLUST ------ \n"))
  # cliques_mclust<-mclust.init(sigma_obs, trueClique[[1]],B, n)
  # goodinit=which(do.call(rbind , lapply(cliques_mclust, function(x) length(x[[1]])))>1)
  # cliques_mclust2=cliques_mclust[goodinit]
  # cat(paste0("   Fitting ",length(cliques_mclust2)," VEM..."))
  # vem_mclust=mclapply(cliques_mclust2, function(init){
  #   run.VEM(clique = (init),counts=counts,sigma_obs, MO=MO, SO=SO, r=r)
  # }, mc.cores=3)
  # 
  # vec.perf<-data.frame(t(sapply(vem_mclust, function(x){
  #   if(length(x)>5){
  #     auc=round(auc(pred = x$Pg, label = ome),4)
  #     h=(p+1):(p+r)
  #     PPVH=accppvtpr(x$Pg,ome,h=h,seuil=0.5)[5]
  #     TPRH=accppvtpr(x$Pg,ome,h=h,seuil=0.5)[8]
  #     J=tail(x$lowbound$J,1)
  #     res=c(auc=auc,PPVH= PPVH,TPRH= TPRH,J=J)
  #   }else{ res=c(auc=NaN,PPVH=NaN,TPRH= NaN,J=NaN)}
  #   return(res)
  # })))
  # 
  # FPN_mclust<-computeFPN(cliques_mclust2, trueClique[[1]], p)
  # mclust.choice=FPN_mclust %>% mutate(crit=FP+FN, crit.rank=rank(crit, ties.method = "min"),
  #                            J=vec.perf$J,sizes=unlist(lapply(cliques_mclust2, function(x) length(x[[1]]))),
  #                            auc=vec.perf$auc, PPVH=vec.perf$PPVH,TPRH=vec.perf$TPRH,
  #                            nbinit=length(cliques_mclust2),  method="mclust")
  # mclust.choice<-data[which.max(data$J),-c(3, 4)]
  # mclust.choice$method="mclust"
  # data3<-data %>% filter(sizes>2) %>% dplyr::select(-crit, -crit.rank)
  #
  # if(nrow(data3)!=0){
  #   mclust.choice.min3<- data3%>% filter(J==max(J, na.rm=TRUE))
  #   mclust.choice.min3$method="mclust.min3"
  #   mclust.choice<-rbind(mclust.choice,mclust.choice.min3)
  # }

  #=== sPCA
  # cat(paste0("\n ------ sPCA ------ \n"))
  # cat(paste0("\n bootstrap \n"))
  # cliques_boot_spca <- boot_FitSparsePCA(scale(MO),B,r=r,minV = 1)$cliqueList
  # cat(paste0("   Fitting ",length(cliques_boot_spca)," VEM..."))
  # boot.spca.choice=run.method(cliques_boot_spca,"boot.sPCA",counts,sigma_obs,MO,SO,MHinit,
  #                             ome,r,trueClique )

  #---
  cat(paste0("\n two axes sPCA \n"))
  cliques_spca<-FitSparsePCA(scale(MO), r=2,minV = 1)$cliques
  complement=lapply(cliques_spca, function(clique){setdiff(1:p,clique)})
  cliques_spca=lapply(c(cliques_spca,complement), function(cl) list(cl))
  cat(paste0("   Fitting ",length(cliques_spca)," VEM..."))
  spca.choice=run.method(cliques_spca,"sPCA",counts,sigma_obs,MO,SO, MHinit,ome,r,
                         trueClique )

  #=== varClust
  # cat(paste0("\n ------ VarClust ------ \n"))
  # cliques_varclust<-varClust.init(sigma_obs)
  # cat(paste0("   Fitting ",length(cliques_varclust)," VEM..."))
  # varclust.choice<-run.method(cliques_varclust,"VarClust",counts,sigma_obs,MO,SO,MHinit,
  #                             ome,r,trueClique )
  # 
  # #---
  # #=== Blockmodels
  # cat(paste0("\n ------ BlockModels ------ \n"))
  # cliques_bm<-tryCatch({init_blockmodels(k=3, counts, sigma_obs, MO, SO)$cliqueList},
  #                      error=function(e){e},finally={})  
  # #if(typeof(cliques_bm[[1]])!="character"){
  #   goodinit=which(do.call(rbind , lapply(cliques_bm, function(x) length(x[[1]])))>1)
  #   cliques_bm2=cliques_bm[goodinit]
  #   cat(paste0("   Fitting ",length(cliques_bm2)," VEM..."))
  #   bm.choice<-run.method(cliques_bm2,"BM",counts,sigma_obs,MO,SO,MHinit,ome,r,
  #                         trueClique )
  # }else{
  #   bm.choice<- data.frame(FP=NaN,FN=NaN, J=NaN,sizes=NaN, 
  #                          auc=NaN,PPVH=NaN,nbinit=NaN,  method="BM")
  #}
  
  
  #---
  # choices<-rbind(boot.spca.choice,spca.choice,bm.choice, mclust.choice , varclust.choice,
  #                data.frame(FP=0,FN=0,crit=0, crit.rank=1, J=oracle.perf$J,sizes=length(trueClique[[1]]), 
  #                           auc=oracle.perf$auc,PPVH=oracle.perf$PPVH,TPRH=oracle.perf$TPRH,nbinit=1,  method="Oracle"))
  choices<-rbind(spca.choice, data.frame(FP=0,FN=0,crit=0, crit.rank=1, 
                                         J=oracle.perf$J,sizes=length(trueClique[[1]]),
                                         auc=oracle.perf$auc,PPVH=oracle.perf$PPVH,
                                         TPRH=oracle.perf$TPRH,FPRH=oracle.perf$FPRH,
                                         FNRH=oracle.perf$FNRH,nbinit=1,  method="Oracle"))
  choices$seed=seed
  t2<-Sys.time()
  runtime=difftime(t2, t1)
  saveRDS(choices, paste0("/Users/raphaellemomal/simulations/selecInit/sPCA_choices_seed",seed,".rds"))
  cat(paste0("\nseed ", seed," in ",round(runtime,3), attr(runtime, "units"),"\n"))
  return(list(choice=choices))
}

seed=1:200
t1<-Sys.time()
res2<-lapply(seed, function(x){ SelectInitClique(x, B=100)})
t2<-Sys.time()
difftime(t2, t1)

saveRDS(object = res2, file = "/Users/raphaellemomal/simulations/Simu/selecInit_sPCA.rds" )
# rank.size<-unlist(lapply(res, function(x){rank(x$perf$clique.size)}))
# rank.crit<-unlist(lapply(res, function(x){ x$perf$crit.rank}))
# plot(rank.crit, rank.size)
# quality<-unlist(lapply(res, function(x){x$quality}))
# quality %>% as_tibble() %>% ggplot(aes(value, y=..prop..))+geom_bar(width=0.5, alpha=0.6, fill="deepskyblue3")+mytheme+labs(x="rank of FN+FP", y="choice with max vBIC")

quality<-do.call(rbind,lapply(res, function(x){x$choice})) %>% as_tibble()


test=quality %>% dplyr::select(-FN, -FP, -auc,-J, -PPVH) %>% 
  spread(method, sizes) %>%dplyr::select(Oracle) 
truthSizes=sort(unique(test$Truth), decreasing = FALSE)
lapply(truthSizes, function(x){
  which(test$Truth==x)
})
#filter by the number of original neighbors
index=which(test$Truth>2)
# FPN
p1<-quality %>% filter(method!="Oracle") %>%  
  gather(key, value, -method, -sizes, -J, -seed,-nbinit) %>%  
  mutate(key = fct_relevel(key,
                             "auc","PPVH","FN","FP")) %>%
  ggplot(aes(as.factor(method), value, color=as.factor(key),fill=as.factor(key)))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  geom_boxplot(alpha=0.3, width=0.6) +theme_light()+labs(x="",y="",title="On 100 scale-free graphs")+
  scale_color_hp_d("",option = "Ravenclaw")+
  scale_fill_hp_d("",option = "Ravenclaw")
# AUC PPVH
quality %>% filter(method!="Oracle") %>% gather(key, value, -method, -sizes, -J,-nbinit, -seed,-FN, -FP) %>%
  ggplot(aes(as.factor(method), value, color=as.factor(key),fill=as.factor(key)))+
  geom_boxplot(alpha=0.3) +theme_light()+labs(x="",y="",title="On 200 scale-free graphs")+
  scale_color_hp_d("",option = "Ravenclaw")+
  scale_fill_hp_d("",option = "Ravenclaw")
test=quality %>% filter(method!="Oracle")%>% dplyr::select(-FN, -FP, -sizes, -J,-nbinit, -PPVH) %>% spread(method, auc)
AUCranks<-t(apply(test[,-1],1,function(x){
  rank(x, ties.method = "max")
}))
AUCranks %>%as_tibble() %>%  gather(key, value) %>% 
  ggplot(aes( x=value, y=..prop..,color=key, fill=key ))+
  geom_bar(position = position_dodge(),width = 0.7)+#stat="binline",bins=14, scale=0.8, draw_baseline=FALSE)
  theme_light()+labs(x="",y="",title="AUC ranks (max is best)")+
  scale_color_hp_d("",option = "Ravenclaw")+ scale_fill_hp_d("",option = "Ravenclaw")

test=quality %>% filter(method!="Oracle")%>% dplyr::select(-FN, -FP, -sizes, -J,-nbinit, -auc) %>% spread(method, PPVH)
PPVHranks<-t(apply(test[,-1],1,function(x){
  rank(x, ties.method = "max")
}))
PPVHranks %>%as_tibble() %>%  gather(key, value) %>% 
  ggplot(aes( x=value, y=..prop..,color=key, fill=key ))+
  geom_bar(position = position_dodge(),width = 0.7)+#stat="binline",bins=14, scale=0.8, draw_baseline=FALSE)
  theme_light()+labs(x="",y="",title="PPVH ranks (max is best)")+
  scale_color_hp_d("",option = "Ravenclaw")+ scale_fill_hp_d("",option = "Ravenclaw")

# Sizes
quality %>%  
  ggplot(aes( sizes, y=..prop.., color=method, fill=method ))+
  geom_bar() +theme_light()+labs(x="",y="",title="Size of found cliques compared to truth")+
  scale_color_hp_d("",option = "Ravenclaw")+scale_fill_hp_d("",option = "Ravenclaw")

p2<-quality %>%  
  ggplot(aes( x=as.numeric(sizes), y=method,color=method, fill=method ))+
  geom_density_ridges()+#stat="binline",bins=14, scale=0.8, draw_baseline=FALSE)
  theme_light()+labs(x="",y="",title="Size of found cliques compared to truth")+
  scale_color_hp_d("",option = "Ravenclaw")+ scale_fill_hp_d("",option = "Ravenclaw")

quality %>% filter(method!="Oracle") %>% 
  ggplot(aes( x=auc, y=method,color=method, fill=method ))+
  geom_density_ridges()+#stat="binline",bins=14, scale=0.8, draw_baseline=FALSE)
  theme_light()+labs(x="",y="",title="PPVH")+
  scale_color_hp_d("",option = "Ravenclaw")+ scale_fill_hp_d("",option = "Ravenclaw")

# Gap in vBIC
test=quality %>% dplyr::select(-FN, -FP, -sizes,-nbinit, -auc, -PPVH) %>% spread(method, J)
Jdiff<-apply(test[,-c(1)],2,function(x){
  x-test$Oracle
})
p3<-Jdiff %>%as_tibble() %>% dplyr::select(-Oracle)%>%  gather(key, value) %>% 
  ggplot(aes( x=key, y=value,color=key, fill=key ))+
  geom_boxplot(width=0.4)+#stat="binline",bins=14, scale=0.8, draw_baseline=FALSE)
  theme_light()+labs(x="",y="",title="J - Oracle.J")+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_hp_d("",option = "Ravenclaw")+ scale_fill_hp_d("",option = "Ravenclaw")+ coord_flip()

# Best vBic in rank
Jranks<-t(apply(test %>% dplyr::select(-seed, -Oracle),1,function(x){
  rank(x, ties.method = "max")
}))
Jranks %>%as_tibble() %>%  gather(key, value) %>% 
  ggplot(aes( x=value, y=..prop..,color=key, fill=key ))+
  geom_bar(position = position_dodge(),width = 0.7)+#stat="binline",bins=14, scale=0.8, draw_baseline=FALSE)
  theme_light()+labs(x="",y="",title="J ranks (max is best)")+
  scale_color_hp_d("",option = "Ravenclaw")+ scale_fill_hp_d("",option = "Ravenclaw")


plot=grid.arrange(p1, p2, p3, ncol=1)
ggsave(filename = "selectInit_final.png", plot = p1, path ="/Users/raphaellemomal/these/R/images/", width = 7, height = 4 )
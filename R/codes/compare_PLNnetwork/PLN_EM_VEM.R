# compare PLN EMtree and VEMtree is smth I wanna do for a long time.
# LEt's see what happens after 1h coding this.
# start 10:10

library(tidygraph)
library(ggraph)
library(tidyverse)
library(EMtree)
library(Matrix)
library(mvtnorm)
library(reshape2)
library(PLNmodels)
library(gridExtra)
library(parallel)
##############
# RUN

sapply(c("erdos","cluster"), function(type){
  type=type
  sapply(1:2, function(numdens){
    cat(paste0("Type ",type,", density ", numdens,"..."))
    t1<-Sys.time()
    lapply(1:100, function(seed){
      file=paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/simu_seed",seed,"_",type,"_dens",numdens,".rds")
      if(!file.exists(file)){
        set.seed(seed) ; n=200 ; p=15  ;O=1:p
        dens=c(3/p,5/p)[numdens]
        missing_data<-missing_from_scratch(n,p,r=0,type,plot=FALSE, dens)
        counts=missing_data$Y
        sigmaO= missing_data$Sigma
        upsilon=missing_data$Upsilon
        PLNfit<-PLN(counts~1, control = list(trace=0))
        MO<-PLNfit$var_par$M  
        SO<-PLNfit$var_par$S  
        sigma_obs=PLNfit$model_par$Sigma
        G=(missing_data$G); 
        
        #-- normalize the PLN outputs
        D=diag(sigma_obs)
        matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
        MO=MO*matsig
        SO=SO*matsig^2
        #RUN VEMtree
        init=initVEM(counts = counts,initviasigma=NULL, cov2cor(sigma_obs),MO,r = 0) 
        Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
        VEMfit<- tryCatch({VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.1,
                                   maxIter=100, plot=FALSE,print.hist=FALSE, verbatim = FALSE,trackJ=FALSE)},
                          error=function(e){e}, finally={})
        
        # RUN PLNnetwork
        network_models<-PLNnetwork(counts~1, control_init = list(nPenalties=200), control_main =list(trace=0, cores=2))
        PLNnet.path=lapply(network_models$penalties, function(pen){
          model_pen <- getModel(network_models, pen) 
          pen*(model_pen$latent_network(type="precision")!=0)  
        })
        K.score <- Reduce("+",PLNnet.path)
        diag(K.score)=0
        PLNscores<- as.matrix(K.score / max(K.score))
        
        # RUN EMtree
        EMfit=EMtree(PLNfit,verbatim = FALSE)
        saveRDS(list(G=G, PLN=PLNscores,EM= EMfit, VEM=VEMfit),
                file = paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/simu_seed",seed,"_",type,"_dens",numdens,".rds"))
      } 
    })
    t2<-Sys.time()
    time=difftime(t2, t1)
    cat(paste0(round(time,3), attr(time, "units"),"\n"))
  })
})


################
### focus on one dataset
seed=6 ; type="erdos" ; numdens=1
listFit=readRDS(paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/simu_seed",seed,"_",type,"_dens",numdens,".rds"))
PLNscores=listFit$PLN
EMprob=listFit$EM$edges_prob
VEMprob=listFit$VEM$Pg
G=listFit$G
# plot image scores and auc
g1<-ggimage(G)+labs(title="Truth")
g2<-draw_network(G, curv=0, layout="kk")$G
g3<-ggimage(PLNscores)+labs(title=paste0("PLNnetwork: AUC=",auc(pred =PLNscores,label = G )))
g4<-ggimage(EMprob)+labs(title=paste0("EMtree: AUC=",auc(pred =EMprob,label = G )))
g5<-ggimage(VEMprob)+labs(title=paste0("VEMtree: AUC=",auc(pred =VEMprob,label = G )))
my_layout <-matrix(c(6,1,1,2,2,2,6,1,1,2,2,2,3,3,4,4,5,5,3,3,4,4,5,5), ncol=6, nrow=4, byrow = TRUE)
grid.arrange(g1, g2, g3, g4, g5, ncol=2,layout_matrix=my_layout)

# plots precrec
list_values=lapply(1:3, function(x){
  method=c("PLNnetwork","EMtree","VEMtree")[x]
  scores=list(PLNscores, EMprob,VEMprob)[[x]]
  values=courbes_seuil(probs = scores,omega = G, h = 15,seq_seuil = seq(0,max(scores),max(scores)/50)) %>% 
    dplyr::select(seuil, PPV, TPR) %>% mutate(method=method)
})
values=do.call(rbind, list_values)
values %>%as_tibble() %>%  
  ggplot(aes(TPR,PPV,color=method))+
  geom_rect(aes(xmin=0.5, xmax=1, ymin=0.5, ymax=1), fill="gray90", color="gray90",alpha=0.1)+
  geom_point()+  geom_line()+facet_wrap(~method)+mytheme.dark("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+guides(color=FALSE)

### pooled results
pooled=do.call(rbind,lapply(c("erdos","cluster"), function(type){
  do.call(rbind,lapply(1:2, function(numdens){
    do.call(rbind, lapply(1:100, function(seed){
      listFit=readRDS(paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/simu_seed",seed,"_",type,"_dens",numdens,".rds"))
      PLNscores=listFit$PLN
      if(sum(!is.na(PLNscores))==0) PLNscores = matrix(0, 15, 15)
      EMprob=listFit$EM$edges_prob
      if(length(listFit$VEM)>5){
        VEMprob=listFit$VEM$Pg
        G=listFit$G
        AUC=data.frame(AUC=c(auc(pred =PLNscores,label = G ),
                             auc(pred =EMprob,label = G ),
                             auc(pred =VEMprob,label = G )),method=c("PLNnetwork","EMtree","VEMtree"))
        
        values=rbind(seq_PPV_TPR(probs = PLNscores, omega=G) %>% mutate(method="PLNnetwork"),
                     seq_PPV_TPR(probs = EMprob, omega=G) %>% mutate(method="EMtree"),
                     seq_PPV_TPR(probs = VEMprob, omega=G) %>% mutate(method="VEMtree"))
        seuils= values %>% group_by(method) %>%filter(PPV>=0.5, TPR>=0.5) %>% 
          summarize(minseuil=min(seuil), maxseuil=max(seuil)) %>% as_tibble() 
        
        max_insquare=values %>% group_by(method) %>%filter(PPV>=0.5, TPR>=0.5) %>% 
          mutate(sumcrit=PPV+TPR) %>%  filter(sumcrit==max(sumcrit, na.rm=TRUE)) %>% 
          summarize(maxPPV=max(PPV, na.rm=TRUE), maxTPR=max(TPR, na.rm=TRUE)) %>% as_tibble()
        
        max_global=values %>% group_by(method) %>% 
          mutate(sumcrit=PPV+TPR) %>%  filter(sumcrit==max(sumcrit, na.rm=TRUE)) %>% 
          summarize(glob.maxPPV=max(PPV, na.rm=TRUE), glob.maxTPR=max(TPR, na.rm=TRUE)) %>% as_tibble()
        maxi=left_join(max_global,max_insquare, by="method")
        values=left_join( maxi,seuils, by="method")  
        values=left_join(AUC,values, by="method") %>%  mutate( seed=seed, numdens=numdens, type=type)
        
      }
      else{
        values= data.frame(AUC=NA,  method=c("PLNnetwork","EMtree","VEMtree"),glob.maxPPV=NA, glob.maxTPR=NA,
                           maxPPV=NA, maxTPR=NA, minseuil=NA, maxseuil=NA, seed=seed, numdens=numdens, type=type)
      }
      
      return(values)
    }))
  }))
}))
length(unique(pooled$seed[which(is.na(pooled$AUC))])) # 29 crashed for VEMtree
pooled %>% as_tibble() %>%
  mutate(method=factor(method, levels=c("PLNnetwork","EMtree","VEMtree"))) %>% 
  gather(key, value, -seed, -numdens, -type,-method,-minseuil,-maxseuil) %>% 
  ggplot(aes(type, value,color=as.factor(numdens), fill=as.factor(numdens)))+
  geom_boxplot(width=0.4, notch = FALSE, alpha=0.6)+
  facet_grid(key~method)+mytheme.dark("")+labs(x="")
pooled %>% as_tibble() %>%
  mutate(method=factor(method, levels=c("PLNnetwork","EMtree","VEMtree"))) %>% 
  ggplot(aes(type, AUC,color=as.factor(numdens), fill=as.factor(numdens)))+
  geom_boxplot(width=0.4, notch = FALSE, alpha=0.6)+
  facet_grid(~method)+mytheme.dark("")+labs(x="")

pooled %>% as_tibble() %>% dplyr::select(-maxseuil, -minseuil) %>% 
  ggplot(aes(maxTPR, maxPPV, color=method))+
  geom_rect(aes(xmin=0.5, xmax=1, ymin=0.5, ymax=1), fill="gray90", color="gray90",alpha=0.1)+
  geom_point()+ facet_grid(numdens~type)+mytheme.dark("")

pooled %>% as_tibble() %>% mutate(rangeSeuil=maxseuil-minseuil) %>% 
  ggplot(aes(rangeSeuil, method, color=method, fill=method))+guides(color=FALSE, fill=FALSE)+
  geom_density_ridges(alpha=0.6)+facet_grid(numdens~type)+mytheme.dark("")+labs(y="")
###########
prec_rec=do.call(rbind,lapply(c("erdos","cluster"), function(type){
  do.call(rbind,lapply(1:2, function(numdens){
    do.call(rbind, lapply(1:100, function(seed){
      listFit=readRDS(paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/simu_seed",seed,"_",type,"_dens",numdens,".rds"))
      PLNscores=listFit$PLN
      if(sum(!is.na(PLNscores))==0) PLNscores = matrix(0, 15, 15)
      EMprob=listFit$EM$edges_prob
      if(length(listFit$VEM)>5){
        VEMprob=listFit$VEM$Pg
        G=listFit$G
        values=rbind(seq_PPV_TPR(probs = PLNscores, omega=G,seq_seuil=seq(0,max(PLNscores),max(PLNscores)/200)) %>% mutate(method="PLNnetwork"),
                     seq_PPV_TPR(probs = EMprob, omega=G,seq_seuil=seq(0,max(EMprob),max(EMprob)/200)) %>% mutate(method="EMtree"),
                     seq_PPV_TPR(probs = VEMprob, omega=G,seq_seuil=seq(0,max(VEMprob),max(VEMprob)/200)) %>% mutate(method="VEMtree"))  %>%  mutate( seed=seed, numdens=numdens, type=type)
        return(values)
      }}))
  }))
}))  
g1<-prec_rec %>% as_tibble() %>% filter(method=="PLNnetwork") %>% group_by(type, numdens) %>% arrange(seuil, by_group=TRUE) %>% 
  ggplot(aes(TPR, PPV))+geom_hline(yintercept=0.5, linetype="dashed", color="gray")+geom_vline(xintercept=0.5, linetype="dashed", color="gray")+
  geom_point(alpha=0.2, size=1)+geom_line()+labs(title="PLNnetwork")+
  facet_grid(type~numdens)+mytheme.dark("")
g2<-prec_rec %>% as_tibble() %>% filter(method=="EMtree") %>% group_by(type, numdens) %>% arrange(seuil, by_group=TRUE) %>% 
  ggplot(aes(TPR, PPV))+geom_hline(yintercept=0.5, linetype="dashed", color="gray")+geom_vline(xintercept=0.5, linetype="dashed", color="gray")+
  geom_point(alpha=0.2, size=1)+geom_line()+facet_grid(type~numdens)+mytheme.dark("")+labs(title="EMtree")
g3<-prec_rec %>% as_tibble() %>% filter(method=="VEMtree")%>% group_by(type, numdens) %>% arrange(seuil, by_group=TRUE) %>% 
  ggplot(aes(TPR, PPV))+geom_hline(yintercept=0.5, linetype="dashed", color="gray")+geom_vline(xintercept=0.5, linetype="dashed", color="gray")+
  geom_point(alpha=0.2, size=1)+geom_line()+facet_grid(type~numdens)+mytheme.dark("")+labs(title="VEMtree")
grid.arrange(g1, g2, g3, ncol=3)
# pour PLNnetwork
model_StARS <- getBestModel(network_models, "StARS")
my_graph <- plot(model_StARS, plot = FALSE)
ggimage(as.matrix(as_adj(my_graph)))
test
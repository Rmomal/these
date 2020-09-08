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
library(ROCR)
##############
# RUN
pal<-c("#FA7E1E","#D62976","#008080")#,"#358CC3"
mytheme.dark <-function(legend){list= list(theme_light(), 
                                           theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                                                 plot.title = element_text(hjust = 0.5)),
                                           scale_color_manual(legend, values=pal),
                                           scale_fill_manual(legend, values=pal),
                                           scale_shape(legend))}
EMtree_optimal<-function(Y,m,cores,difficulty){
  cond.tol<-switch(as.character(difficulty),"1"=1e-12,"2"=1e-6)
  resample<-ResampleEMtree(counts=Y, covar_matrix=m, S=20, maxIter=50,
                           cond.tol=cond.tol,cores=cores)
  pmat<-resample$Pmat # rect. matrix gathering edges probabilities obtained for all sub-samples
  Freqs=freq_selec(pmat, Pt=2/ncol(Y))
  inf<-1*( Freqs>= 0.8) # optimal final network after thresholding to 80%
  return(list(freq=Freq, net=inf))
}
##############
sapply(c("erdos","cluster"), function(type){
  type=type
  sapply(1:2, function(numdens){
    cat(paste0("Type ",type,", density ", numdens,"..."))
    t1<-Sys.time()
    lapply(1:100, function(seed){
      file=paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/simu_seed",seed,"_",type,"_dens",numdens,".rds")
     # if(!file.exists(file)){
        set.seed(seed) ; n=200 ; p=15  ;O=1:p
        dens=c(3/p,5/p)[numdens]
        missing_data<-missing_from_scratch(n,p,r=0,type,plot=FALSE, dens)
        counts=missing_data$Y
        # sigmaO= missing_data$Sigma
        # upsilon=missing_data$Upsilon
        # PLNfit<-PLN(counts~1, control = list(trace=0))
        # MO<-PLNfit$var_par$M
        # SO<-PLNfit$var_par$S
        # sigma_obs=PLNfit$model_par$Sigma
        # G=(missing_data$G);
        # 
        # #-- normalize the PLN outputs
        # D=diag(sigma_obs)
        # matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
        # MO=MO*matsig
        # SO=SO*matsig^2
        #RUN VEMtree
        # init=initVEM(counts = counts,initviasigma=NULL, cov2cor(sigma_obs),MO,r = 0)
        # Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
        # VEMfit<- tryCatch({VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.1,
        #                            maxIter=100, plot=FALSE,print.hist=FALSE, verbatim = FALSE,trackJ=FALSE)},
        #                   error=function(e){e}, finally={})

        # RUN PLNnetwork
        network_models<-PLNnetwork(counts~1, control_init = list(min.ratio=1e-3), control_main =list(trace=0, cores=2))
        # PLNnet.path=lapply(network_models$penalties, function(pen){
        #   model_pen <- getModel(network_models, pen)
        #   pen*(model_pen$latent_network(type="precision")!=0)
        # })
        # K.score <- Reduce("+",PLNnet.path)
        # diag(K.score)=0
        # PLNscores<- as.matrix(K.score / max(K.score))
        model_StARS <-1*(getBestModel(network_models, "StARS")$latent_network(type="precision")!=0)
        
       # plot(network_models)
        #RUN EMtree
      #  EMfit=EMtree(PLNfit,verbatim = FALSE)
        # saveRDS(list(G=G, PLN=PLNscores,EM= EMfit, VEM=VEMfit),
        #         file = paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/simu_seed",seed,"_",type,"_dens",numdens,".rds"))
      #ResamEMfit=EMtree_optimal(Y=counts,m=NULL,cores =3,difficulty = numdens )
      saveRDS(list(PLNstars=model_StARS),
                      file = paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/PLNstars_simu_seed",
                                    seed,"_",type,"_dens",numdens,".rds"))

        #} 
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
      ResampEMFreq=readRDS(paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/ResampEM_simu_seed",seed,"_",type,"_dens",numdens,".rds"))$ResampEM
      PLNscores=listFit$PLN
      if(sum(!is.na(PLNscores))==0) PLNscores = matrix(0, 15, 15)
      EMprob=listFit$EM$edges_prob
      if(length(listFit$VEM)>5){
        VEMprob=listFit$VEM$Pg
        G=listFit$G
        res=data.frame(F_Sym2Vec(VEMprob),F_Sym2Vec(EMprob),F_Sym2Vec(PLNscores),F_Sym2Vec(G), type, numdens, seed)
        colnames(res)=c("nestor","EMtree","PLN-network","G","type","numdens","seed")
   
        # AUC=data.frame(AUC=c(auc(pred =PLNscores,label = G ),
        #                      auc(pred =EMprob,label = G ),
        #                      auc(pred =ResampEMFreq,label = G ),
        #                      auc(pred =VEMprob,label = G )),method=c("PLNnetwork","EMtree","ResampEM","VEMtree"))
        # 
        # values=rbind(seq_PPV_TPR(probs = PLNscores, omega=G) %>% mutate(method="PLNnetwork"),
        #              seq_PPV_TPR(probs = EMprob, omega=G) %>% mutate(method="EMtree"),
        #              seq_PPV_TPR(probs = ResampEMFreq, omega=G) %>% mutate(method="ResampEM"),
        #              seq_PPV_TPR(probs = VEMprob, omega=G) %>% mutate(method="VEMtree"))
        # seuils= values %>% group_by(method) %>%filter(PPV>=0.5, TPR>=0.5) %>% 
        #   summarize(minseuil=min(seuil), maxseuil=max(seuil)) %>% as_tibble() 
        # 
        # max_insquare=values %>% group_by(method) %>%filter(PPV>=0.5, TPR>=0.5) %>% 
        #   mutate(sumcrit=PPV+TPR) %>%  filter(sumcrit==max(sumcrit, na.rm=TRUE)) %>% 
        #   summarize(maxPPV=max(PPV, na.rm=TRUE), maxTPR=max(TPR, na.rm=TRUE)) %>% as_tibble()
        # 
        # max_global=values %>% group_by(method) %>% 
        #   mutate(sumcrit=PPV+TPR) %>%  filter(sumcrit==max(sumcrit, na.rm=TRUE)) %>% 
        #   summarize(glob.maxPPV=max(PPV, na.rm=TRUE), glob.maxTPR=max(TPR, na.rm=TRUE)) %>% as_tibble()
        # maxi=left_join(max_global,max_insquare, by="method")
        # values=left_join( maxi,seuils, by="method")  
        # values=left_join(AUC,values, by="method") %>%  mutate( seed=seed, numdens=numdens, type=type)
        #
      }
      else{
        # values= data.frame(AUC=NA,  method=c("PLNnetwork","EMtree","ResampEM","VEMtree"),glob.maxPPV=NA, glob.maxTPR=NA,
        #                    maxPPV=NA, maxTPR=NA, minseuil=NA, maxseuil=NA, seed=seed, numdens=numdens, type=type)
        res=data.frame(NA, NA, NA, NA, type, numdens, seed)
        colnames(res)=c("nestor","EMtree","PLN-network","G","type","numdens","seed")
         }
      
     # return(values)
      return(res)
    }))
  }))
})) %>% as_tibble()
test=pooled %>% filter(is.na(AUC), method=="VEMtree") %>% dplyr::select(numdens, type)
length(unique(pooled[which(is.na(pooled$AUC)),])) # 29 crashed for VEMtree
pooled %>% as_tibble() %>%
  mutate(method=factor(method, levels=c("PLNnetwork","EMtree","VEMtree"))) %>% 
  gather(key, value, -seed, -numdens, -type,-method,-minseuil,-maxseuil) %>% 
  ggplot(aes(type, value,color=as.factor(numdens), fill=as.factor(numdens)))+
  geom_boxplot(width=0.4, notch = FALSE, alpha=0.6)+
  facet_grid(key~method)+mytheme.dark("")+labs(x="")
pooled %>% as_tibble() %>% mutate(dens=ifelse(numdens==1,"3/p","5/p")) %>%
  mutate(method=factor(method, levels=c("PLNnetwork","EMtree","VEMtree"))) %>% 
  ggplot(aes(type, AUC,color=as.factor(dens), fill=as.factor(dens)))+
  geom_boxplot(width=0.4, notch = TRUE, alpha=0.8)+
  facet_grid(~method)+mytheme.dark("Density:")+labs(x="")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

g=pooled %>% as_tibble() %>% mutate(dens=ifelse(numdens==1,"3/p","5/p")) %>%
  mutate(method=factor(method, levels=c("PLNnetwork","EMtree","ResampEM","VEMtree"))) %>% 
  mutate(method=fct_recode(method, `PLN-network`="PLNnetwork",nestor="VEMtree")) %>% 
  ggplot(aes(as.factor(dens),AUC,color=as.factor(method), fill=as.factor(method)))+
  geom_boxplot(width=0.4, notch = TRUE, alpha=0.6)+
  facet_grid(~type)+mytheme.dark("")+labs(x="")
ggsave(plot=g,filename = "AUC_PLN_EM_VEM.png", path =  "/Users/raphaellemomal/these/R/images",
       width=7, height=4 )


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


All=do.call(rbind, lapply(c("erdos","cluster"), function(type){
  do.call(rbind, lapply(1:2, function(numdens){
    data=do.call(rbind, lapply(1:100, function(seed){
      listFit=readRDS(paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/simu_seed",seed,"_",type,"_dens",numdens,".rds"))
      ResampEMFreq=readRDS(paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/ResampEM_simu_seed",seed,"_",type,"_dens",numdens,".rds"))$ResampEM
      PLNstars=readRDS(paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/PLNstars_simu_seed",seed,"_",type,"_dens",numdens,".rds"))$PLNstars
      PLNscores=listFit$PLN
      if(sum(!is.na(PLNscores))==0) PLNscores = matrix(0, 15, 15)
      PLNscores=F_Sym2Vec(PLNscores)
      EMprob=F_Sym2Vec(listFit$EM$edges_prob)
      ResampEM=F_Sym2Vec(ResampEMFreq)
      G=F_Sym2Vec(listFit$G)
      return(data.frame(PLNscores=PLNscores, EMprob=EMprob,ResampEM=ResampEM,PLNstars=F_Sym2Vec(PLNstars),G=G))
    }))
    data_VEM=do.call(rbind, lapply(1:100, function(seed){
      listFit=readRDS(paste0("/Users/raphaellemomal/simulations/compar_PLNnet_EM_VEM/simu_seed",seed,"_",type,"_dens",numdens,".rds"))
      if(length(listFit$VEM)>5){
        VEMprob=F_Sym2Vec(listFit$VEM$Pg)
        G=F_Sym2Vec(listFit$G)
        return(data.frame(VEMprob=VEMprob,G=G))
      }}))
    all_PLN=data.frame(seq_PPV_TPR(probs =data$PLNscores, omega = data$G,
                                   seq_seuil = seq(0,max(data$PLNscores),max(data$PLNscores)/200)),
                       type=type, numdens=numdens, method="PLNnetwork")
    all_PLNstars=data.frame(seq_PPV_TPR(probs =data$PLNstars, omega = data$G,
                                   seq_seuil = seq(0,max(data$PLNscores),max(data$PLNscores)/200)),
                       type=type, numdens=numdens, method="PLNstars")
    all_EM=data.frame(seq_PPV_TPR(probs =data$EMprob, omega = data$G,
                                  seq_seuil = seq(0,max(data$EMprob),max(data$EMprob)/200)),
                      type=type, numdens=numdens, method="EMtree")
    all_REM=data.frame(seq_PPV_TPR(probs =data$ResampEM, omega = data$G,
                                  seq_seuil = seq(0,max(data$ResampEM),max(data$ResampEM)/200)),
                      type=type, numdens=numdens, method="ResampEM")
    all_VEM=data.frame(seq_PPV_TPR(probs =data_VEM$VEMprob, omega = data_VEM$G,
                                   seq_seuil = seq(0,max(data_VEM$VEMprob),max(data_VEM$VEMprob)/200)),
                       type=type, numdens=numdens, method="VEMtree")
    all_netVEM=data.frame(seq_PPV_TPR(probs=1*(data_VEM$VEMprob>=0.5), omega = data_VEM$G,
                                   seq_seuil = seq(0,1,1/200)),
                       type=type, numdens=numdens, method="netNestor")
    return(rbind(all_PLN,all_PLNstars, all_EM,all_REM,all_VEM,all_netVEM))
  }))
})) %>% as_tibble()

g2=All %>% as_tibble()  %>%  
  mutate(method=fct_recode(method, `PLN-network`="PLNnetwork",nestor="VEMtree")) %>%
  mutate(dens=ifelse(numdens==1,"3/p","5/p")) %>%  group_by(type, numdens) %>% arrange(seuil, by_group=TRUE) %>% 
  ggplot(aes(TPR, PPV, color=method, shape=method))+geom_hline(yintercept=0.5, linetype="dashed", color="gray")+geom_vline(xintercept=0.5, linetype="dashed", color="gray")+
  geom_point(alpha=0.6)+geom_line(alpha=0.4)+facet_grid(dens~type)+mytheme.dark("")+
  scale_shape_manual(values=c(1,1,1,2,3,4))+
  labs(x="Recall", y="Precision")
ggsave(plot=g2,filename = "precrec_PLN_EM_REM_VEM.png", path =  "/Users/raphaellemomal/these/R/images",
       width=7, height=5 )
g1<-prec_rec %>% as_tibble() %>% filter(method=="PLNnetwork") %>% group_by(type, numdens) %>% arrange(seuil, by_group=TRUE) %>% 
  ggplot(aes(TPR, PPV))+geom_hline(yintercept=0.5, linetype="dashed", color="gray")+geom_vline(xintercept=0.5, linetype="dashed", color="gray")+
  geom_point(alpha=0.2, size=1)+labs(title="PLNnetwork")+
  facet_grid(type~numdens)+mytheme.dark("")
g2<-prec_rec %>% as_tibble() %>% filter(method=="EMtree") %>% group_by(type, numdens) %>% arrange(seuil, by_group=TRUE) %>% 
  ggplot(aes(TPR, PPV))+geom_hline(yintercept=0.5, linetype="dashed", color="gray")+geom_vline(xintercept=0.5, linetype="dashed", color="gray")+
  geom_point(alpha=0.2, size=1)+facet_grid(type~numdens)+mytheme.dark("")+labs(title="EMtree")
g3<-prec_rec %>% as_tibble() %>% filter(method=="VEMtree")%>% group_by(type, numdens) %>% arrange(seuil, by_group=TRUE) %>% 
  ggplot(aes(TPR, PPV))+geom_hline(yintercept=0.5, linetype="dashed", color="gray")+geom_vline(xintercept=0.5, linetype="dashed", color="gray")+
  geom_point(alpha=0.2, size=1)+facet_grid(type~numdens)+mytheme.dark("")+labs(title="VEMtree")
grid.arrange(g1, g2, g3, ncol=3)
# pour PLNnetwork
model_StARS <- getBestModel(network_models, "StARS")
my_graph <- plot(model_StARS, plot = FALSE)
ggimage(as.matrix(as_adj(my_graph)))
test



############
# compar proba
pooled %>% ggplot(aes(nestor, EMtree))+geom_point(size=1, alpha=0.2)+facet_grid(type~numdens)+theme_light()
pooled %>% mutate(numdens=ifelse(numdens==1,"3/p","5/p")) %>%  gather(key, value,-G,-type,-numdens,-seed,-`PLN-network`) %>% 
  ggplot(aes(value, color=key, fill=key))+geom_histogram(aes(y=..density..), alpha=0.3, bins=50,position="identity")+
  coord_cartesian(xlim=c(0,1))+
  facet_grid(numdens~type)+mytheme.dark("")

Data=pooled %>% filter(type=="cluster", numdens==2, nestor<1.000001)
pmain=  ggplot(Data,aes( EMtree, nestor))+geom_point(size=1,alpha=0.2)+theme_light()


xdens <- axis_canvas(pmain, axis = "x")+
  geom_histogram(data = Data, aes(x = EMtree, y = ..density..),
                 alpha = 0.5, size = 0.2,bins=40, color="#FA7E1E", fill="#FA7E1E")  
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_histogram(data = Data, aes(x = nestor, y = ..density..),
                 alpha = 0.5, size = 0.2,bins=40, color="#D62976", fill="#D62976")+ coord_flip()  
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

ggsave(plot=p2,filename = "probEM_nestor.png", path =  "/Users/raphaellemomal/these/R/images",
       width=5, height=4 )

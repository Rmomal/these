#####
# tableaux de comparaison TPFPTNFN avec deux niveaux de difficulté
####

library(ROCR)
library(PLNmodels)
library(SpiecEasi)
library(MInt)
library(parallel)
library(tidyverse)
library(ggbeeswarm)
library(grid)
library(gridExtra)
library(xtable)
source('~/simulations/codes/codegCoda.R')
source('~/simulations/codes/fonctions.R')
source('~/simulations/codes/FunctionsInference.R')
source('~/simulations/codes/Kit_pour_VEM_EM.R')
##########

eval_store_mint<-function(Y,covar,path){
  #  browser()
  data<-data_for_MInt(Y,covar,path)
  x <- data[["x"]]
  y <- data[["y"]]
  T1<-Sys.time()
  m <- mint(y,x,fmla = ~feature1 + feature2)
  m <- estimate(m)
  T2<-Sys.time()
  temps<-difftime(T2,T1)
  pred<-m$param$P
  #obs<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,".rds"))$omega
  #roc_curve(pred,obs),
  # res=list(pred=pred,temps= temps)
  return(pred)
}
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"


#############################################################################################
#########-----------------                  DATA                -----------------############                  
#############################################################################################
# covar<-data.frame(X1=round(runif(n)*5),X2=rnorm(n,0,2),
#                   X3=(runif(n)))
# saveRDS(covar,paste0(path,"mint_covarn50.rds"))
types<-c("erdos","cluster","scale-free")
diff=c("easy","hard")
for( difficulty in diff){
  # n<-switch(difficulty,"easy"=100,"hard"=50)
  # nbspecies<-switch(difficulty,"easy"=20,"hard"=30)
  # covar<-data.frame(X1=round(runif(n)*5),X2=rnorm(n,0,2),  X3=(runif(n)))
  # saveRDS(covar,paste0(path,"TPFN/Data/covar_",difficulty,".rds"))
  # 
  for(type in types){
    cat(type,difficulty,"\n")
    obj<-mclapply(1:100,function(nbgraph){
      cat("\n",nbgraph,":")
      # omegaOriginal <- readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLN.2.0/",type,
      #                                 "/d/Sets_param/Graph",nbgraph,"_",nbspiecies,".rds"))$omega
      # edgesOrigin<-ifelse(abs(F_Sym2Vec(omegaOriginal))<1e-16,0,1)
      # Y<-data_from_stored_graphs(type, variable, nbgraph, nbspecies, covar=covar,path,fixe=FALSE)
      dat <- data_from_scratch(type, p=nbspecies, r=10, covar=covar, signed=TRUE)
      saveRDS(dat,paste0(path,"TPFN/Data/Signed.Data_",type,"_",difficulty,nbgraph,".rds"))
      
    }, mc.cores=1)
  }
}

#############################################################################################
#########--------------                  DISAG RATE                --------------############                  
#############################################################################################
res<-do.call(rbind, lapply(diff, function(difficulty){
  obj<- lapply(c(2,10,20,50,100,150), function(B){
    EMconverged<-readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLN.2.0/TPFN/Results/GhatEMtree_erdos_TFPN_",difficulty,"150.rds"))
    tocompare<-readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLN.2.0/TPFN/Results/GhatEMtree_erdos_TFPN_",difficulty,B,".rds"))
    obj<-lapply(1:100, function(x){
      tpfn<-c(table(tocompare[[x]],EMconverged[[x]]))
      res<-data.frame("TN"=tpfn[1],"FP"=tpfn[2],"FN"=tpfn[3],"TP"=tpfn[4],S=B,
                      type=type,difficulty=difficulty,graph=x)
      return(res)
    })
    res<-as_tibble(do.call(rbind,obj))
  })
  res<-as_tibble(do.call(rbind,obj))%>%
    mutate(TN=as.numeric(as.character(TN)),
           FN=as.numeric(as.character(FN)),
           TP=as.numeric(as.character(TP)),
           FP=as.numeric(as.character(FP)),
           sum=TN+TP+FP+FN,FDR=FP/(TP+FP))
  
  res<-res%>% mutate(FN=ifelse(is.na(sum),FP,FN),FP=ifelse(is.na(sum),0,FP),
                     TP=ifelse(is.na(sum),0,TP))
  
  
}))
saveRDS(res,paste0(path,"/TPFN/results/",method,"_",type,"_TFPN_disagreementRate.rds"))

### plot disag rate
p1<-res %>% filter(difficulty=="easy") %>%
  ggplot(aes( y=FDR, x=as.factor(S), color=as.factor(S)))+geom_quasirandom()+theme_minimal()+
  scale_color_manual(values=c(brewer.pal(n = 9, name = "YlOrRd")[-c(1:3)],"black"))+labs(x="")+
  guides(color=FALSE)+
  geom_hline(yintercept = 0.05,linetype="dashed")+
  stat_summary(aes(y=FDR),fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               fill="white",color="black",pch=22)+ labs(title="Easy",y="Disagreement rate",x="S")+
  theme(plot.title = element_text(hjust = 0.5))
p2<-res %>% filter(difficulty=="hard") %>%
  ggplot(aes( y=FDR, x=as.factor(S), color=as.factor(S)))+geom_quasirandom()+theme_minimal()+
  scale_color_manual(values=c(brewer.pal(n = 9, name = "YlOrRd")[-c(1:3)],"black"))+labs(x="")+
  guides(color=FALSE)+
  geom_hline(yintercept = 0.05,linetype="dashed")+
  
  stat_summary(aes(y=FDR),fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               fill="white",color="black",pch=22)+ labs(title="Hard",y="Disagreement rate",x="S")+
  theme(plot.title = element_text(hjust = 0.5))
p<-grid.arrange(p1,p2,nrow=1,ncol=2)
ggsave(paste0(path,"/TPFN/DisagRateErdos.png"),plot=p,height=5,width=7)

################################################################################################################
################################                  COMPUTE TPFN                  ################################                  
################################################################################################################

TPFN_compute<-function(methods,diffs,types, cores=3,B=100,S=10){
  for(method in methods){
    for( difficulty in diffs){#c(1e-3,1e-2,1e-1,1,10
      #   for(landrat in c(1e-4,1e-3,1e-2,1e-1)){
      #  for(S in c(2,5,10,20)){
      
      cond.tol<-switch(difficulty,"easy"=1e-12,"hard"=1e-6)
      for(type in types){
        T1<-Sys.time()
        cat(method,", ",difficulty,type,": ")
        
        
        obj<-mclapply(1:B,function(nbgraph){# attention signed data !! Signed.Data_
          ########################
          # récupérer les données et covariables associées générées précédemment, possibilités graphs signés
          dat<-readRDS(paste0(path,"TPFN/Data/Signed.Data_",type,"_",difficulty,nbgraph,".rds"))
          Y<-dat[[1]]
          edgesOrigin<-ifelse(abs(F_Sym2Vec(dat[[2]]))<1e-16,0,1)
          covar <- readRDS(paste0(path,"TPFN/Data/covar_",difficulty,".rds"))
          m<- model.matrix(~X1+X2+X3,covar)# l'intercept est enlevé dans le PLN de resample
          p=ncol(Y)
          ########################
          # inférence du réseau par les différentes méthodes
          if(method=="MRFcov"){
            T1<-Sys.time()
            mrfres<- MRFcov(data = cbind(Y,m[,-1]), n_nodes = p, family = 'poisson', 
                            symmetrise = "mean",n_cores=3, n_covariates = 3)$graph
            inf=F_Sym2Vec(1*(mrfres!=0))
            T2<-Sys.time()
            time<-difftime(T2,T1)
          }
          if(method=="MInt"){
            T1<-Sys.time()
            inf <- F_Sym2Vec(eval_store_mint(Y,covar,path)>1e-16)*1
            T2<-Sys.time()
            time<-difftime(T2,T1)
          }
          if(method=="EMtree"){
            T1<-Sys.time()
            resample<-ResampleEMtree(counts=Y, covar_matrix=m, S=S, maxIter=200,  cond.tol=cond.tol,cores=3)
            if(length(resample)!=3) browser()
            pmat<-resample$Pmat
            # inf<- 1*(ifelse(colSums(ifelse(pmat<2/p,0,1))/B >0.8,1,0))
            inf<-F_Sym2Vec(1*(freq_selec(pmat, Pt=0.5)>0.5))
            T2<-Sys.time()
            time<-difftime(T2,T1)
          }
          if(method=="gCoda"){
            T1<-Sys.time()
            inf <- F_Sym2Vec(gcoda(Y, counts=T, covar=covar, nlambda=30, lambda.min.ratio =1)$opt.icov>1e-16)*1
            T2<-Sys.time()
            time<-difftime(T2,T1)
          }
          if(method=="SpiecEasi"){
            U<-t(clr.matrix(Y,mar=1))
            model<-lm(U~m)
            T1<-Sys.time()
            inf<-spiec.easi(model$residuals, method='mb', lambda.min.ratio=1e-3, nlambda=100,
                            icov.select.params=list(rep.num=2 ))
            inf<-F_Sym2Vec(as.matrix(inf$refit[[1]]))
            T2<-Sys.time()
            time<-difftime(T2,T1)
          }
          if(method=="ecoCopula"){
            T1<-Sys.time()
            my_mod=manyglm(Y~m, family="negativ.binomial")
            inf = F_Sym2Vec(cgr(my_mod)$best_graph$graph)
            T2<-Sys.time()
            time<-difftime(T2,T1)
          }
          ########################
          # comparaison de inf à edgesOrigin
          tpfn<-c(table(inf,edgesOrigin))
          res<-data.frame("TN"=tpfn[1],"FP"=tpfn[2],"FN"=tpfn[3],"TP"=tpfn[4] ,method=method,
                          times=time,unit=attr(time, "units"), type=type,difficulty=difficulty,
                          graph=nbgraph)
          return(res)
        },mc.cores=cores)
        ########################
        # obj est une liste de taille B qui contient les TPFN pour une méthode*type de graph*difficulté
        res<-as_tibble(do.call(rbind,obj))%>%
          mutate(TN=as.numeric(as.character(TN)),
                 FN=as.numeric(as.character(FN)),
                 TP=as.numeric(as.character(TP)),
                 FP=as.numeric(as.character(FP)),
                 times=ifelse(unit=="mins",60*as.numeric(as.character(times)),as.numeric(as.character(times))),
                 sum=TN+TP+FP+FN,FDR=FP/(TP+FP))
        
        res<-res%>% mutate(FN=ifelse(is.na(sum),FP,FN),FP=ifelse(is.na(sum),0,FP),
                           TP=ifelse(is.na(sum),0,TP))
        ########################
        # res est un tibble qi contient pour B graphs d'un type*difficulté, les TPFN corrigés, times et FDR pour une méthode
        saveRDS(res, paste0(path,"/TPFN/Results/",method,"_",type,"_TFPN_",difficulty,".rds"))
        #   saveRDS(res,paste0(path,"/TPFN/results/",method,"_",type,"_TFPN_",difficulty,".rds"))
        #   saveRDS(lapply(obj, function(x){x[[2]]}), paste0(path,"/TPFN/results/Ghat",method,"_",type,"_TFPN_",difficulty,".rds"))
        T2<-Sys.time()
        cat(difftime(T2,T1),attr(difftime(T2,T1), "units"),"\n")
        # times<-c(times,difftime(T2,T1))
      }
    }
    #  }
    #    }
  }
}

#######@######################################@@@@@@@@@@####################################
#########--------------     S EFFECT ON EMTREE/SPIECEASI PERF     --------------############                  
#######@#####################################@@@@@@@@@@#####################################

res=readRDS(paste0(path,"/TPFN/results/Pmat",method,"_",type,"_TFPN_",difficulty,"2.rds"))
# S effect
type<-"erdos"
difficulty<-"easy"
method<-"SpiecEasi"
files=c()


for(landrat in c(1e-4,1e-3,1e-2,1e-1)){
  files=c(files,paste0(path,"/TPFN/Results/",method,"_",type,"_TFPN_",difficulty,"_",landrat,".rds"))
}


reseasy<-do.call(rbind, lapply(files,function(x){
  split=strsplit(x ,split="_")[[1]]
  
  landrat=substr(split[5],1,nchar(split[5])-4)
  cbind(readRDS(x), iterStars=S,
        landrat=landrat)
  
}))

reshard<-do.call(rbind, lapply(B,function(x){
  readRDS(paste0(path,"TPFN/results/",method,"_",type,"_TFPN_hard",x,".rds"))
}))
plotFDRratio<-function(res,diff){
  
  res<- res %>% mutate( FDR=FP/(TP+FP) ,densPred=(TP+FP)/(TP+FN))
  ncolors=length(unique(res$landrat))
  p <-   ggplot(res,aes( x=as.factor(landrat), color=as.factor(landrat)))+
    scale_color_manual(values=c(brewer.pal(n = ncolors+2, name = "YlOrRd")[-c(1,2)],"black"))+labs(x="")+
    guides(color=FALSE)+
    theme_minimal()+
    theme(axis.text.x = element_text( hjust = 1))
  
  
  p1<-p+geom_quasirandom(aes(y=FDR))+stat_summary(aes(y=FDR),fun.ymin = function(z) { quantile(z,0.25) },
                                                  fun.ymax = function(z) { quantile(z,0.75) },
                                                  fun.y = median,
                                                  fill="white",color="black",pch=22)+ labs(y="FDR %",x="lambda.min.ratio")
  #  coord_cartesian(ylim=c(0,0.8))
  
  p2<-p+geom_quasirandom(aes(y=log(densPred+0.01)))+stat_summary(aes(y=log(densPred+0.01)),fun.ymin = function(z) { quantile(z,0.25) },
                                                                 fun.ymax = function(z) { quantile(z,0.75) },
                                                                 fun.y = median,
                                                                 fill="white",color="black",pch=22)+
    labs(y="log(ratio + c)", x="lambda.min.ratio")+
    scale_y_continuous(breaks = c(log(0.01),log(0.1),log(0.5),log(1), log(2), log(4)),
                       labels=c("log(c)","log(1e-1)","log(5e-1)","log(1)", "log(2)","log(4)"))+
    geom_hline(yintercept = 0 ,linetype="dashed")
  #  coord_cartesian(ylim=c(log(3e-1),log(3)))
  
  grid.arrange(p1,p2,nrow=2,ncol=1, top=diff)
}

g1<-plotFDRratio(reseasy,"Easy (nlambda=100)")
ggsave("landrat_effect_gCoda.png",plot=g1, path=paste0(path,"/images"), width=6.5, height=7)
g2<-plotFDRratio(reshard,"Hard")
p<-grid.arrange(g1,g2,nrow=1,ncol=2)

ggsave("S_effect.png",plot=p,  width=7, height=5.5)

#### times per S

times<- res%>% select(method,times) %>%
  group_by_at(vars(-times)) %>%
  mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
  spread(key=method, value=times) %>%    # spread
  select(EMtree)
res%>% select(method,FDR)
group_by_at(vars(-FDR)) %>%
  mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
  spread(key=method, value=FDR) %>%    # spread
  select(EM1,EMTree) %>% mutate(diff=EMTree-EM1,time=times$EMTree) %>%
  ggplot(aes(time,diff))+geom_point()+theme_minimal()+geom_hline(yintercept = 0,col="red")

#######@#########################################################################################
#########--------------     TIMES AND EMPTY FOR GCODA AND SPIECEASI    --------------############                  
#######@#########################################################################################

difficulty = c("easy","hard")
res<-do.call(rbind, lapply(difficulty, function(difficulty){
  type<- c("erdos","cluster","scale-free")
  obj<-do.call(rbind, lapply(type,function(type){
    method <- c("EMtree")
    cols<-c("FP","TP","method" ,"times","unit","type", "difficulty", "medIter",
            "thirdQuartile","maxIter")
    res<-do.call(rbind, lapply(method, function(meth){
      if(meth=="EMtree"){
        obj<-readRDS(paste0(path,"TPFN/results/",meth,"_",type,"_TFPN_", difficulty,"150.rds")) %>% select(cols)
      }else{
        obj <-readRDS(paste0(path,"TPFN/results/",meth,"_",type,"_TFPN_", difficulty,".rds")) %>% select(cols)
      }
      
    }))
    res <- res%>%  mutate(detect=as.numeric(as.character(TP))+as.numeric(as.character(FP)),
                          times=as.numeric(as.character(times)))
  }))
}))

res %>% group_by(difficulty,method,type) %>%filter(method%in%c("gCoda","SpiecEasi")) %>%
  summarise(nul=length(which(detect==0))) %>%
  spread(method,nul, drop = TRUE)  %>% group_by(difficulty) %>%
  summarise(mgcoda=mean(gCoda),mspiec=mean(SpiecEasi))


res %>% group_by(difficulty,method) %>%
  summarise(running_times=paste0(round(median(times),2)," (",round(sd(times),2),")")) %>%
  spread(method,running_times)


xtable(res %>% group_by(difficulty,method) %>%
         summarise(running_times=paste0(round(median(times),2)," (",round(sd(times),2),")")) %>%
         spread(method,running_times) )

resErdos<-rbind(reseasy, reshard)
xtable(resErdos %>% group_by(B, difficulty) %>%
         summarise(med=paste0(round(median(times),2),"(", sd=round(sd(times),2),")")) %>%
         spread(B,med))

#######@#####################################################################################
#########-------------------     Ft EFFECT ON EMTREE PERF    --------------------############                  
#######@#####################################################################################

type<-"erdos";difficulty<-"easy";method<-"EMtree"; fseq=seq(0.8,1,0.02)
TPFNseuil<-function(pmat,edgesOrigin,p,f,graph){
  inf<- 1*(ifelse(colSums(ifelse(pmat<2/p,0,1))/100>f,1,0))
  tpfn<-c(table(inf,edgesOrigin))
  if(length(tpfn)!=4){
    res<-data.frame("TN"=tpfn[1],"FP"=0,"FN"=tpfn[2],"TP"=0 ,difficulty=difficulty,graph=graph, seuil=f)
  }else{
    res<-data.frame("TN"=tpfn[1],"FP"=tpfn[2],"FN"=tpfn[3],"TP"=tpfn[4] ,difficulty=difficulty,graph=graph, seuil=f)
  }
  return(res)
}
infVariaSeuil<-function(difficulty){
  p<-switch(difficulty,"easy"=20,"hard"=30)
  obj<-readRDS(paste0(path,"/TPFN/results/PmatEMtree_erdos_TFPN_",difficulty,".rds"))
  do.call(rbind, lapply(1:100,function(x){
    pmat<-obj[[x]][[1]]
    edgesOrigin<-obj[[x]][[2]]
    res<-do.call(rbind,lapply(fseq,function(f){TPFNseuil(pmat,edgesOrigin,p,f,x)}))
    return(res)
  }))
}
reseasy=infVariaSeuil("easy")
reashard=infVariaSeuil("hard")
plotFDRratio<-function(res,diff){
  res<- res %>% mutate( FDR=FP/(TP+FP) ,densPred=(TP+FP)/(TP+FN)) %>% filter(seuil!=1)
  p <- ggplot(res,aes( x=as.factor(seuil)))+
    guides(color=FALSE)+ theme_minimal()+ labs(x="Threshold")+theme(axis.text.x = element_text( hjust = 1))
  
  p1<-p+geom_quasirandom(aes(y=FDR), color="deepskyblue3")+stat_summary(aes(y=FDR),fun.ymin = function(z) { quantile(z,0.25) },
                                                                        fun.ymax = function(z) { quantile(z,0.75) },
                                                                        fun.y = median,
                                                                        fill="white",color="black",pch=22)+ labs(y="FDR %")#+ coord_cartesian(ylim=c(0,0.8))
  
  p2<-p+geom_quasirandom(aes(y=log(densPred+0.01)), color="deepskyblue3")+stat_summary(aes(y=log(densPred+0.01)),fun.ymin = function(z) { quantile(z,0.25) },
                                                                                       fun.ymax = function(z) { quantile(z,0.75) },
                                                                                       fun.y = median,
                                                                                       fill="white",color="black",pch=22)+
    labs(y="log(ratio + c)")+
    scale_y_continuous(breaks = c(log(0.01),log(0.1),log(0.5),log(1), log(2), log(4)),
                       labels=c("log(c)","log(1e-1)","log(5e-1)","log(1)", "log(2)","log(4)"))+
    geom_hline(yintercept = 0 ,linetype="dashed")#+ coord_cartesian(ylim=c(log(3e-1),log(3)))
  
  grid.arrange(p1,p2,nrow=2,ncol=1, top=diff)
}

g1<-plotFDRratio(reseasy,"Easy")
g2<-plotFDRratio(reashard,"Hard")
p<-grid.arrange(g1,g2,nrow=1,ncol=2)

ggsave("Threshold_effect.png",plot=p,  width=14, height=8)

#######@#####################################################################################
#########--------------------------        TPFN PLOTS        --------------------############                  
#######@#####################################################################################

build_TFPN_plots<-function(type,difficulty,method,factors=methods,colors,FDR=FALSE, TPFN=FALSE,FOR=FALSE, senspe=FALSE){
  cols<-c("TN","FP" ,"FN"  ,"TP","method" ,"times","unit","type", "difficulty" ,"graph",  "sum","FDR")
  res<-do.call(rbind, lapply(method, function(meth){
    obj <-readRDS(paste0(path,"TPFN/Results/",meth,"_",type,"_TFPN_", difficulty,".rds")) %>% dplyr::select(cols)
    obj$method=meth
    
    return(obj)
  }))
  res<-res%>% mutate(method=fct_relevel(method,factors),# %>% mutate(method = fct_recode(method, "EMtree-10" = "EM"))
                     FDR=FP/(TP+FP),  densPred=(TP+FP)/(TP+FN), FOR=FN/(FN+TN),
                     sens=TP/(TP+FN), spe=TN/(TN+FP))
  p <-ggplot(res,aes( x=method, color=method))+
    scale_color_manual(values=colors)+labs(x="")+
    guides(color=FALSE)+ theme_minimal()+ theme(axis.text.x = element_text(angle = 25, hjust = 1))
  if (FDR) {
    p=p+geom_quasirandom(aes(y=FDR))+
      stat_summary(aes(y=FDR),fun.ymin = function(z) { quantile(z,0.25) },
                   fun.ymax = function(z) { quantile(z,0.75) },
                   fun.y = median,
                   fill="white",color="black",pch=22)+
      labs(y="FDR %",x="")
    return(p)
  }else{
    p+geom_quasirandom(aes(y=log(densPred+0.01)))+
      stat_summary(aes(y=log(densPred+0.01)),fun.ymin = function(z) { quantile(z,0.25) },
                   fun.ymax = function(z) { quantile(z,0.75) },
                   fun.y = median,
                   fill="white",color="black",pch=22)+
      scale_y_continuous(breaks = c(log(0.01),log(0.1),log(0.5),log(1), log(2), log(4)),
                         labels=c("c",".1",".5","1", "2","4"))+
      labs(y="Density ratio")+geom_hline(yintercept = 0 ,linetype="dashed")
  }
  
  # if(TPFN){
  #   p1=p+geom_quasirandom(aes(y=FP))+
  #     stat_summary(aes(y=FP),fun.ymin = function(z) { quantile(z,0.25) },
  #                  fun.ymax = function(z) { quantile(z,0.75) },
  #                  fun.y = median,
  #                  fill="white",color="black",pch=22)+
  #     labs(y="FP ",x="")
  #   p2=p+geom_quasirandom(aes(y=FN))+
  #     stat_summary(aes(y=FN),fun.ymin = function(z) { quantile(z,0.25) },
  #                  fun.ymax = function(z) { quantile(z,0.75) },
  #                  fun.y = median,
  #                  fill="white",color="black",pch=22)+
  #     labs(y="FN ",x="")
  #   grid.arrange(p1,p2,nrow=1,ncol=2)
  # }
  # if(FOR){
  #   p1=p+geom_quasirandom(aes(y=FDR))+
  #     stat_summary(aes(y=FDR),fun.ymin = function(z) { quantile(z,0.25) },
  #                  fun.ymax = function(z) { quantile(z,0.75) },
  #                  fun.y = median,
  #                  fill="white",color="black",pch=22)+
  #     labs(y="FDR %",x="")
  #   p2= p+geom_quasirandom(aes(y=FOR))+
  #     stat_summary(aes(y=FOR),fun.ymin = function(z) { quantile(z,0.25) },
  #                  fun.ymax = function(z) { quantile(z,0.75) },
  #                  fun.y = median,
  #                  fill="white",color="black",pch=22)+
  #     labs(y="FOR %",x="")
  #   grid.arrange(p1,p2,nrow=1,ncol=2)
  # }
  # if(senspe){
  #   p1=p+geom_quasirandom(aes(y=sens))+
  #     stat_summary(aes(y=sens),fun.ymin = function(z) { quantile(z,0.25) },
  #                  fun.ymax = function(z) { quantile(z,0.75) },
  #                  fun.y = median,
  #                  fill="white",color="black",pch=22)+
  #     labs(y="Sensitivity",x="")
  #   p2= p+geom_quasirandom(aes(y=spe))+
  #     stat_summary(aes(y=spe),fun.ymin = function(z) { quantile(z,0.25) },
  #                  fun.ymax = function(z) { quantile(z,0.75) },
  #                  fun.y = median,
  #                  fill="white",color="black",pch=22)+
  #     labs(y="Specificity",x="")
  #   grid.arrange(p1,p2,nrow=1,ncol=2)
  # }
}

gridLine<-function(type,method,factors,colors,top=TRUE){
  label<-switch(type,"erdos"="Erdös \n \n","cluster"="Cluster \n \n",
                "scale-free"="Scale-free \n \n")
  if(top){
    top1="Easy (n=100, p=20)\n"; top2="Hard (n=50, p=30)\n"
  }else{
    top1="";top2=""
  }
  p1<-build_TFPN_plots(type = type,difficulty ="easy",method=method,factors=factors,colors=colors,
                       FDR = TRUE)

  p2<-build_TFPN_plots(type = type,difficulty ="easy",method=method,factors=factors,colors=colors,
                       FDR = FALSE)
  g1<-grid.arrange(p1,p2, ncol=2, nrow=1,top=top1,right="\n")
  p1<-build_TFPN_plots(type = type,difficulty ="hard",method=method,factors=factors,colors=colors,
                       FDR = TRUE)
  p2<-build_TFPN_plots(type = type,difficulty ="hard",method=method,factors=factors,colors=colors,
                       FDR = FALSE)
  g2<-grid.arrange(p1,p2, ncol=2, nrow=1,top=top2, right=label)
  
  grid.arrange(g1,g2,nrow=1,ncol=2)
}

################################################################################################
#################################        EASY RUN       ########################################
################################################################################################

TPFN_compute(methods = c("SpiecEasi","gCoda"),
             diffs=c("easy","hard"),
             types=c("erdos","cluster","scale-free"),
             B=100, cores=3)

#methods:
# c("EMtree","MInt","gCoda","SpiecEasi")
# c("MRFcovSigned","MRFcovmean1")
# c("ecoCopula","MRFcov1","MRFcov2","MInt","EMtree")
# c("MRFcovmin1","MRFcovmin2","MRFcovmean1","MRFcovmean2","MRFcovmax1","MRFcovmax2","EMtree")

#colors:
#c("indianred2","orange2","#6abd35","steelblue4")
#c("#24babf","#aa2061","#e67caf","#6abd35","steelblue4")
#c("deepskyblue1","deepskyblue2","olivedrab2","olivedrab3","orangered","orangered3","steelblue4")
type=c("erdos","cluster","scale-free")
method=c("SpiecEasi","gCoda")
factors=method
colors=c("indianred2","orange2")
p1=build_TFPN_plots(type = type,difficulty =c("easy","hard"),method=method,factors=factors,colors=colors,FDR = TRUE)
p2=build_TFPN_plots(type =type ,difficulty = c("easy","hard"),FDR = FALSE)
grid.arrange(p1,p2, nrow=1, ncol=2)

E<-gridLine("erdos",method,factors,colors,TRUE)
C<-gridLine("cluster",method,factors,colors,FALSE)
SF<-gridLine("scale-free",method,factors,colors,FALSE)
grid.arrange(E,C,SF,nrow=3,ncol=1)
p<-arrangeGrob(E,C,SF,nrow=3,ncol=1)
ggsave("panel_TPFN-ecocopula.png", plot=p, width=10, height=9)



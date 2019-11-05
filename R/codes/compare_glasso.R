# compare EMtree and glasso
library(mvtnorm)
library(ROCR)
library(parallel)
library(glasso)
library(tidyverse)
library(EMtree)
library(Matrix)
theme_set(theme_light())
data=FALSE
# à charger aussi : diagnost_auc et vec_obs_pred de "fonctions.R"
##################
##################
EMtreeGaussien<-function(Y, maxIter=30, cond.tol=1e-10, verbatim=FALSE, plot=FALSE){
  CorY=cov2cor(cov(Y))
  p = ncol(CorY)
  alpha.psi = Psi_alpha(CorY, n, cond.tol=cond.tol)
  psi = alpha.psi$psi
  n=nrow(Y)
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)
  
  FitEM = FitBetaStatic(beta.init=beta.unif, psi=psi, maxIter = maxIter,
                        verbatim=verbatim, plot=plot)
  
  return(FitEM)
}
ResampleEMtreeGaussien<-function(data,S, cores=1 ,v=0.8){
  n=nrow(data);V = round(v * n)
  obj<-mclapply(1:S,function(b){
    cat("\nS=",b," ")
    set.seed(b)
    
    sample = sample(1:n, V, replace = F)
    data.sample = data[sample,]
    EMtreeGaussien(data.sample)
    
  }, mc.cores=cores)
  Pmat<-do.call(rbind,lapply(obj,function(x){F_Sym2Vec(x$edges_prob)}))
  return(Pmat)
}
diagnost_auc<-function(obs, pred){
  obs_pred<-vec_obs_pred(obs,pred)
  prediction<-prediction(obs_pred[[1]],obs_pred[[2]])
  # Run the AUC calculations
  ROC_auc <- performance(prediction,"auc")
  res<-round(ROC_auc@y.values[[1]],digits=3)
  return(res)
}

vec_obs_pred<-function(obs, pred){
  nvar<-ncol(obs)
  obs[which(abs(obs)<1e-16)]<-0
  indices_nuls<-which(obs==0)
  label<-matrix(1,nrow=nrow(obs),ncol=ncol(obs))
  label[indices_nuls]<-0
  
  vec_pred<-as.vector(pred[upper.tri(pred)])
  vec_obs<-as.vector(label[upper.tri(label)])
  
  return(list(vec_pred,vec_obs))
}
##################
# data
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
types<-c("scale-free")
#generator param est redéfini dans gener_data dans pacage EMtree
#attention, charger le gener_param de Kit_pour_VEM
origin_graph=function(type, p=20, r=10, covar=NULL,prob=log(p)/p,dens=log(p)/p, v=0.01,signed=FALSE, n=NULL){
  graph<- generator_graph(graph=type,p=p,prob=prob,dens=dens,r=r)
  if(signed){
    param<-generator_param(as.matrix(graph), signed=TRUE,v=v)
  }else{
    param<-generator_param(as.matrix(graph), signed=FALSE,v=v)
  }
  return(param$omega)
}
vecnoise=c(TRUE)
if(data){
  for( noise in vecnoise){
    rep.noise=ifelse(noise,"","lessNoise")
    v=ifelse(noise,1,0.01)
    for(type in types){
      mclapply(1:100,function(nbgraph){
        cat("\n",nbgraph,":")
        omega= origin_graph(type,p=25, signed=TRUE,v=v)
        
        saveRDS(omega,paste0(path,"compare_glasso/graphs_",rep.noise,"/",type,"_",nbgraph,".rds"))
        
      }, mc.cores=3)
    }
  }
  
}
##################
#simulation



compare_glasso<-function(methods,seqn=seq(10,40,10),types, cores=3,B=100,S=10, maxpen=FALSE,noise=FALSE){
  rep.noise=ifelse(noise,"","lessNoise")
  for(method in methods){
    for( n in seqn){
      path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
      
      for(type in types){
        T1<-Sys.time()
        cat(method,type,", n=",n,":")
        
        
        obj<-mclapply(1:B,function(nbgraph){
          cat("\ngraph ",nbgraph)
          ########################
          # récupérer les graphs signés et générer les données gaussiennes
          omega<-readRDS(paste0(path,"compare_glasso/graphs_",rep.noise,"/",type,"_",nbgraph,".rds"))
          sigma=cov2cor(solve(omega))
          Y<- rmvnorm(n, rep(0,nrow(sigma)), sigma)
          edgesOrigin<-ifelse(abs(F_Sym2Vec(omega))<1e-16,0,1)
          p=ncol(Y)
          ########################
          # inférence du réseau par les différentes méthodes
          
          if(method=="EMtree_noresamp"){
            T1<-Sys.time()
            # browser()
            #  resample<-ResampleEMtreeGaussien(data = Y,S = 100,cores=1)
            preinf<-EMtreeGaussien(Y = Y)
            #  inf<-freq_selec(resample, Pt=2/p)
            inf<-preinf$edges_prob
            T2<-Sys.time()
            time<-difftime(T2,T1)
          }
          if(method=="glasso_maxpen"){
            T1<-Sys.time()
            pathglass <- glassopath(s = cov(Y))
            res=lapply(seq_along(pathglass$rholist),function(x){
              pathglass$rholist[x]*(1*(pathglass$wi[,,x]!=0))
            })
            K.score <- F_Sym2Vec(Reduce("+",res))
            
            # pénalité qui tue l'arête
            cumrho=cumsum(pathglass$rholist)
            for( i in 1:length(K.score)){
              if(K.score[i]%in%cumrho){
                index=which(cumrho==K.score[i])
                K.score[i]<-pathglass$rholist[index]
              }
              inf=F_Vec2Sym(K.score)
            }
            T2<-Sys.time()
            time<-difftime(T2,T1)
          }
          if(method=="glasso_cumsum"){
            T1<-Sys.time()
            pathglass <- glassopath(s = cov(Y))
            res=lapply(seq_along(pathglass$rholist),function(x){
              pathglass$rholist[x]*(1*(pathglass$wi[,,x]!=0))
            })
            
            K.score <-Reduce("+",res)
            # somme cumulée des pénalités
            inf<- K.score / max(K.score)
            T2<-Sys.time()
            time<-difftime(T2,T1)
          }
          #compute auc
          
          auc= diagnost_auc(F_Vec2Sym(edgesOrigin),inf)
          res=data.frame(method=method,type=type,auc=auc,graph=nbgraph,n=n)
          return(res)
        },mc.cores=cores)
        
        # save results
        
        res= do.call(rbind,obj)
        saveRDS(res,file=paste0(path,"compare_glasso/results_",rep.noise,"/",method,"_",type,"_",n,".rds"))
        T2<-Sys.time()
        cat(difftime(T2,T1),attr(difftime(T2,T1), "units"),"\n")
      }
    }
  }
}
# pour glasso ne pas lancer scale-free en même temps que erdos et cluster, va savoir pourquoi
#si scale-free prend plus de 3s, relancer
# sur scale-free noise FALSE, glasso peine à n=10
compare_glasso(c("glasso_maxpen"),types=c("scale-free"), seqn=10, cores=3, noise=FALSE)
files=c()
for(method in c("EMtree_noresamp","glasso_cumsum","glasso_maxpen")){
  for( type in c("cluster","erdos","scale-free")){
    for (n in seq(10,60,10)){
      files=c(files,(paste0(path,"compare_glasso/results/",method,"_",type,"_",n,".rds")))
      files=c(files,(paste0(path,"compare_glasso/results_lessNoise/",method,"_",type,"_",n,".rds")))
      
    }
  }
}

list=lapply(1:length(files),function(x){
  print(x)
  item=strsplit(files[x],split="_")[[1]][4]
  n=substr(item,1,nchar(item)-4)
  
  noise=ifelse(strsplit(files[x],split="/")[[1]][8]=="results","noise +","noise -")
  if(ncol(readRDS(files[x]))!=5){
    cbind(readRDS(files[x]),n=n, noise=noise)
  }else{
    cbind( readRDS(files[x]), noise=noise)
  }
})

list=lapply(list, function(x){
  x[,"n"] = as.factor(x[,"n"])
  return(x)
})
list=do.call(rbind,list)
list=list %>% as_tibble() %>% mutate(method=case_when(method=="EMtree_noresamp"~"EMtree",
                                                      method== "glasso_cumsum"~"glassoSomme",
                                                      method== "glasso_maxpen"~"glassoPen"))
plot=list  %>% ggplot(aes(x=n,y=auc, color=method))+geom_boxplot(width=0.5)+
  facet_grid(type~noise)+theme_light()+
  scale_color_brewer("Method",palette="Dark2")+labs(title="25 nodes",
                                                    x="number of samples",y="AUC")+
  theme( strip.text.x = element_text( size = 12 ) ,
         strip.background.x=element_rect(fill = "gray50") )

ggsave(plot=plot, filename = "compare_glasso_noise_SF_false10.png",path ="/Users/raphaellemomal/these/R/images",
       height=8.5, width=9)

# comments:
# sur problème difficile:
#   des scores équivalents entre EMtree_freq (S=100) et glasso_cumsum.
#   des scores moins bons pour glasso_maxpen
#   EMtree_prob meilleur que EMtree_freq mais reste encore très proche de glasso_maxpen
#   AUC <0.6 pour n=10, AUC<0.85 pour n=60
#
# sur problème facile
#   écart de 5% entre EMtree_prob et glasso
#   meilleurs scores généraux
#   AUC<0.7 (+0.12 EMtree_prob) pour n=10, AUC>0.9 (+0.08 EMtree_prob) pour n=60
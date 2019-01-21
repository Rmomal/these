#####
# tableaux de comparaison TPFPTNFN avec deux niveaux de difficulté
####
# niveau facile : les graphs sont déjà faits
library(ROCR)
library(PLNmodels)
library(SpiecEasi)
library(MInt)
library(parallel)
source('~/simulations/codes/codegCoda.R')
source('~/simulations/codes/fonctions.R')
source('~/simulations/codes/FunctionsInference.R')
source('~/simulations/codes/Kit_pour_VEM_EM.R')
##########
type <-"cluster" ; variable <- "d";  valeur <-0.25 
eval_store_mint<-function(Y,covar,path){
  #  browser()
  data<-data_for_MInt(Y,covar,path)
  x <- data[["x"]]
  y <- data[["y"]]
  T1<-Sys.time()
  m <- mint(y,x,fmla = ~feature1 + feature2+ feature3)
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
covar<-data.frame(X1=round(runif(n)*5),X2=rnorm(n,0,2),
                  X3=(runif(n)))
saveRDS(covar,paste0(path,"mint_covar.rds"))
difficulty="easy"
##########
res<-lapply(1:100,function(nbgraph){
  omegaOriginal <- readRDS(paste0("/Users/raphaellemomal/simulations/Simu/PLN.2.0/",type,"/",
                                  variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds"))$omega
  edgesOrigin<-ifelse(abs(F_Sym2Vec(omegaOriginal))<1e-16,0,1)
  Y<-data_from_stored_graphs(type, variable, nbgraph, valeur, covar=covar,path,fixe=FALSE)
  saveRDS(Y,paste0(path,"Ymint_",variable,"_",valeur,".rds"))
 # Y<-readRDS(paste0(path,"Ymint_",variable,"_",valeur,".rds"))[[1]]
  Y=Y[[1]]
  m<- model.matrix(~X1+X2+X3,covar)# on choisi de mettre la constante dans PLN 
  p=ncol(Y)
  T1<-Sys.time()
  infMInt <- F_Sym2Vec(eval_store_mint(Y,covar,path)>1e-16)*1
  T2<-Sys.time()
  timeMInt<-difftime(T2,T1)
  
  T1<-Sys.time()
  pmat<-F_ResampleTreePLN(Y, X=m, O=matrix(0,nrow=100,ncol=20), v=0.8, B=150, maxIter=10, cond.tol=1e-12)
  infEM<- 1*(ifelse(colSums(ifelse(pmat<2/p,0,1))/100 >0.8,1,0))
  T2<-Sys.time()
  timeEM<-difftime(T2,T1)
  
  T1<-Sys.time()
  infgcoda <- F_Sym2Vec(gcoda(Y, counts=T, covar=covar)$opt.icov>1e-16)*1
  T2<-Sys.time()
  timegCoda<-difftime(T2,T1)
  
  U<-t(clr.matrix(Y,mar=1))

  model<-lm(U~m)
  T1<-Sys.time()
  infSpiec<-spiec.easi(model$residuals, method='mb', lambda.min.ratio=1e-2, nlambda=30,icov.select.params=list(rep.num=150))
  infSpiec<-F_Sym2Vec(as.matrix(infSpiec$refit[[1]]))
  T2<-Sys.time()
  timeSpiec<-difftime(T2,T1)

  res<-data.frame(cbind( rbind(c(table(infMINT,edgesOrigin)),
                    c(table(infEM,edgesOrigin)),
                    c(table(infgcoda,edgesOrigin)),
                    c(table(infSpiec,edgesOrigin))),method=c("MInt","EM","gCoda","SpiecEasi"),times=c(timeMInt,timeEM,timegCoda,timeSpiec),
              type=type,difficulty=difficulty,graph=nbgraph))
  
  return(res)
})


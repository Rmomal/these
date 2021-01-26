library(EMtree)
library(nestor)
library(PLNmodels)
library(tidyverse)
distance<-function(P1,B1,P2,B2){
  sum(F_Sym2Vec(log(B1/B2)*(P1-P2)/2))
}

n=200 ; p=100 ; S=200
types=c("erdos","scale-free","tree","cluster")
types="cluster"
some_fits=lapply(types, function(type){
  lapply(1:5, function(i){
    data=data_from_scratch(type,p=p, n=n,signed = TRUE,dens = 8/p, r=10, draw=TRUE)
    PLN_Y = PLNmodels::PLN(data$data~1, control=list(trace=0))
    EM=new_EMtree(PLN.Cor=PLN_Y, plot=FALSE)
    out=list(P=EM$edges_prob, B=EM$edges_weight)
  })
})
str(some_fits)
names(some_fits)=types
some_fits=do.call(rbind,some_fits)

mat_dist=sapply(1:20, function(i){
  sapply(1:20, function(j){
    distance( some_fits[[i]]$P,some_fits[[i]]$B,
              some_fits[[j]]$P,some_fits[[j]]$B)
  })
})

ggimage(mat_dist)
distance( some_fits[[1]]$P,some_fits[[1]]$B,
          some_fits[[18]]$P,some_fits[[18]]$B)

#-----------------
# Oak data
data.dir = '/Users/raphaellemomal/these/Data_Oak_remmoa/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))
theme_set(theme_bw())
# Parms
Y = as.matrix(Data$count)
colnames(Y)<-unlist(lapply(strsplit(colnames(Y),split="_"),function(x){paste0(toupper(x[1]),x[2])}))
colnames(Y)[48]<-"EA"
models=c("Null","Tree","Tree + D1 + D2 + D3")
# Models
null = PLNmodels::PLN(Y~1, control=list(trace=0))
EM_null=new_EMtree(PLN.Cor=null, plot=TRUE,
                   maxIter = 20, unif=TRUE,eps1 = 1e-3,
                   eps2 = 1e-3, cond.tol = 1e-6)
hist((EM_null$edges_prob), breaks=30)
nestor::ggimage(EM_null$edges_prob)
tree=PLNmodels::PLN(Y~Data$covariates$tree, control=list(trace=0))
EM_tree=new_EMtree(PLN.Cor=tree, plot=TRUE,
                   maxIter = 20, unif=TRUE,eps1 = 1e-3,
                   eps2 = 1e-3, cond.tol = 1e-6)
all=PLNmodels::PLN(Y~Data$covariates$tree+
                     Data$covariates$distTObase+
                     Data$covariates$distTOground+
                     Data$covariates$distTOtrunk, control=list(trace=0))
EM_all=new_EMtree(PLN.Cor=all, plot=TRUE,
                  maxIter = 20, unif=TRUE,eps1 = 1e-3,
                  eps2 = 1e-3, cond.tol = 1e-6)
hist(log(EM_null$edges_prob))

# distances
D1=distance(EM_null$edges_prob,EM_null$edges_weight,
         EM_tree$edges_prob,EM_tree$edges_weight)

D2=distance(EM_all$edges_prob,EM_all$edges_weight,
         EM_tree$edges_prob,EM_tree$edges_weight)

D3=distance(EM_all$edges_prob,EM_all$edges_weight,
         EM_null$edges_prob,EM_null$edges_weight)
D3
D2+D1

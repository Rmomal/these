# tests PLN sur REMMOA
library(PLNmodels)
library(EMtree)
library(tidyverse)
library(parallel)
library(tidygraph)
library(ggraph)
# data 
setwd(WorkDir <- "/Users/raphaellemomal/these/Data_Oak_remmoa/REMMOA")
ShapeDir <- paste(WorkDir, "shape", sep = "/")
OutDir <- paste(WorkDir, "output", sep = "/")

load(paste(OutDir, "20180620_PLN_REMMOANC.RData", sep = "/"))
colnames(X)[7:8]<-c("Xcart","Ycart")

# model naif
Y=as.matrix(Y)
f<-0.8
model<-PLN(Y ~ 1)
resample_output<-ResampleEMtree(Y,"1", S=20, maxIter=15,cond.tol=1e-8, cores=1)
str(resample_output)


df<-freq_selec(resample_output$Pmat,p=ncol(Y),f=f)
g_naif<-draw_network(df,"Naif", names=colnames(Y), layout="nicely")

# comptage median par détection pour chaque espèce
Y_corr=as.matrix(Y/N)
Y_corr[!is.finite(Y_corr)]<-0
coeff_troupeau<- apply( Y_corr, 2, function(x){
  median(x[x!=0])
})
Y_corr=lapply(seq_along(dim(Y_corr)[2]),function(x){
  Y_corr[,x]/coeff_troupeau[x]
})

# model données corrigé
resample_corr<-ResampleEMtree(Y_corr,"1", S=20, maxIter=15,cond.tol=1e-8, cores=1) # long ++
df_corr<-freq_selec(resample_output$Pmat,p=ncol(Y),f=f)
g_corr<-draw_network(df_corr,"comptages corrigés", pal="dodgerblue3", names=colnames(Y_corr), layout="nicely")

# model matrice offset coeff troupeau

Offset<-outer(rep(1,nrow(Y)),coeff_troupeau)
resample_offset<-ResampleEMtree(Y,"1", O=Offset,S=20, maxIter=15,cond.tol=1e-8, cores=1) # 20s par s
df_offset<-freq_selec(resample_output$Pmat,p=ncol(Y),f=f)
g_offset<-draw_network(df_offset,"offset", pal="dodgerblue3", names=colnames(Y), layout="nicely")

grid.arrange(g_naif,g_effort,g_offset, nrow=1, ncol=3)


# model covariables + offset d'effort
# null
mat_effort=outer(X$Effort,rep(1,32))
resample_effort<-ResampleEMtree(Y,"1", O=mat_effort,S=20, maxIter=15,cond.tol=1e-8, cores=1) 
df_effort<-freq_selec(resample_effort$Pmat,p=ncol(Y),f=f)
g_effort<-draw_network(df_effort,"effort", pal="dodgerblue3", names=colnames(Y_corr), layout="nicely")

t1<-Sys.time()
resample_covar_coord<-ResampleEMtree(Y,vec_covar=c("Xcart","Ycart","SEA_STATE","Depth","SUBJECTIVE"), 
                                     data_covar=data.frame(as.matrix(X)),O=mat_effort,S=20, maxIter=15,cond.tol=1e-8, cores=1) 
df_covar_coord<-freq_selec(resample_covar_coord$Pmat,p=ncol(Y),f=f)
g_covar_coord<-draw_network(df_covar_coord,"comptages corrigés", pal="dodgerblue3", names=colnames(Y_corr), layout="nicely")

t2<-Sys.time()
difftime(t2,t1)



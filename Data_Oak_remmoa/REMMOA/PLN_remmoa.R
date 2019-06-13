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
null_index=which(rowSums(Y)==0)
Y=as.matrix(Y[-null_index,])
X=as.matrix(X[-null_index,])
N=as.matrix(N[-null_index,])

###############
# essai prise en compte effet troupeau : 
#finalement la constante de PLN semble suffire (réseau invariant)
#################################################
# # comptage median par détection pour chaque espèce
# Y_corr=as.matrix(Y/N)
# Y_corr[!is.finite(Y_corr)]<-0
# coeff_troupeau<- apply( Y_corr, 2, function(x){
#   median(x[x!=0])
# })
# Y_corr=lapply(seq_along(dim(Y_corr)[2]),function(x){
#   Y_corr[,x]/coeff_troupeau[x]
# })
# 
# # model données corrigé
# resample_corr<-ResampleEMtree(Y_corr,"1", S=20, maxIter=15,cond.tol=1e-8, cores=1) # long ++
# df_corr<-freq_selec(resample_output$Pmat,p=ncol(Y),f=f)
# g_corr<-draw_network(df_corr,"comptages corrigés", pal="dodgerblue3", names=colnames(Y_corr), layout="nicely")
# 
# # model matrice offset coeff troupeau
# 
# Offset<-outer(rep(1,nrow(Y)),coeff_troupeau)
# resample_offset<-ResampleEMtree(Y,"1", O=Offset,S=20, maxIter=15,cond.tol=1e-8, cores=1) # 20s par s
# df_offset<-freq_selec(resample_output$Pmat,p=ncol(Y),f=f)
# g_offset<-draw_network(df_offset,"offset", pal="dodgerblue3", names=colnames(Y), layout="nicely")
#################################################

mat_effort=outer(X$Effort,rep(1,32))
X=data.frame(X)
X[,c(6:10,14:20)]<-apply(X[,c(6:10,14:20)], 2, function(x) as.numeric(x))


#################################################
# Etude des modèles ajustés de PLN selon les covariables
#################################################
attach(X)


get_model<-function(data, vec){
  t1<-Sys.time()
  string<-paste(deparse(substitute(data)), paste(vec, collapse=" + "), sep=" ~ ")
  formula<-as.formula(string)
  mat = as.matrix(lm(formula, x=T)$x)
  model<-PLN(data ~ -1+mat)
  t2<-Sys.time()
  print(difftime(t2,t1))
  return(model)
}
vec1=c("SEA_STATE","Depth","SUBJECTIVE") 

m0<-get_model(Y,"1") #22s
m1<-get_model(Y, "SEA_STATE") #30s
m2<-get_model(Y, "Depth") #26s, glm.fit: fitted rates numerically 0 occurred 
m3<-get_model(Y, "SUBJECTIVE") #28s
m4<-get_model(Y, vec1[-1])# 28s, In doTryCatch(return(expr), name, parentenv, handler) :display list redraw incomplete
m5<-get_model(Y, vec1[-2]) #56s, glm.fit: fitted rates numerically 0 occurred 
m6<-get_model(Y, vec1[-3]) #33s, glm.fit: fitted rates numerically 0 occurred 
m7<-get_model(Y, vec1) # 57s, glm.fit: fitted rates numerically 0 occurred 
m8<-get_model(Y, vec2) # 2.6min
m9<-get_model(Y, c(vec2,"Depth")) # 2.6min, glm.fit: fitted rates numerically 0 occurred 
m10<-get_model(Y, c(vec1,vec2)) # 56s, glm.fit: fitted rates numerically 0 occurred
save(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,file="models.RData")
# get characteristics and plot comparison
#load("models.RData")
Lmodels<-list(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
names(Lmodels)=0:10
criteria<-lapply(seq_along(Lmodels),function(x){
  c(Lmodels[[x]]$criteria,model=names(Lmodels)[x])
})
criteria=data.frame(apply(do.call(rbind,criteria),2,  as.numeric))
criteria %>% gather(crit,value,-model) %>% 
  ggplot(aes(reorder(as.factor(model),-value), value, color=crit))+
  geom_point()+facet_grid(crit~., scales = "free")+theme_light()+guides(color=FALSE)+
  labs(x="model",y="criteria value")
# Dans PLN on veut des grandes valeurs pour BIC, ICL et loglik
# les modèles 2, 4, et 6 comprenne "depth" qui semble essentielle
# les coordonnées spatiales seules (modèle 8) montrent de belles performances
# cependant le temps d'ajustement est très conséquent (2.6min contre 28 à 56s pour les autres)
# le gain en ICL de 10 par rapport à 7 est négligeable 
# depth+sea_state est meilleur que depth+spatial

#################################################
# Comparons les réseaux obtenus avec m0, m2 et m7 (pas d'ajustement spatial donc)
#################################################

t1<-Sys.time()
R0<-ResampleEMtree(Y,vec_covar = "1",O=NULL,S=20, maxIter=15,cond.tol=1e-8, cores=1) 
t2<-Sys.time()
R_effort<-ResampleEMtree(Y,vec_covar = "1",O=mat_effort,S=20, maxIter=15,cond.tol=1e-8, cores=1) 
t3<-Sys.time()
R_depth<-ResampleEMtree(Y,vec_covar = "Depth",O=mat_effort,S=20, maxIter=15,cond.tol=1e-8, cores=1) 
t4<-Sys.time()
R_full<-ResampleEMtree(Y,vec_covar = vec1,O=mat_effort,S=20, maxIter=15,cond.tol=1e-8, cores=1) 
t5<-Sys.time()
save(R0,R_effort,R_depth,R_full, file="Reseaux.RData")
timeR0=difftime(t2,t1) # 5min, alpha mean = 0.6959584
timeR_effort=difftime(t3,t2) # 5.3min alpha mean =  0.710248
timeR_depth=difftime(t4,t3) # 5.9min alpha mean = 0.7966224
timeR_full=difftime(t5,t4) # 15.3min alpha mean = 0.9731524

get_network<-function(resample_output,title, data=Y){
  df<-freq_selec(resample_output$Pmat,p=ncol(data),f=f)
  graph<-draw_network(df,title, pal="dodgerblue3", 
                              names=colnames(data), layout="nicely", size=3, curv=0.1)
  return(graph)
}
#load("Reseaux.RData")
graph0<-get_network(R0,"Naif")
graph_effort<-get_network(R_effort,"Effort")
graph_depth<-get_network(R_depth,"~ Depth")
graph_full<-get_network(R_full,"~ Depth+sea_state+subjective")

grid.arrange(graph0,graph_effort,graph_depth,graph_full,nrow=1, ncol=4)

# avec Xcart et Ycart error: inv_sympd(): matrix is singular or not positive definite
# Error in optimizer(unlist(par0), Y, X, O, w, ctrl) : nlopt failure
# refaire test avec Xcart et Ycart seuls

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


## nouvel offset : détectabilité
detec=data.frame(matrix(c(250,	"SMADEL",
                          285,	"BIGDEL",
                          340,	"SMAGLO",
                          300,	"BIGGLO",
                          310,	"BALSPP",
                          470,	"PHYMAC",
                          270,	"KOGIA",
                          310,	"ZIPHIUS",
                          260,	"DUGONG",
                          170,	"CHESPP",
                          430,	"RHITYP",
                          215,	"MARTEAU",
                          210,	"SHARKS",
                          170,	"DAURADE",
                          340,	"THON",
                          230,	"BILLFISH",
                          240,	"MANBIR",
                          210,	"MOBSPP",
                          145,	"DASSPP",
                          165,	"AETNAR",
                          160,	"RAIESP",
                          200,	"PHAETON",
                          200,	"GRETER",
                          200,	"BROTER",
                          200,	"GREPET",
                          200,	"BROPET",
                          200,	"LARNOV",
                          200,	"ANOSPP",
                          200,	"GYGALB",
                          200,	"FREGATE",
                          200,	"OCEANITE",
                          200,	"FOUS"), ncol=2, byrow=TRUE))
colnames(detec)=c("detec", "species" )
mat_detec=outer(rep(1,nrow(Y)), as.numeric(as.character(detec$detec)))
mat_offset=outer(X$Effort, 2*mat_detec[1,]/100)

## Réseaux totaux
t0<-Sys.time()
R_brut<-ResampleEMtree(Y,vec_covar = "1",S=20, maxIter=15,
                       cond.tol=1e-8, cores=1) 
t1<-Sys.time()
R_offset<-ResampleEMtree(Y,vec_covar = "1",O1=mat_offset,S=20, maxIter=15,
                         cond.tol=1e-8, cores=1) 
t2<-Sys.time() # 5.6 mins

R_depth<-ResampleEMtree(Y,vec_covar = "Depth",O1=mat_offset,S=20, maxIter=15,
                        cond.tol=1e-8, cores=1) 
t3<-Sys.time() # 6.4 mins
R_seastate<-ResampleEMtree(Y,vec_covar = "SEA_STATE",O1=mat_offset,S=20, maxIter=15,
                           cond.tol=1e-8, cores=1) 
t4<-Sys.time()#  10mins
R_depth_seastate<-ResampleEMtree(Y,vec_covar = c("Depth","SEA_STATE"),O1=mat_offset,S=20, maxIter=15,
                                 cond.tol=1e-8, cores=1) 
t5<-Sys.time()# 10mins
R_regions<-ResampleEMtree(Y,vec_covar = c("Region.Label"),O1=mat_offset,S=20, maxIter=15,
                          cond.tol=1e-8, cores=1) 
t6<-Sys.time()
save(R_brut,R_offset,R_depth,R_seastate,R_depth_seastate,R_regions, file="results/ReseauxEntiers.RData")


## Comparaisons
### Deux endroits de l'île
####  North
indicesNorth=which(X$Ycart > 7.76e6 & X$Xcart > 2e5 & X$Xcart<6e5)

Xfiltre=X[indicesNorth,]
YNorth=Y[indicesNorth,]
null_colnorth=which(colSums(YNorth)==0)
YNorth=Y[indicesNorth,-null_colnorth]
mat_offsetnorth=mat_offset[indicesNorth,-null_colnorth]
attach(Xfiltre)
t1<-Sys.time()
R_North<-ResampleEMtree(YNorth,vec_covar = "1",O1=mat_offsetnorth,S=20, maxIter=15,
                        cond.tol=1e-8, cores=1) 
t2<-Sys.time() # 1.9mins
t3<-Sys.time()
R_north_depth<-ResampleEMtree(YNorth,vec_covar = "Depth",O1=mat_offsetnorth,S=20, maxIter=15,
                              cond.tol=1e-8, cores=1) 
t4<-Sys.time() # 2.9mins
R_north_full<-ResampleEMtree(YNorth,vec_covar = c("SEA_STATE","Depth"),O1=mat_offsetnorth,S=20, maxIter=15,
                             cond.tol=1e-8, cores=1) 
t5<-Sys.time() # 2.8 mins

save(R_North,R_north_depth,R_north_full, file="ReseauxNorth.RData")

####  South
indicesSouth=which(X$Ycart < 7.55e6 & X$Ycart > 7.4e6 & X$Xcart > 5.5e5 & X$Xcart<9.5e5)
Xfiltre=X[indicesSouth,]
YSouth=Y[indicesSouth,]
null_colsouth=which(colSums(YSouth)==0)
YSouth=Y[indicesSouth,-null_colsouth]
mat_offsetsouth=mat_offset[indicesSouth,-null_colsouth]
attach(Xfiltre)
t6<-Sys.time()
R_South<-ResampleEMtree(YSouth,vec_covar = "1",O1=mat_offsetsouth,S=20, maxIter=15,
                        cond.tol=1e-8, cores=1) 
t7<-Sys.time() # 41 s
R_south_depth<-ResampleEMtree(YSouth,vec_covar = "Depth",O1=mat_offsetsouth,S=20, maxIter=15,
                              cond.tol=1e-8, cores=1) 
t8<-Sys.time() # 46 s
R_south_full<-ResampleEMtree(YSouth,vec_covar = c("SEA_STATE","Depth"),O1=mat_offsetsouth,S=20, maxIter=15,
                             cond.tol=1e-8, cores=1) 
t9<-Sys.time() # 57s

save(R_South,R_south_depth,R_south_full, file="ReseauxSouth.RData")

################################
# Comparaison par Region.Label #
################################

region_data<-function(region, vec="1", S=20){
  #data : filtrer pour la région et supression des colonnes vides

  Yregion=Y %>% as_tibble() %>% mutate(reg=X$Region.Label) %>% 
    filter(reg==region) %>% dplyr::select(-reg)
  offset=mat_offset %>% as_tibble() %>% mutate(reg=X$Region.Label) %>% 
    filter(reg==region) %>% dplyr::select(-reg)
  cols=which(colSums(Yregion)==0)
  if(length(cols)!=0){
    Yregion=Yregion %>% dplyr::select(-cols)
    offset=offset%>% dplyr::select(-cols)
  }
  offset=as.matrix(offset)
  dim(Yregion)
  return(list(Y=Yregion, O=offset))
}

region_network<-function(region, vec="1", S=20, data){
  #data : filtrer pour la région et supression des colonnes vides

  data=data[[region]]
  Yregion=data$Y

  offset=data$O
  Xfiltre=X %>% as_tibble() %>% filter(Region.Label==region)
  attach(Xfiltre)
  t1<-Sys.time()
  R<-ResampleEMtree(Yregion,vec_covar = vec,O1=offset,S=S, maxIter=30,
                    cond.tol=1e-8, cores=1)
  t2<-Sys.time()
  cat("\n",region, difftime(t2,t1),"\n")
  return(R)
}
DataRegion<-lapply(unique(X$Region.Label), function(x) region_data(x))
names(DataRegion)=unique(X$Region.Label)
save(DataRegion, file="results/DataRegion.RData")

ReseauxRegion<-lapply(unique(X$Region.Label), function(x) region_network(region=x, data=DataRegion))
names(ReseauxRegion)=unique(X$Region.Label)
save(ReseauxRegion,file="results/ReseauxRegion.RData")


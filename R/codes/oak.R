rm(list=ls()); par(pch=20);
#devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
library(PLNmodels);
library(igraph)
library(RColorBrewer)
library(tidyverse)
library(parallel)
source('/Users/raphaellemomal/simulations/codes/FunctionsMatVec.R')
source('/Users/raphaellemomal/simulations/codes/FunctionsTree.R')
source('/Users/raphaellemomal/simulations/codes/FunctionsInference.R')
#source('TreeMixture-mac.R')
#source('/Users/raphaellemomal/simulations/codes/TreeMixture-RML.R')
source('/Users/raphaellemomal/simulations/codes/fonctions.R')
# Data
data.dir = '/Users/raphaellemomal/these/Data/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))
theme_set(theme_bw())
# Parms
Y = as.matrix(Data$count);
O = Data$offset; X = Data$covariates
# Selection des especes
# YO = Y/O
# Rank = rank(colSums(YO))
# plot(cumsum(sort(colSums(YO))))
# Seuil = 20
# Ycum = colSums(Y); Order = order(Ycum)
# #plot(cumsum(Ycum[Order]), col = 1+(Rank[Order]<Seuil))
# 
# Y = Y[, Rank > Seuil];
# O = O[, Rank >Seuil];
n = nrow(Y); p = ncol(Y)
colnames(Y)<-unlist(lapply(strsplit(colnames(Y),split="_"),function(x){paste0(toupper(x[1]),x[2])}))
colnames(Y)[48]<-"EA"
#################################################################################
#################################################################################

# Inférences
PLN.vide<-PLN(Y ~ 1 + offset(log(O)))
PLN.tree<-PLN(Y ~ 1 + X$tree+  offset(log(O)))
PLN.orient<-PLN(Y ~ 1 + X$orientation+  offset(log(O)))
PLN.dist1<-PLN(Y ~ 1 + X$distTOground+  offset(log(O)))
PLN.dist2<-PLN(Y ~ 1 + X$distTObase+  offset(log(O)))
PLN.dist3<-PLN(Y ~ 1 + X$distTOtrunk+  offset(log(O)))
PLN.distOrient<-PLN(Y ~ 1 + X$distTOground+ X$orientation+  offset(log(O)))
PLN.infect<-PLN(Y ~ 1 +  X$pmInfection+  offset(log(O)))
PLN.infectDist<-PLN(Y ~ 1 +  X$pmInfection + X$distTOground+  offset(log(O)))
PLN.infectAllDist<-PLN(Y ~ 1 +  X$pmInfection +X$distTObase+X$distTOtrunk+X$distTOground+   offset(log(O)))
PLN.infectTree<-PLN(Y ~ 1 +  X$pmInfection + X$tree+  offset(log(O)))
PLN.infectOrient<-PLN(Y ~ 1 +  X$pmInfection + X$orientation+  offset(log(O)))
PLN.infectTreeOrient<-PLN(Y ~ 1 +  X$pmInfection + X$tree+X$orientation+  offset(log(O)))
PLN.Alldist<-PLN(Y ~ 1 +  X$distTObase+X$distTOtrunk+X$distTOground+ offset(log(O)))
PLN.treeOrient<-PLN(Y ~ 1 + X$tree+ X$orientation+  offset(log(O)))
PLN.treeAlldist<-PLN(Y ~ 1 +  X$tree+ X$distTObase+X$distTOtrunk+X$distTOground+ offset(log(O)))

PLN.treeDist<-PLN(Y ~ 1 +
                    X$tree+
                    X$distTOground+
                    offset(log(O)))
PLN.treeallDist<-PLN(Y ~ 1 +
                       X$tree+
                       X$distTObase+X$distTOtrunk+X$distTOground+
                       offset(log(O)))
PLN.allcovariates<-PLN(Y ~ 1 +
                         X$tree+
                         X$distTObase+X$distTOtrunk+X$distTOground+
                         X$pmInfection+X$orientation+ 
                         offset(log(O)))

datamodels<-data.frame(cbind(rbind(PLN.vide$criteria,
                                   PLN.tree$criteria,
                                   PLN.orient$criteria,
                                   PLN.dist1$criteria,
                                   PLN.dist2$criteria,
                                   PLN.dist3$criteria,
                                   PLN.treeDist$criteria,
                                   PLN.treeOrient$criteria,
                                   PLN.Alldist$criteria,
                                   PLN.treeDist$criteria,
                                   PLN.allcovariates$criteria,
                                   PLN.infect$criteria,
                                   PLN.infectDist$criteria,
                                   PLN.infectAllDist$criteria,
                                   PLN.infectTree$criteria,
                                   PLN.infectOrient$criteria,
                                   PLN.infectTreeOrient$criteria,
                                   PLN.treeAlldist$criteria,
                                   PLN.distOrient$criteria),
                             model=c("vide","tree","orient","dist1","dist2","dist3","treeDist",
                                     "treeOrient","alldist","treeDist","allcovariates",
                                     "infect","infectDist","infectAllDist", "infectTree", "infectOrient",
                                     "infectTreeOrient","TreeAllDist","distOrient")))

datamodels[,1:5]<-apply(datamodels[,1:5],2,function(x) as.numeric(as.character(x)))
p1<-ggplot(datamodels, aes(x = reorder(model, loglik), y = loglik)) +
  geom_point()+
  # geom_point(aes(y=ICL),color="red")+
  # geom_point(aes(y=BIC),color="blue")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2<-ggplot(datamodels, aes(x = reorder(model, -ICL), y = ICL)) +
  geom_point(color="red")+
  # geom_point(aes(y=loglik))+
  # geom_point(aes(y=BIC),color="blue")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p3<-ggplot(datamodels, aes(x = reorder(model, -BIC), y = BIC)) +
  geom_point(color="blue")+
  # geom_point(aes(y=ICL),color="red")+
  # geom_point(aes(y=loglik))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(p1,p2,p3,nrow=3,ncol=1)
png(filename="~/these/fig.png", width=4, height=7, units="in",res=300)
grid.arrange(p1,p2,p3,nrow=3,ncol=1)
dev.off()
# datamodels<-datamodels %>% 
#   gather(key, value, -"model")
# ggplot(datamodels,aes(model))+
#   geom_col(aes(y=value,fill=key))

Z.vide = PLN.vide$model_par$Sigma
Z.tree=PLN.tree$model_par$Sigma
Z.treeDist=PLN.treeAlldist$model_par$Sigma



# trace de la diff, de la loglik. affichages des histogrammes de beta et de leur summary (save des summary)
# 1. trajectoires vraisemblances EM
track_convergence<-function(inf.vide,title){
  loglik<-data.frame(logL=inf.vide$L,diff=inf.vide$diff_beta)
  
  p1<-loglik %>% 
    mutate(index=1:nrow(loglik)) %>% 
    ggplot(aes(x=index,y=logL))+
    geom_point(color="deepskyblue")+
    geom_line()+
    theme_bw()
  p2<-loglik %>% 
    mutate(index=1:nrow(loglik)) %>% 
    ggplot(aes(x=index,y=diff))+
    geom_point(color="deepskyblue")+
    geom_line()+
    theme_bw()
  p3<-loglik %>% 
    mutate(index=1:nrow(loglik)) %>% 
    ggplot(aes(x=index,y=diff))+
    geom_point(color="deepskyblue")+
    geom_line()+
    coord_cartesian(ylim = c(0,2e-4))+
    geom_abline(intercept = 1e-4,slope=0,color="deepskyblue",size=1.3)+
    geom_abline(intercept = 1e-6,slope=0,color="deepskyblue3")+
    theme_bw()+
    labs(y="diff zoom1")
  p4<-loglik %>% 
    mutate(index=1:nrow(loglik)) %>% 
    ggplot(aes(x=index,y=diff))+
    geom_point(color="deepskyblue")+
    geom_line()+
    coord_cartesian(ylim = c(0,2e-6))+
    geom_abline(intercept = 1e-6,slope=0,color="deepskyblue3", size=1.3)+
    theme_bw()+
    labs(y="diff zoom2")
  grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2,top=title)
}

T1<-Sys.time()
inf.vide<-TreeGGM(cov2cor(Z.vide),print=FALSE,step="FALSE",maxIter = 150, n=nrow(Y) )
# itermax=154, 13.3 min // alpha=0.1034, 6.13 mins itermax=22 mais NAN // 1.34e_17 // NAN et Error in while ((diff > eps) & (iter < maxIter)) { : 
# missing value where TRUE/FALSE needed
# tolérance trop faible ?
#// 5.7 min, itermax=24
T2<-Sys.time()
difftime(T2,T1)
T1<-Sys.time()
inf.sansnorm<-TreeGGM(cov2cor(Z.tree),print=FALSE,step="FALSE",maxIter = 150,n=nrow(Y),
                      cond.tol=1e-12)
# itermax = 243 20min // alpha = 0.224, 2.56 min, itermax=20, et again NaN 1.497385e-17 NaN 
# Error in while ((diff > eps) & (iter < maxIter)) { : 
#     missing value where TRUE/FALSE needed
#// 2.6 mins, itermax=21
T2<-Sys.time()
difftime(T2,T1)

#inf.betamin1, 40, 500.26
# inf.betamin2, 150, 500.26
# inf.norm1, 37, 488.6
# inf.norm2, 42, 489

# impact sur les valeurs de proba ?
# attendu : différence entre betalin1 et betamin2, mais pas entre norm1 et norm2
save(inf.betamin1,inf.betamin2,inf.norm1,inf.norm2,inf.norm3,file="testsBetas")
summary(diff(F_Sym2Vec(inf.betamin1$P)-F_Sym2Vec(inf.betamin2$P)))
summary(diff(F_Sym2Vec(inf.norm1$P)-F_Sym2Vec(inf.norm2$P)))
summary(diff(F_Sym2Vec(inf.norm3$P)-F_Sym2Vec(inf.norm2$P)))
#############################
data.frame(inf.tree$L )%>% 
  mutate(index=1:length(inf.tree$L)) %>% 
  ggplot(aes(index,inf.tree.L))+
  geom_point(color="deepskyblue")+
  geom_line() + 
  facet_zoom(y= inf.tree.L >488.6)+
  theme_bw()


T1<-Sys.time()
inf.treeDist<-TreeGGM(cov2cor(Z.treeDist),print=FALSE,step="FALSE",maxIter = 150,n=nrow(Y))
# itermax=183, 13 min // alpha=0.25, 6.52min, itermax=29, pas d'erreur.
T2<-Sys.time()
difftime(T2,T1)

track_convergence(inf.vide,"Vide")
track_convergence(inf.tree,"Tree")
track_convergence(inf.treeDist,"Tree and distances")

save(inf.tree,inf.vide, inf.treeDist,file="Res_inf_maxlikelihood.RData")
############################################################
#compare beta1 beta2 inf finale
vide.beta1<-readRDS("seuil1vide.rds")
vide.beta2<-readRDS("seuil2vide.rds")
vide.finale<-inf.vide$P
plot(F_Sym2Vec(vide.beta1),F_Sym2Vec(vide.finale))
# 4. Coudes selon les modèles
compar_coudes<-data.frame(vide=sort(F_Sym2Vec(inf.vide$probaCond)), tree=sort(F_Sym2Vec(inf.tree$probaCond)),
                          treeDist=sort(F_Sym2Vec(inf.treeDist$probaCond)))
compar_coudes %>% 
  mutate(index=1:nrow(.)) %>% 
  gather(model,Edgeproba,-index) %>% 
  ggplot(aes(x=index,y=Edgeproba,color=model))+
  geom_point(size=0.3)+
  coord_cartesian(xlim=c(3800,nrow(compar_coudes)),ylim=c(0,0.5)) # !! zoom !!

# 5. jitter des degrés selon les modèles
maxDegrees<-1*(colSums(inf.vide$probaCond) %in% tail(sort(colSums(inf.vide$probaCond))))+
  2*(colSums(inf.tree$probaCond) %in% tail(sort(colSums(inf.tree$probaCond))))+
  3*(colSums(inf.treeDist$probaCond) %in% tail(sort(colSums(inf.treeDist$probaCond))))

compar_degres<-data.frame(vide=colSums(inf.vide$probaCond), tree=colSums(inf.tree$probaCond),
                          treeDist=colSums(inf.treeDist$probaCond),color=rep(c(1,2,1,3,1),c(13,1,29,1,50)),
                          color2= maxDegrees)
# color according to EA/F19/else
compar_degres %>% 
  mutate(color=as.factor(color),color2=as.factor(color2),index=1:nrow(compar_degres)) %>% 
  gather(model,degree,-color,-color2,-index) %>% 
  ggplot(aes(x=reorder(model,degree),y=degree,color=color, size=color))+
  scale_color_manual(values=c("lightslategrey","royalblue1","chartreuse3"),breaks=c(2,3,1),
                     labels=c("F19","EA","else"))+
  scale_size_manual(values=c(0.5,2.5,2.5),breaks=c(2,3,1),
                    labels=c("F19","EA","else"))+
  geom_jitter(width=0.2)+
  labs(x="",y="nodes degree estimates")

# color according to max degree
compar_degres %>% 
  mutate(color=as.factor(color),color2=as.factor(color2),index=1:nrow(compar_degres)) %>% 
  gather(model,degree,-color,-color2,-index) %>% 
  ggplot(aes(x=reorder(model,degree),y=degree,color=color2, size=color2))+
  scale_color_manual(name="max in:",values=c("lightslategrey","hotpink2","deepskyblue2","orange2"),breaks=c(6,5,1,0),
                     labels=c("all models","models with \ncovariates","empty","else"))+
  scale_size_manual(name="max in:",values=c(0.5,2.5,2.5,2.5),breaks=c(6,5,1,0),
                    labels=c("all models","models with \ncovariates","empty","else"))+
  geom_jitter(width=0.2)+
  labs(x="",y="nodes degree estimates")

# 6. ridges
compar_degres %>% 
  mutate(color=as.factor(color),index=1:nrow(compar_degres)) %>% 
  gather(model,degree,-color,-index) %>% 
  ggplot(aes(y=reorder(model,degree),x=degree,fill=factor(..quantile..))) +
  stat_density_ridges(scale=1.5,geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      quantiles = c(0.025, 0.975)) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
  )





############################################################
T1<-Sys.time()
inf.allcovariates<-TreeGGM(cov2cor(Z.allcovariates),print=FALSE,
                           step="FALSE",maxIter = 150)# 8 min avec eps=1e-6, maxier=79
T2<-Sys.time()
difftime(T2,T1)

saveRDS(inf.vide,"inf_vide.rds")
saveRDS(inf.tree,"inf_tree.rds")
saveRDS(inf.treeDist,"inf_treeDist.rds")
saveRDS(inf.allcovariates,"inf_allcovariates.rds")
inf.vide<-readRDS("inf_vide.rds")
inf.tree<-readRDS("inf_tree.rds")
inf.allcovariates<-readRDS("inf_allcovariates.rds")
scores<-data.frame(vide=c(inf.vide$P[upper.tri(inf.vide$P,diag = FALSE)]),
                   tree=c(inf.tree$P[upper.tri(inf.tree$P,diag = FALSE)]),
                   treeDist=c(inf.treeDist$P[upper.tri(inf.treeDist$P,diag = FALSE)]),
                   allcov=c(inf.allcovariates$P[upper.tri(inf.allcovariates$P,diag = FALSE)]),
                   index=1:length(c(inf.vide$P[upper.tri(inf.vide$P,diag = FALSE)])))
scores<-scores %>% mutate(color=vide) %>% 
  gather(model,proba,-"index",-"color") 
ggplot(scores,aes(model,proba,group=index,color=color))+
  geom_line()+ 
  scale_x_discrete(limits=c("vide","tree","treeDist","allcov"))

#################################################################################
#################################################################################
# Stability selection
# F_ResampleTreePLN <- function(Y, X, O, v = 0.8, B = 1e2){
#   n = nrow(Y)
#   p = ncol(Y)
#   P = p * (p - 1) / 2
#   V = round(v * n)
#   # Pmat = matrix(0, B, P)
#   
#   # browser()
#   # if(is.null(vecVar)){
#   #    fmla <- as.formula(paste("Y.sample ~ -1 + offset(log(O.sample))")) 
#   #  }else{
#   #    fmla <- as.formula(paste("Y.sample ~ -1 + offset(log(O.sample))+", 
#   #                             paste(paste0("X.sample$",vecVar), collapse= "+")))
#   #  }
#   res<-mclapply(1:B,function(x){
#     sample = sample(1:n, V, replace = F)
#     Y.sample = Y[sample,]
#     X.sample = X[sample,]
#     O.sample = O[sample,]
#     # PLN.sample = switch(model,"1"=PLN(Y.sample ~ 1 + offset(log(O.sample))),
#     #                     "2"=PLN(Y.sample ~ 1 + X.sample$tree+  offset(log(O.sample))),
#     #                     "3"=PLN(Y.sample ~ 1 + X.sample$tree+ X.sample$orientation+  offset(log(O.sample))),
#     #                     "4"=PLN(Y.sample ~ 1 + X.sample$tree+ X.sample$orientation+ X.sample$distTObase+  offset(log(O.sample))))
#     try({PLN.sample=PLN(Y.sample ~ 1 + X$tree+ offset(log(O.sample)))},silent=TRUE)
#     if(exists("PLN.sample")){
#       Sigma.sample = PLN.sample$model_par$Sigma
#     inf<-TreeGGM(cov2cor(Sigma.sample), "FALSE", FALSE, maxIter = 150)
#     res<-list(F_Sym2Vec(inf$P), F_Sym2Vec(inf$probaCond),inf$L)
#     }else{
#       res<-list(NA,NA,NA)
#     }
#     
#     return(res)
#   }, mc.cores=3)
#   return(res)
# }
if(Tree.res){
  #  Stab.sel = list()
  covar = list(); 
  covar$null = matrix(1, n, 1); 
  covar$tree = as.matrix(lm(Y ~ X$tree, x=T)$x)
  covar$treeDist = as.matrix(lm(Y ~ X$tree+ X$distTObase+ X$distTOtrunk+X$distTOground, x=T)$x)
  T1<-Sys.time()
  
  StabselOak<- lapply(seq_along(covar),function(x){
    cat(" ",names(covar)[x])
    co=covar[[x]]
    ResampleEMtree(Y,vec_covar =NULL, covariate=co,log(O), S=100,maxIter=300, cond.tol=1e-8, cores=1)
  })
  
  T2<-Sys.time()
  difftime(T2,T1)
  save(StabselOak, file = paste0(data.dir, data.name, '-StabselOak_allSpecies.Rdata'))
}

load(file = paste0(data.dir, data.name, '-StabSeliter_null.Rdata'))
mat<-matrix(do.call(c,lapply(Stab.sel[[1]], function(x){x[[2]]})),nrow =B.resample)

# Stab.sel est une liste de taille 5, qui récupère listPmat, les temps, L, alpha et maxIter
# listPmat est de taille 3, chaque elmt est la Pmat à l'itération demandée (1,3, et 10 ou max)

# STEP 1 : mettre en jolie forme pour avoir envie de travailler. lol.
indices=c(1,3,10)
freqs<-do.call(rbind,lapply(seq_along(Stab.sel[[1]]), function(x) {
  proba<-Stab.sel[[1]][[x]]
  colScores<-colMeans(1 * (proba> 2 / p))
  data<-tibble(meanProb=colScores,Maxiter=indices[x],edge=seq_along(colScores))
  return(data)
}))
maxItAlpha<-tibble(Alpha=Stab.sel[[4]],Maxiter=Stab.sel[[5]])

indices_inf<-which(Stab.sel[[5]]!=10)
times<-tibble(time=Stab.sel[[2]], iter=rep(c(1,3,10),40)[-3*indices_inf])
times<-  times %>% mutate(time=ifelse(iter==10,time*60,time))
times %>% group_by(iter) %>% summarise(med=median(time))
# likelihoods<-do.call(rbind,lapply(seq_along(Stab.sel), function(x) {
#   if(x>1){
#     if(x==4){
#       like <- Stab.sel[[x]][[3]]
#       indices<-seq(329*4+1,330*4-1)
#       
#       L<-data.frame(matrix(like[-indices],ncol=x, byrow = TRUE))
#       
#     }else{L<-data.frame(matrix(Stab.sel[[x]][[3]],ncol=x, byrow = TRUE))
#     }
#     
#     colnames(L)<-1:x
#     data<-L %>% mutate(ligne=1:nrow(L)) %>% gather(key=Iter,value=likelihood,-ligne)
#     data$Maxiter=x
#   }else{
#     L<-Stab.sel[[x]][[3]] 
#     data<-tibble(Iter=x,likelihood=L, Maxiter=x, ligne=1:length(L))
#   }
#   return(data)
#}))
iters<-c()
sapply(Stab.sel[[5]],function(x){
  iters<<-c(iters,1:x)
})
likelihoods<-tibble(L=Stab.sel[[3]],iter=iters) 
likelihoods%>% ggplot(aes(x=as.factor(iter),y=L))+
  geom_boxplot()

# STEP2 etude
timesAlpha<-timesAlpha %>% mutate(Time=as.numeric(Time), Maxiter=as.factor(Maxiter)) %>% filter(Time<3000) 
timesAlpha%>% 
  ggplot(aes(x=Time, y=Maxiter, fill=Alpha))+
  geom_density_ridges()

likelihoods<-likelihoods %>% mutate(Iter=as.numeric(Iter)) # % variation de 3 à 4
likelihoods %>% filter(Maxiter==4,ligne==2) %>% 
  ggplot(aes(x=Iter,y=likelihood))+
  geom_point()

likelihoods %>% arrange(ligne) %>% filter(Maxiter==4) %>% spread(Iter,likelihood) %>% 
  mutate(pctdiff3=100*(`4`-`3`)/(`4`-`1`),pctdiff2=100*(`4`-`2`)/(`4`-`1`)) %>% 
  select(pctdiff3,pctdiff2) %>% gather(diff,value) %>% 
  ggplot(aes(x=value,color=diff,fill=diff))+
  geom_density(alpha=0.8)



freqs %>% arrange(edge) %>% spread(Maxiter,meanProb) %>% 
  gather(iter,freq,-edge,-`10`) %>% 
  ggplot(aes(x=`10`, y=freq, color=iter))+geom_point()+
  geom_abline()+
  coord_cartesian(xlim=c(0.4,1),ylim=c(0.4,1))+ geom_hline(yintercept=0.8)+geom_vline(xintercept = 0.8)


freqlarge<-freqs %>% arrange(edge) %>% spread(Maxiter,meanProb)
freqlarge<-1*(freqlarge[,-1]>0.8)
table(freqlarge[,2],freqlarge[,3]) # en rate 1 et ajoute 3 
table(third=freqlarge[,2],last=freqlarge[,3]) # en rate 5 et en ajoute 83

############################################@

load(paste0(data.dir, data.name, '-StabSel.Rdata'))

library(ggraph)
library(tidygraph)
M<-4
p<-ncol(Y)
p<-94
f=0.95
Stabdata=StabselOak
models=c("Null","Tree","Tree + D1 + D2 + D3")
countEdges<-function(Stabdata,f,p){
  mat<-data.frame(freq_selec_list(Stabdata,1,p,f))
  # models=names(covar)
  allNets<-tibble(P = list(mat), models =models )  %>%
    mutate(P=map( seq_along(P), function(x) {
      df<-freq_selec_list(Stabdata,x,p,f)
      df[lower.tri(df, diag = TRUE)]<-0
      df<-data.frame(df)
      colnames(df)<-1:ncol(df)
      df
    })) %>%
    mutate(P = map(P,~rownames_to_column(.) %>%
                     gather(key, value , -rowname))) %>%
    unnest()
  allNets<-allNets[,c(3,2,1,4)]
  # res<-allNets %>% group_by(models) %>% summarise(sum=sum(value))
  # res$seuil=f
  # return(res)
  return(allNets)
}

compar_graphs(oakNets)
countsEdges<-do.call(rbind,lapply(seq(0,1,0.01), function(x) countEdges(StabselOak, x)))
plot<-countsEdges %>% as_tibble() %>% mutate(models=rename_factor(as.factor(models), 
                                                                  "null"="~ 1","tree"="~ tree","treeDist"="~ tree + D1 + D2 + D3"
)) %>% 
  ggplot(aes(seuil, sum, color=models, shape=models))+ theme_minimal()+geom_line(alpha=0.4)+
  labs(x="Selection threshold", y="Quantity of edges")+ theme(axis.text=element_text(size=12),
                                                              axis.title=element_text(size=12), legend.position="bottom",
                                                              legend.text = element_text(size=12))+
  scale_colour_manual("",values=pal[c(2,1,4,3)])+ scale_shape_manual("",values=c(15,16,17))
legend=g_legend(plot+geom_point(size=2.2,alpha=0.8))

plot=plot+geom_point(alpha=0.8)+guides(color=FALSE, shape=FALSE)

plot<-grid.arrange(plot,legend, nrow=2, ncol=1,heights=c(9,1))


ggsave("QuantEdgeThreshold.png",plot=plot, width=7.3, height=7.3)



oak_plot<-function(data,seuil,models, plot=FALSE){
  p=nrow(Y)
  allNets=countEdges(data,seuil,p)
  groupes<-as.factor(substr(colnames(Y),1,1))
  set_graph_style()
  mods<-models
  spliT<-data.frame(allNets) %>% 
    split(allNets$models) %>% 
    tibble(P=map(.,function(x){
      model<-x$models[1]

      
      res<- as_tbl_graph(x, directed=FALSE) %>%
        activate(edges) %>%
        filter(value !=0) %>% 
        activate(nodes) %>% 
        mutate( groupe=groupes, btw=centrality_betweenness(normalized = TRUE),imp=centrality_degree(),
                clo=centrality_closeness(normalized = TRUE),
                lab=colnames(Y),auth=centrality_authority(),
                keyplayer = node_is_keyplayer(k=6), model=model,EA=(name%in%c("48")),
                boolbtw=(btw>sort(btw, decreasing = TRUE)[7]), boolimp=(imp>=sort(imp, decreasing = TRUE)[7])) 
      nodes =res %N>% as.tibble()
      
      print(nodes$btw[48])
    print(pemp(nodes$btw[48],nodes$btw, discrete = FALSE))
  #    print(nodes$clo[48])
 
      print(nodes$auth[48])
      print(pemp(nodes$auth[48],nodes$auth, discrete = FALSE))
    #  print(nodes$imp[48])
      res<- res %>% activate(edges) %>% 
        mutate(neibs=edge_is_incident(which(.N()$EA)), model=model) %>% 
        activate(nodes) %>% 
        mutate(neibEA=name%in%unique(c(.E()$from[which(.E()$neibs)],.E()$to[which(.E()$neibs)]))) %>% 
        mutate(label=ifelse(boolbtw|EA|neibEA,lab,""))# %>% 
      #filter(importance!=0)
      
      return(res)
      
    }))
  datatest=spliT %>% 
    mutate( neibors=map(P,function(x){
      edges =x %N>% as.tibble()
      res=edges[edges$neibEA,c(9,14)]
      return(res) }))
  join=full_join(datatest$neibors[[1]],datatest$neibors[[2]], by="label")
  joinjoin=full_join(join,datatest$neibors[[3]], by="label")
  joinjoin=full_join(joinjoin,data.frame(label=colnames(Y)), by="label")
  
  if(plot){
    pal <- viridisLite::viridis(5, option = "C")
    pal2<-c("gray15","#9de1e1")
    coords<-create_layout(spliT$P[3][[1]],layout="star",center=48)[,c("x","y")]
    plot=spliT$P[1][[1]] %>%
      bind_graphs(spliT$P[2][[1]] )%>%
      bind_graphs(spliT$P[3][[1]] )%>%
      activate(nodes) %>% 
      mutate(model=factor(model,levels=mods), x=rep(coords$x,3),y=rep(coords$y,3)) %>% 
      ggraph(layout="auto")+
      geom_edge_arc(aes(color=model,alpha=neibs),curvature=0.1,show.legend=FALSE)+ 
      scale_edge_alpha_manual(values=c(0.2,1))+
      geom_node_point(aes(color=groupe, size=boolbtw), show.legend=FALSE)+
      scale_edge_colour_manual("model",values=pal[c(2,1,4,3)], labels=mods)+
      scale_color_manual(values=c("indianred1","steelblue4","orange2"))+
      scale_size_manual(values=c(1.5,6))+
      geom_node_text(aes(label = label),color="black")+
      facet_nodes(~model, scales="free")+
      th_foreground(border=FALSE)+
      theme(strip.background = element_rect(fill="white",color="white"),
            strip.text = element_text(color="black",size=14))
    res=list(plot,joinjoin)
  }
  datajak=data.frame(label=jakuschlabels, jak="jak")
  final_join= full_join(joinjoin,datajak, by="label")
  final_join=final_join[,c(2,1,3,4,5)]
  final_join= final_join %>% mutate(EMtree=ifelse(is.na(model.x) & is.na(model.y) & is.na(model),0,1),
                                    jak=ifelse(is.na(jak),0,1))
  
  
  compare=data.frame(table(EMtree=final_join$EMtree,jak=final_join$jak))
  res=spread(unite(compare,"EMtree_jak",c(EMtree,jak),sep=""), EMtree_jak, Freq)
  
  
  # num=length(intersect( jakuschlabels, joinjoin$label))
  # pct1=num/length(joinjoin$label)
  # pct2=num/length(jakuschlabels)
  # browser()
  #    res=data.frame("Emtree"=pct1,"jakusch"=pct2, "count"=length(joinjoin$label))
  return(res)
}

oak_plot(StabselOak,0.9, models)
datapct=lapply(seq(0,0.99,0.02),function(x){
  return(cbind(oak_plot(StabselOak,x, models),"seuil"=x))
})
datapct=do.call(rbind,datapct)
dataplot=datapct %>% as_tibble() %>% 
  mutate(sumEMtree=`10`+`11`, sumjak=`01`+`11`, FDR=`10`/sumEMtree, PPVEMtree=`11`/sumEMtree,
         PPVjak=`11`/sumjak, =) 



plot=dataplot %>%  select(-PPVEMtree,-PPVjak,-FDR,-sumjak,-sumEMtree) %>% gather(type, value, -seuil) %>% 
  ggplot( aes(seuil, value, col=type))+geom_point()+geom_line()+theme_minimal()+
  labs(x="Selection threshold",y="Number of OTUs",color = "EMtree vs.\nJakuschkin \net al.:")+
  geom_vline(xintercept = 0.9, linetype="dashed")+scale_color_brewer(palette="Dark2")
ggsave("curves_OTUS.png",plot=plot,height=5,width=7.5 )


brewer.pal(8,"Dark2")

plot2=dataplot %>%  select(PPVEMtree,seuil) %>% gather(type, value, -seuil) %>% 
  ggplot( aes(seuil, value, col=type))+geom_point()+geom_line()+theme_minimal()+
  scale_color_manual(name="",labels="(OTUs in common)/\n(EMtree discoveries)", values="#66A61E")+
  geom_vline(xintercept = 0.9, linetype="dashed")+
  labs(x="Selection threshold",y="%")

plot3=grid.arrange(plot,plot2, nrow=1, ncol=2)
ggsave("propOTUS.png",plot=plot2,height=4,width=6)



ggsave("Oaknet_90_btw.png", width=11, height=4.5)
################
# compare tables accross models
load("/Users/raphaellemomal/these/Data/Oaks-CVacher/oaks-StabSel.Rdata")
null1<-Stab.sel[[1]][[1]][[1]] # comparaison des proba non seuillées pour 1 sous échantillon
null2<-Stab.sel[[1]][[2]][[1]]
null3<-Stab.sel[[1]][[3]][[1]]
null4<-Stab.sel[[1]][[4]][[1]]
datafreq<-tibble(it1=colMeans(null1),it2=colMeans(null2),it3=colMeans(null3),it4=colMeans(null4))
datafreq %>% 
  ggplot(aes(x=it4))+
  geom_point(aes(y=it2),color="orange", size=0.4)+
  geom_point(aes(y=it1),color="deepskyblue", size=0.4)+
  geom_point(aes(y=it3),color="red", size=0.6)+
  geom_abline(linetype="dashed")+
  geom_hline(yintercept = 2/p)+
  geom_vline(xintercept = 2/p)+
  coord_cartesian(xlim=c(0,0.08), ylim=c(0,0.16))




tree<-1*(colMeans( 1*(Stab.sel[[2]]$Pmat>2/p))>freq.sel.thres)
all<-1*(colMeans( 1*(Stab.sel[[3]]$Pmat>2/p))>freq.sel.thres)

table(null,tree)
table(null,all)
table(tree,all)
table(tree,all,null) # 8 en commun aux trois
table(site, date_Site)



summary(Stab.sel[[1]]$iter)
summary(Stab.sel[[2]]$iter)
summary(Stab.sel[[3]]$iter)

deg1<-colSums(F_Vec2Sym(null))
deg2<-colSums(F_Vec2Sym(tree))
deg3<-colSums(F_Vec2Sym(all))

summary(deg3-deg1)

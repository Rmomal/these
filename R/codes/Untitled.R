rm(list=ls()); par(pch=20);
#devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
library(PLNmodels); library(sna);
library(igraph)
library(RColorBrewer)
library(ggplot2)
source('/home/momal/Git/these/pack1/R/codes/FunctionsMatVec.R')
source('/home/momal/Git/these/pack1/R/codes/FunctionsTree.R')
source('/home/momal/Git/these/pack1/R/codes/FunctionsInference.R')
#source('TreeMixture-mac.R')
source('/home/momal/Git/these/pack1/R/codes/TreeMixture-RML.R')
source('/home/momal/Git/these/pack1/R/codes/fonctions.R')
# Data
data.dir = '/home/momal/Git/these/Data/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))
theme_set(theme_bw())
# Parms
Y = as.matrix(Data$count);
O = Data$offset; X = Data$covariates
# Selection des especes
YO = Y/O
Rank = rank(colSums(YO))
#plot(cumsum(sort(colSums(YO))))
Seuil = 20
Ycum = colSums(Y); Order = order(Ycum)
#plot(cumsum(Ycum[Order]), col = 1+(Rank[Order]<Seuil))
Y = Y[, Rank > Seuil];
O = O[, Rank >Seuil];
n = nrow(Y); p = ncol(Y)

#################################################################################
#################################################################################

# Inférences
PLN.vide<-PLN(Y ~ 1 + offset(log(O)))
PLN.tree<-PLN(Y ~ 1 + X$tree+  offset(log(O)))
PLN.treeDist<-PLN(Y ~ 1 +
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
                                   PLN.treeDist$criteria,
                                   PLN.allcovariates$criteria),model=c("vide","tree","treeDist","allcovariates")))

datamodels[,1:5]<-apply(datamodels[,1:5],2,function(x) as.numeric(as.character(x)))
ggplot(datamodels, aes(x = reorder(model, -ICL), y = ICL)) +
  geom_point()+
  geom_point(aes(y=BIC),color="red")+
  geom_point(aes(y=loglik),color="blue")

# datamodels<-datamodels %>% 
#   gather(key, value, -"model")
# ggplot(datamodels,aes(model))+
#   geom_col(aes(y=value,fill=key))

Z.vide = PLN.vide$model_par$Sigma
Z.tree=PLN.tree$model_par$Sigma
Z.treeDist=PLN.treeDist$model_par$Sigma
Z.allcovariates = PLN.allcovariates$model_par$Sigma
T1<-Sys.time()
inf.vide<-TreeGGM(cov2cor(Z.vide),print=FALSE,step="FALSE",maxIter = 150)# 3 minutes
T2<-Sys.time()
difftime(T2,T1)
T1<-Sys.time()
inf.tree<-TreeGGM(cov2cor(Z.tree),print=FALSE,step="FALSE",maxIter = 150)# 3 minutes
T2<-Sys.time()
difftime(T2,T1)
T1<-Sys.time()
inf.treeDist<-TreeGGM(cov2cor(Z.treeDist),print=FALSE,step="FALSE",maxIter = 150)# 3 minutes
T2<-Sys.time()
difftime(T2,T1)
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
X2<-X[,c(6,8)]
Pmat2 = F_ResampleTreePLN(Y2, X2, O2, v=0.5, B=5e2)
Y2<-Y[,1:10]
O2<-O[,1:10]
F_ResampleTreePLN <- function(Y, X, O, v = 0.8, B = 1e2) {
  n = nrow(Y)
  p = ncol(Y)
  P = p * (p - 1) / 2
  V = round(v * n)
  # Pmat = matrix(0, B, P)
  
  res<-mclapply(1:B,function(x){
    cat('\n', b, '')
    sample = sample(1:n, V, replace = F)
    Y.sample = Y[sample,]
    X.sample = X[sample,]
    O.sample = O[sample,]
    
    fmla <- as.formula(paste("Y.sample ~ -1 +offset(log(O.sample))+", 
                             paste(paste0("X.sample$",colnames( X.sample)), collapse= "+")))
    
    PLN.sample = PLN(fmla)
    Sigma.sample = PLN.sample$model_par$Sigma
    inf<-TreeGGM(cov2cor(Sigma.sample), "FALSE", FALSE, maxIter = 150)$P
    return(F_Sym2Vec(inf))
    
  }, mc.cores=2)
  return(res)
}

# Stability selection
Pfreq = 1*(Pmat > 2/p) # Threshold at tree density
Psel = (colMeans(Pfreq) > .5) # Keep edges selected more than half of the time
table(Gvec, Psel)
cat(sum(Gvec), sum(Psel))
Porder = order(Pvec); plot(Pvec[Porder], col=1+Gvec[Porder]); 
abline(h=2/p); abline(h=quantile(Pvec, prob=1-sum(Psel)/P), col=2)



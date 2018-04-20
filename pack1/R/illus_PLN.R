rm(list=ls()); par(pch=20);
#devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
library(PLNmodels); library(sna);
source('/home/momal/Git/these/pack1/R/FunctionsMatVec.R')
source('/home/momal/Git/these/pack1/R/FunctionsTree.R')
source('/home/momal/Git/these/pack1/R/FunctionsInference.R')
source('/home/momal/Git/these/pack1/R/TreeMixture-RML.R')
# Data
data.dir = '/home/momal/Git/these/Data/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))

# Parms
Y = as.matrix(Data$count); n = nrow(Y); p = ncol(Y)
O = Data$offset; X = Data$covariates
dev.off()
plot(X$distTObase^2,X$distTOground^2+X$distTOtrunk^2)
abline(0,1,col="red")
summary(lm(X$distTObase^2 ~ X$distTOground^2+X$distTOtrunk^2))

# PLN models
PLN.offset = PLN(Y ~ 1 + offset(log(O)))
PLN.tree = PLN(Y ~ 1 + X$tree + offset(log(O)))
PLN.tree.base = PLN(Y ~ 1 + X$tree + X$distTObase + offset(log(O)))
PLN.tree.base.infect = PLN(Y ~ 1 + X$tree + X$distTObase + X$pmInfection + offset(log(O)))

# BIC
Crit = rbind(PLN.offset$criteria, PLN.tree$criteria, PLN.tree.base$criteria, PLN.tree.base.infect$criteria)
Crit
apply(Crit, 2, which.max)

# inférences
Z.offset = PLN.offset$model_par$Sigma
inf.offset<-TreeGGM(cov2cor(Z.offset),print=TRUE,step="FALSE")

Z.tree = PLN.tree$model_par$Sigma
inf.offset<-TreeGGM(cov2cor(Z.tree),print=TRUE,step="FALSE")

Z.tree.base = PLN.tree.base$model_par$Sigma
inf.offset<-TreeGGM(cov2cor(Z.tree.base),print=TRUE,step="FALSE")

Z.tree.base.infect = PLN.tree.base.infect$model_par$Sigma
inf.offset<-TreeGGM(cov2cor(Z.tree.base.infect),print=TRUE,step="FALSE")

# Listmodels<-list(PLN.tree,PLN.tree.base,PLN.tree.base.infect)
# ListZSigma<-lapply(Listmodels,function(x) x$model_par$Sigma)
#
# ListInf<-mclapply(ListZSigma,function(x) TreeGGM(cov2cor(x),step="FALSE",print=TRUE),mc.cores=3)

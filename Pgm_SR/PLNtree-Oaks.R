rm(list=ls()); par(pch=20); 
#devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
library(PLNmodels); library(sna); 
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('/home/robin/PgmR/Network/FunctionsTree.R')
source('Functions/FunctionsInference.R')

# Data
data.dir = '/home/momal/Git/these/Data/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))

# Parms
Y = as.matrix(Data$count); n = nrow(Y); p = ncol(Y)
O = Data$offset; X = Data$covariates

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

# PLNnetworks
# PLNnet.tree = PLNnetwork(Y ~ 1 + X$tree + offset(log(O)))
load('/home/robin/Bureau/CountPCA/sparsepca/Pgm/PLNnetwork/output_oaks/network_oak.RData')
models_nocov$plot()
PLNnet_nocov = models_nocov$getBestModel()
G_nocov = 1*(PLNnet_nocov$model_par$Omega!=0)
Pos = gplot(G_nocov, gmode='graph', label=colnames(Y))

models$plot()
PLNnet_cov = models$getBestModel()
G_cov = 1*(PLNnet_cov$model_par$Omega!=0)
gplot(G_cov, gmode='graph', label=colnames(Y), coord=Pos)

plot(colSums(G_nocov), colSums(G_cov))

# # EM beta
# beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)
# P.unif = Kirshner(beta.unif)$P
# alpha = 1
# Cor = cov2cor(PLN.tree$model_par$Sigma)
# log.phi = -.5*n*alpha*log(1-Cor^2); 
# log.phi = log.phi - mean(log.phi[upper.tri(log.phi)]); diag(log.phi) = 0
# phi = exp(log.phi); diag(phi) = 0
# PedgeY.1step = KirshnerRM(beta.unif*phi)$P; sum(PedgeY.1step)
# FitEM = FitBetaStatic(beta.init=beta.unif, phi=phi)
# image(FitEM$beta)

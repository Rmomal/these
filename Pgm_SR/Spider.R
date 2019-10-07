rm(list=ls())
library(mvabund); library(PLNmodels); library(qgraph); library(corpcor); library(EMtree)
library(parallel); 
source('/home/robin/PgmR/General/FunctionsMatVec.R')

# Data
data("spider"); 
Y = as.matrix(spider$abund); X = as.matrix(spider$x)
p = ncol(Y)

# PLN
PLNo = PLN(Y ~ X)
corZ = cov2cor(PLNo$model_par$Sigma)
pCorZ = cor2pcor(corZ)

# PLN network
lambda = exp(seq(2, -4, length.out=50))
PLNnet = PLNnetwork(Y ~ X, penalties=lambda)
PLNnet$plot()
PLNbest = PLNnet$getBestModel()
sparsePCorR = as.matrix(wi2net(PLNbest$model_par$Omega)); diag(sparsePCorR) = 1

# EM tree
EMt = EMtree(PLNo)

# EM tree + resampling
EMres = ResampleEMtree(counts=Y, covariate=X, S=100)
probThresh = 2/p; EMsel = 1*(EMres$Pmat > probThresh)
hist(colMeans(EMsel), breaks=p)
EMfreq = F_Vec2Sym(colMeans(EMsel))
freqThresh = .5; EMnet = 1*(EMfreq > freqThresh)
par(mfrow=c(1, 2), mex=.6); 
# PLNbest$plot_network()
draw_network(EMnet, title='', nodes_label=colnames(Y))

par(mfrow=c(2, 2), mex=.6); 
image(1:p, 1:p, abs(corZ), zlim=c(0, 1)); 
image(1:p, 1:p, abs(pCorZ), zlim=c(0, 1)); 
image(1:p, 1:p, abs(sparsePCorR), zlim=c(0, 1)); 
image(1:p, 1:p, EMfreq, zlim=c(0, 1))


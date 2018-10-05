# PLNtree for Barents fish data

# rm(list=ls()); par(pch=20)
library(PLNmodels); library(sna); library(ade4)
source('../pack1/R/codes/FunctionsMatVec.R')
source('../pack1/R/codes/FunctionsTree.R')
source('../pack1/R/codes/FunctionsInference.R')
source('Functions/FunctionEMtreeResampling.R')
options(error=NULL)

# Data Aravo
data.dir = '../Data_SR/'
data(aravo); data.name = 'Aravo'
Data = list(count=as.matrix(aravo$spe), covariates=aravo$env)
Data$covariates$Aspect = as.factor(Data$covariates$Aspect)
Data$covariates$Form = as.factor(Data$covariates$Form)
Data$covariates$ZoogD = as.factor(Data$covariates$ZoogD)
colnames(Data$count) = aravo$spe.names

# Function
# Algo parms
iterMax = 3e2; 

# Test
n = nrow(Data$count); p = ncol(Data$count); species = colnames(Data$count)
# VEM = PLN(count ~ 1, data=Data); CorY = cov2cor(VEM$model_par$Sigma); save(CorY, file='CorY.Rdata')
# VEM = PLN(count ~ covariates$Aspect + covariates$Slope + covariates$Form + covariates$PhysD + covariates$ZoogD + covariates$Snow, data=Data); CorY = cov2cor(VEM$model_par$Sigma); save(CorY, file='CorY.Rdata')
load('CorY.Rdata')
IsSymm(CorY)
alpha.n = F_AlphaN(CorY, n)
phi = alpha.n$phi
print(alpha.n$alpha)
beta.init = matrix(1, p, p); 
FitEM = FitBetaStatic(Y, beta.init=beta.init, phi=phi, iterMax=iterMax)
plot(FitEM$logpY)

dim(beta.init)

# PLNtree for Barents fish data

rm(list=ls()); par(pch=20)
library(PLNmodels)
source('../pack1/R/codes/FunctionsMatVec.R')
source('../pack1/R/codes/FunctionsTree.R')
source('../pack1/R/codes/FunctionsInference.R')

# Data
data.dir = '../Data_SR/'
data.name = 'BarentsFish_Group'
load(paste0(data.dir, data.name, '.Rdata'))

# Functions
TreeGGM <- function(CorY){
   p = ncol(CorY);
   phi = 1/sqrt(1 - CorY^2); diag(phi) = 0
   beta.unif = matrix(1, p, p); 
   FitEM = FitBetaStatic(Y, beta.init=beta.unif, phi=phi)
   return(list(P=Kirshner(FitEM$beta)$P, L=FitEM$logpY))
}

# Fit VEM-PLN with 1 = no covariates, 2 = depth+temp, 3 = all covariates
# VEM = list()
# VEM[[1]] = PLN(Data$count ~ 1)
# VEM[[2]] = PLN(Data$count ~ Data$covariates[, 1:2])
# VEM[[3]] = PLN(Data$count ~ Data$covariates)
# save(VEM, file=paste0(data.dir, data.name, '-VEM.Rdata'))
load(paste0(data.dir, data.name, '-VEM.Rdata'))
     
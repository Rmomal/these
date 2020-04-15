library(EMtree)
library(PLNmodels)
library(LITree)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(mvtnorm)
library(useful)
library(mclust)
library(MASS)
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(harrypotter)
library(sparsepca)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-exactDet.R")

ExactSumTree<-function(W){
  mat=Laplacian(W)[-1, -1]
  exactDet=det.fractional(mat)
  return(exactDet)
}


NA_SumTree <- readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/NA_SumTree.rds")

SumTree(NA_SumTree$Wg)
ExactSumTree(NA_SumTree$Wg)


SumTree(VEM_1$Wg)
ExactSumTree(VEM_1$Wg)

NA_W<- readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/NA_W.rds")
SumTree(NA_W)
log(ExactSumTree(NA_W))


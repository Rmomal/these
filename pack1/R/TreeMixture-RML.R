# Edge probability in a mixture of trees : GGM
# beta = coef for 'prior' edge probabilities : p(T) = prod_{ij \in T} beta_ij / B
# phi = coef for conditional distribution of : P(Y | T) = prod_{ij \in T} phi_ij

setwd("/home/momal/Git/these/pack1/R")
rm(list=ls()); par(pch=20, mfrow=c(2, 2), mex=3/4);
source('/home/momal/Git/these/pack1/R/FunctionsMatVec.R')
source('/home/momal/Git/these/pack1/R/FunctionsTree.R')
source('/home/momal/Git/these/pack1/R/FunctionsInference.R')
library(sna);
library(readxl)

# Data
Y<-as.matrix(read_excel("~/Documents/codes/Data/Data Files/1. cd3cd28.xls"))[1:200,]

TreeGGM<-function(Y,step){
  p = ncol(Y); n = nrow(Y)
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)

  Cor = cor(Y); alpha = 1; if(n > 50){alpha = 100/n}
  log.phi = -.5*n*alpha*log(1-Cor^2);
  log.phi = log.phi - mean(log.phi[upper.tri(log.phi)]); diag(log.phi) = 0
  phi = exp(log.phi); diag(phi) = 0

  FitEM = switch(step,"FALSE"=FitBetaStatic(beta.init=beta.unif, phi=phi),
                 "TRUE"=FitBeta1step(beta.init=beta.unif, phi=phi))

  return(list(P=Kirshner(FitEM$beta)$P,L=FitEM$logpY))
}


P<-TreeGGM(Y)$P
L<-TreeGGM(Y)$L

## Plots
sum(P)
par(mfrow=c(2,2))
hist(P)
image(P)
plot(L, type='b')
G = matrix(1, p, p); diag(G) = 0;
gplot(G, gmode='graph', edge.lwd=10*P, edge.col="chartreuse4")








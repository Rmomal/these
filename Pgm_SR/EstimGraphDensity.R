source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('/home/robin/PgmR/Network/FunctionsTree.R')
source('Functions/FunctionsInference.R')
source('Functions/FunctionsSimul.R')
source('Functions/FunctionsDivers.R')
library(mvtnorm)

p = 20; n = 50; d = .1
G = erdos(p, d)
G = SimCluster(p, 3, d, 10)
lambda = 1; Omega = diag(rep(lambda, p)) + G; 
Gdiag = G+diag(rep(1, p))
   
while (min(eigen(Omega)$values) < 0){lambda = 1.1*lambda; Omega = diag(rep(lambda, p)) + G}
Sigma = solve(Omega)

Y = rmvnorm(n, sigma=Sigma); 
# Cor = cor(Y); Cov = cov(Y); O = solve(Cov)
Y<-X
seuil<-0.5
# S = 
estimDensity<-function(Y,seuil,plot=TRUE){
  p <- ncol(Y)
  P = matrix(0, p, p)
  for (j in 1:p) {
    LM = lm(Y[, j] ~ -1 + Y[,-j])
    P[j,-j] = summary(LM)$coef[, 4]
    # S[j, -j] = summary(LM)$coef[, 3]
  }
  # par(mfrow=c(2, 2))
  if(plot) hist(as.vector(P), breaks=p)
  # boxplot(as.vector(S) ~ as.vector(Gdiag))
  # qqnorm(as.vector(S[which(Gdiag==0)])); abline(0, 1)
  M = sum(P > seuil)/(1-seuil)
  res<-1-(M)/(p*(p-1))
  print(c(M))
  return(M)
   
}

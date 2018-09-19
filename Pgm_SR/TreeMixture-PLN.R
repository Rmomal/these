# Edge probability in a mixture of trees : PLN model
# beta = coef for 'prior' edge probabilities : p(T) = prod_{ij \in T} beta_ij / B
# alpha = coef for conditional distribution of : P(Y | T) = prod_{ij \in T} alpha_ij

rm(list=ls()); par(pch=20, mfrow=c(2, 2), mex=3/4); 
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('Functions/FunctionsInference.R')
source('Functions/FunctionsSimul.R')
source('Functions/FunctionsTree.R')
source('Functions/FunctionsDivers.R')
library(sna); library(mvtnorm); library(ROCR); library(glasso); library(vegan); library(PLNmodels)

# Parms
p = 10; n = 5e1; B = 5e1
beta = matrix(1, p, p)

# Graph + Sigma
G = SimTree(p)
G = SimErdos(p, 2/p); 
G = SimErdos(p, 4/p)
# while (min(colSums(G)) == 0){G = SimErdos(p, 2/p)}
gplot(G, gmode='graph', label=1:p)
lambda = 1.1; Omega = diag(rep(lambda, p)) + G; 
while (min(eigen(Omega)$values) < 0){lambda = 1.1*lambda; Omega = diag(rep(lambda, p)) + G}
Sigma = solve(Omega)
P0 = EdgeProb(beta)

Perf1 = Perf2 = PerfGL = rep(0, B)
P1 = P2 = matrix(0, B, p*(p-1)/2)
for (b in 1:B){
   cat(b, '')
   Z = rmvnorm(n, sigma=Sigma); 
   Y = matrix(rpois(n*p, exp(Z)), n, p)
   PLN = PLN(Y ~ 1);
   Cov = PLN$model.par$Sigma; Cor = cov2cor(Cov)
   alpha = 1/sqrt(1 - Cor^2); diag(alpha) = 0
   PedgeY1 = EdgeProb(beta*alpha)
   P1[b, ] = F_Sym2Vec(PedgeY1)
   FitEM = FitBetaStatic(beta.init=beta, alpha=alpha)
   # Pb = (min(diff(FitEM$logpY)) < 0)
   beta.hat = FitEM$beta
   PedgeY2 = EdgeProb(beta.hat*alpha)
   P2[b, ] = F_Sym2Vec((PedgeY2))
   # plot(F_Sym2Vec(PedgeY1), F_Sym2Vec(PedgeY2), col=1+F_Sym2Vec(G), xlab='', ylab='')
   Pred1 = prediction(F_Sym2Vec(PedgeY1), F_Sym2Vec(G)); 
   Perf1[b] = performance(Pred1, measure='auc')@y.values[[1]][1]
   Pred2 = prediction(F_Sym2Vec(PedgeY2), F_Sym2Vec(G)); 
   Perf2[b] = performance(Pred2, measure='auc')@y.values[[1]][1]
   PredGL = prediction(F_Sym2Vec(abs(FitGlasso(cov(Y)))), F_Sym2Vec(G)); 
   PerfGL[b] = performance(PredGL, measure='auc')@y.values[[1]][1]
}
boxplot(cbind(Perf1, Perf2, PerfGL), log='y')
boxplot(cbind(as.vector(P1), as.vector(P2)), ylim=c(0, 1))

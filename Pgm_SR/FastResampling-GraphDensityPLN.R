# Fast resampling for PL-network density estimaition
# Based on Schur's trick used in https://arxiv.org/pdf/1505.07281.pdf, p23

rm(list=ls())
source('Functions/FunctionsSimul.R')
source('Functions/FunctionGGMpvalResampling.R')
library(mvtnorm); library(PLNmodels); library(sna)

# Sim parms
p = 30; n = 100; d = .1; P = p*(p-1)

# Model parms
G = SimErdos(p, d); 
D = rep(1, p); lambda = 1; Omega = lambda*D + G; 
while (min(eigen(Omega)$values) < 1e-2){lambda = 2*lambda; Omega = lambda*diag(D) + G}
Sigma = solve(Omega)
O = matrix(rnorm(n*p), n, p)

# Data simulation
Z = rmvnorm(n, sigma=Sigma); 
Y = matrix(rpois(n*p, exp(O+Z)), n, p); 

# PLN inference
PLN = PLN(Y ~ -1 + offset(O)); 
M = PLN$var_par$M; S = PLN$var_par$S; 

# Test statistics + resampled p-values
Stat = F_StatAll(t(M)%*%M + diag(colSums(S)), n, p)
Pval = F_FastResamplePval(M, S, Stat, B=1e3)

# Density estimation
P0 = 2*sum(Pval >= .5, na.rm=T); P1 = P - P0
cat(P0, P1, '/', sum(G), '\n')

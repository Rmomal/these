source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('/home/robin/PgmR/Network/FunctionsTree.R')
source('Functions/FunctionsInference.R')
source('Functions/FunctionsSimul.R')
source('Functions/FunctionsDivers.R')
library(mvtnorm)

p = 20; n = 50; d = .1
G = SimErdos(p, d)
G = SimCluster(p, 3, d, 10)
lambda = 1; Omega = diag(rep(lambda, p)) + G; 
Gdiag = G+diag(rep(1, p))
   
while (min(eigen(Omega)$values) < 0){lambda = 1.1*lambda; Omega = diag(rep(lambda, p)) + G}
Sigma = solve(Omega)

Y = rmvnorm(n, sigma=Sigma); 
Cor = cor(Y); Cov = cov(Y); O = solve(Cov)

S = P = matrix(0, p, p)
for (j in 1:p){
   LM = lm(Y[, j] ~ -1 + Y[, -j])
   P[j, -j] = summary(LM)$coef[, 4]
   S[j, -j] = summary(LM)$coef[, 3]
}
par(mfrow=c(2, 2))
hist(as.vector(P), breaks=p)
boxplot(as.vector(S) ~ as.vector(Gdiag))
qqnorm(as.vector(S[which(Gdiag==0)])); abline(0, 1)
M = 2*sum(P>.5); 
print(c(M, sum(G==0)))

source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('/home/robin/PgmR/Network/FunctionsTree.R')
source('Functions/FunctionsInference.R')
source('Functions/FunctionsSimul.R')
source('Functions/FunctionsDivers.R')
library(mvtnorm)

p = 20; n = 100; d = .1
G = SimErdos(p, d)
G = SimCluster(p, 3, d, 10)
lambda = 1; Omega = diag(rep(lambda, p)) + G; 
Gdiag = G+diag(rep(1, p))
   
while (min(eigen(Omega)$values) < 0){lambda = 1.1*lambda; Omega = diag(rep(lambda, p)) + G}
Sigma = solve(Omega)

Y = rmvnorm(n, sigma=Sigma); 
Cor = cor(Y); Cov = cov(Y); O = solve(Cov)

Stat = Pval = matrix(NA, p, p)
for (j in 1:p){
   LM = lm(Y[, j] ~ -1 + Y[, -j])
   Pval[j, -j] = summary(LM)$coef[, 4]
   Stat[j, -j] = summary(LM)$coef[, 3]
}
par(mfrow=c(2, 2))
hist(as.vector(Pval), breaks=p)
boxplot(as.vector(Stat) ~ as.vector(Gdiag))
qqnorm(as.vector(Stat[which(Gdiag==0)])); abline(0, 1)
M0 = 2*sum(Pval>.5, na.rm=T); 
N0 = sum(G==0)-p
print(c(M0, N0))

P = p*(p-1)
M0 = 2*sum(Pval>.5, na.rm=T); M1 = P - M0
N0 = sum(G==0)-p; N1 = P - N0
print(c(P, N0, N1, N0/P, N1/P))
print(c(P, M0, M1, M0/P, M1/P))


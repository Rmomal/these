# Edge probability in a mixture of trees : GGM
# beta = coef for 'prior' edge probabilities : p(T) = prod_{ij \in T} beta_ij / B
# phi = coef for conditional distribution of : P(Y | T) = prod_{ij \in T} phi_ij

rm(list=ls()); par(pch=20); 
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('Functions/FunctionsInference.R')
source('Functions/FunctionsSimul.R')
source('Functions/FunctionsTree.R')
library(sna); library(mvtnorm); library(gtools)
par(mfrow=c(2, 2))

# Parms
p = 5; n = 1e2; d = 3/p
P = p*(p-1)/2
beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)

# Graph + Sigma
G = SimErdos(p, d); while (min(colSums(G)) == 0){G = SimErdos(p, d)}
# G  = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
gplot(G, gmode='graph', label=1:p)
lambda = 1.1; Omega = diag(rep(lambda, p)) + G; 
while (min(eigen(Omega)$values) < 0){lambda = 1.1*lambda; Omega = diag(rep(lambda, p)) + G}
Sigma = solve(Omega)
P.unif = Kirshner(beta.unif)$P

# Data sampling
Y = rmvnorm(n, sigma=Sigma); Cor = cor(Y); Cov = cov(Y)
phi = 1/sqrt(1 - Cor^2); diag(phi) = 0

if (p <= 4){
   # Grid maximization
   # beta.grid = .5*UnitSimplex(P, 10)$V
   beta.step = .01
   beta.list = .5*qbeta(seq(beta.step, 1-beta.step, by=beta.step), .5, .5)
   beta.grid = list(length=P-1); sapply(1:(P-1), function(pp){beta.grid[[pp]]<<-beta.list})
   beta.grid = expand.grid(beta.grid);
   beta.grid = beta.grid[-which(rowSums(beta.grid)>.5), ]
   beta.grid = as.matrix(cbind(beta.grid, .5-rowSums(beta.grid)))
   colnames(beta.grid) = c(); beta.nb = nrow(beta.grid)
}else{
   # Monte Carlo
   beta.nb = 1e5
   beta.grid = rdirichlet(beta.nb, rep(1/P, P))/2
}

logpY.grid = rep(0, beta.nb)
invisible(sapply(1:beta.nb, function(b){
   beta.tmp = matrix(0, p, p)
   beta.tmp[upper.tri(beta.tmp)] = beta.grid[b, ]
   beta.tmp = beta.tmp + t(beta.tmp)
   logpY.grid[b] <<- log(SumTree(beta.tmp*phi)) - log(SumTree(beta.tmp))
   if(b%%round(beta.nb/log(beta.nb))==0){cat(b, '')}
}))
beta.opt = beta.grid[which.max(logpY.grid), ]

FitEM = FitBetaStatic(beta.init=beta.unif, phi=phi)
plot(FitEM$logpY); abline(h = max(FitEM$logpY)); log(SumTree(FitEM$beta*phi)) - log(SumTree(FitEM$beta))
# beta.tmp = matrix(0, p, p); beta.tmp[upper.tri(beta.tmp)] = beta.opt; beta.tmp = beta.tmp + t(beta.tmp)
# FitEM = FitBetaStatic(beta.init=beta.tmp, phi=phi)
# plot(FitEM$logpY); abline(h = max(FitEM$logpY)); log(SumTree(FitEM$beta*phi))

beta.EM = as.vector(FitEM$beta[upper.tri(FitEM$beta)])
plot(beta.opt, beta.EM); abline(0, 1)
print(rbind(beta.opt, beta.EM))
print(c(logpY.grid[which.max(logpY.grid)], FitEM$logpY[length(FitEM$logpY)]))

# # Prior edge proba
# p.list = 5*(1:20); p.nb = length(p.list); p.slope = rep(0, p.nb)
# for (pp in 1:p.nb){
#    p = p.list[pp];
#    beta = matrix(runif(p^2), p, p)
#    beta[lower.tri(beta, diag=T)] = 0; beta = beta+t(beta); beta = beta/sum(beta)
#    prob = Kirshner(beta)$P
#    plot(beta, prob); abline(0, b=2*p-4, col=2)
#    p.slope[pp] = lm(as.vector(prob) ~ as.vector(beta))$coef[2]
# }
# plot(p.list, p.slope);
# LM = lm(p.slope ~ p.list); print(LM)
# abline(a=LM$coef[1], b=LM$coef[2])
# abline(a=-4, b=2, col=2)

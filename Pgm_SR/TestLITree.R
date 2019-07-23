rm(list=ls())
library(huge); library(saturnin); # library(LITree); 
library(sna); library(MASS); library(mvtnorm); library(Matrix); library(pracma)
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('LITree/LITree/R/functions-estimation.R')

seed = 1; set.seed(seed)
n = 50; d = 20; prob = 3/d

# Simul G
G = F_Vec2Sym(F_Sym2Vec(matrix(rbinom(d^2, 1, prob), d, d)))
while(min(colMeans(G)) == 0){G = F_Vec2Sym(F_Sym2Vec(matrix(rbinom(d^2, 1, prob), d, d)))}
indexMiss = which.max(colSums(G)); k = 1; nodeCol = rep(1, d); nodeCol[indexMiss] = 2; 
gplot(G, gmode='graph', vertex.col=nodeCol)

# Parms
c = 1.1; Omega = diag(c*rep(1, d)) + G
while(min(eigen(Omega)$values) < 1e-6){cat(cond(Omega), ''); c = 1.1*c; Omega = diag(c*rep(1, d)) + G}
Sigma = solve(Omega)
eps = 1; Omega[, indexMiss] = eps*Omega[, indexMiss]; Omega[indexMiss, ] = eps*Omega[indexMiss, ]; Omega[indexMiss, indexMiss] = Omega[indexMiss, indexMiss]/eps; 
SNR = norm(Omega[-indexMiss, indexMiss]%*%solve(Omega[indexMiss, indexMiss])%*%Omega[indexMiss, -indexMiss], '2')^2 / 
   norm(Omega[-indexMiss, -indexMiss], '2')^2

# Simul X
X = rmvnorm(n, sigma=Sigma)

# Satunin
w = lweights_gaussian(X); edgeP = edge.prob(w)
boxplot((edgeP) ~ G)

# Huge
lambda = exp(seq(0, -5, by=-.2)); net = huge(X, lambda=lambda)
score = as.matrix(Reduce('+', net$path)); score = score / max(score)
boxplot(score ~ G)

# Missing edge
# Xmiss = X[, indexMiss]; Xobs = X[, -indexMiss]

# Spatial effect
V = 2*sqrt(max(diag(Sigma)))*rnorm(n); Xobs = X + V%o%rnorm(d); k = 1; Xmiss = V

# Barents
# load('/home/robin/RECHERCHE/ECOLOGIE/CountPCA/sparsepca/Data/BBVA/BarentsFish.Rdata')
# Xobs = Data$count

# Infer missing edge
init <- initEM(Xobs, cliquelist = findCliques(Xobs, k+1)[1:k])
res <- treeAgr.EM(S = cov(Xobs), k=k, K0 = init$K0,
                                Sigma0 = init$Sigma0, pii=0.5, n=nrow(Xobs))
# plot(Sigma[-indexMiss, -indexMiss], res$Sigma[-d, -d]); abline(0, 1)
# plot(Omega[-indexMiss, -indexMiss], res$K[-d, -d]); abline(0, 1)
invSxx = solve(res$Sigma[-d, -d]); Syx = res$Sigma[d, -d]; Syy = res$Sigma[d, d]
# invSxx = solve(Sigma[-indexMiss, -indexMiss]); Syx = Sigma[indexMiss, -indexMiss]; Syy = Sigma[indexMiss, indexMiss]
condEsp = as.vector(Syx%*%invSxx); condVar = (Syy - Syx%*%invSxx%*%Syx)[1, 1]
Xpred = as.vector(condEsp%*%t(Xobs))
plot(Xpred, Xmiss, main=round(cor(Xpred, Xmiss), 2)); abline(0, 1)
# plot(Xpred, Data$covariates[, 4], xlim=c(-100, 100))

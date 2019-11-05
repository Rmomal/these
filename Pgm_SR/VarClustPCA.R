# Variable clustering from hierachical PCA

rm(list=ls()); par(mfrow=c(1, 1))
library(mvtnorm); library(sna)
source('Functions/FunctionsSimul.R')
source('Functions/FunctionsMatVec.R')
source('Functions/FunctionVarClustPCA.R')

# Parms
n = 50; p = 20; seed = 4; set.seed(seed); traceS = TRUE

# # Covariance
# sigma = sqrt(rgamma(p, 2, 1))
# Rho = exp(-as.matrix(dist(matrix(rnorm(2*p), p, 2), diag=TRUE)))
# Sigma = diag(sigma)%*%Rho%*%diag(sigma)
# G = SimErdos(p+1, log(p+1)/(p+1)); while(min(colSums(G))==0){G = SimErdos(p+1, log(p+1)/(p+1))}
# Gcenter = which.max(colSums(G))[1]; centerNeigbors = which(G[Gcenter, -Gcenter]==1)
# vertexCol = rep(2, (p+1)); vertexCol[Gcenter] = 4
# gplot(G, gmode='graph', vertex.col=vertexCol)
# Gsign = F_Vec2Sym(F_Sym2Vec(1 - 2*matrix(rbinom((p+1)^2, 1, .5), (p+1), (p+1))))
# Gsign = G
# lambda = 1.1; Omega = lambda*diag(colSums(G)) - Gsign
# while(min(eigen(Omega)$values) < 1e-6){lambda = 1.1*lambda; Omega = lambda*diag(colSums(G)) - Gsign}
# Sigma = solve(Omega); Sigma = Sigma[-Gcenter, -Gcenter];
#  
# # Data
# Y = rmvnorm(n, sigma=Sigma)
# S = cov(Y); S = cor(Y)

# Barents
library(PLNmodels)
load('../Data_SR/BarentsFish.Rdata')
# PLNall = PLN(Data$count ~ scale(Data$covariates))
# PLNnone = PLN(Data$count ~ 1)
# save(PLNall, PLNnone, file='../Data_SR/BarentsFish-PLN.Rdata')
load('../Data_SR/BarentsFish-PLN.Rdata')
p = ncol(Data$count)
Sall = PLNall$model_par$Sigma * nrow(Data$count)
Snone = PLNnone$model_par$Sigma * nrow(Data$count)
S = cov2cor(Snone)

# Variable clustering
par(mfrow=c(2, 2))
vcPCA = F_VarClustPCA(S)
plot(vcPCA$clustPath$cost, type='b', ylab='', xlab=''); 
plot(sum(diag(S)) - vcPCA$clustPath$height, type='b', ylim=c(0, sum(diag(S))), ylab='', xlab=''); 
abline(h = vcPCA$lastCost, col=2)

# Comparison with ClustOfVar
library(ClustOfVar)
HCV = hclustvar2(S)
plot(HCV)
plot(cumsum(HCV$height), type='b'); points(vcPCA$clustPath$height, col=2, type='b')
cbind(vcPCA$clustMerge, HCV$merge)

# # Old plots
# plot(vcPCA$clustPath$height, vcPCA$clustMat[1, ], ylim=c(0, 2*p), type='s', col=1)
# sapply(2:p, function(j){lines(vcPCA$clustPath$height, vcPCA$clustMat[j, ], type='s', col=j)})
# abline(v=vcPCA$clustPath$height[step], col=2)

# Results
image(t(PLNall$model_par$Theta[HCV$order, -1]), main='beta')
image(PLNnone$model_par$Sigma[HCV$order, HCV$order], main='Sigma none')
image(PLNall$model_par$Sigma[HCV$order, HCV$order], main='Sigma all')


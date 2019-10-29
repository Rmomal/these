# Variable clustering from hierachical PCA

rm(list=ls())
library(mvtnorm); library(sna)
source('Functions/FunctionsSimul.R')
source('Functions/FunctionsMatVec.R')
source('Functions/FunctionVarClustPCA.R')

# Parms
n = 50; p = 20; seed = 4; set.seed(seed); traceS = TRUE
par(mfrow=c(5, 5), mex=.6)

# Covariance
sigma = sqrt(rgamma(p, 2, 1))
Rho = exp(-as.matrix(dist(matrix(rnorm(2*p), p, 2), diag=TRUE)))
Sigma = diag(sigma)%*%Rho%*%diag(sigma)
G = SimErdos(p+1, log(p+1)/(p+1)); while(min(colSums(G))==0){G = SimErdos(p+1, log(p+1)/(p+1))}
Gcenter = which.max(colSums(G))[1]; centerNeigbors = which(G[Gcenter, -Gcenter]==1)
vertexCol = rep(2, (p+1)); vertexCol[Gcenter] = 4
gplot(G, gmode='graph', vertex.col=vertexCol)
Gsign = F_Vec2Sym(F_Sym2Vec(1 - 2*matrix(rbinom((p+1)^2, 1, .5), (p+1), (p+1))))
Gsign = G
lambda = 1.1; Omega = lambda*diag(colSums(G)) - Gsign
while(min(eigen(Omega)$values) < 1e-6){lambda = 1.1*lambda; Omega = lambda*diag(colSums(G)) - Gsign}
Sigma = solve(Omega); Sigma = Sigma[-Gcenter, -Gcenter]; 

# Data
Y = rmvnorm(n, sigma=Sigma)
S = cov(Y); # S = cor(Y)

# Barents
library(PLNmodels)
load('../Data_SR/BarentsFish.Rdata')
p = ncol(Data$count)
PLN = PLN(Data$count ~ Data$covariates)
S = PLN$model_par$Sigma * nrow(Data$count)

# Variable clustering
vcPCA = F_VarClustPCA(S)

# Results
plot(sum(diag(S)) - vcPCA$clustPath$height, type='b', ylim=c(0, sum(diag(S))), ylab='', xlab='')
abline(h = vcPCA$lastCost)

# Clustering matrix
plot(vcPCA$clustPath$height, vcPCA$clustMat[1, ], ylim=c(0, 2*p), type='s', col=1)
sapply(2:p, function(j){lines(vcPCA$clustPath$height, vcPCA$clustMat[j, ], type='s', col=j)})

# # Comparison with ClustOfVar
# library(ClustOfVar)
# plot(hclustvar(Y)$height, type='b'); points(clustPath[, 8], col=2, type='b')

# Sampling from the tree distribution

rm(list=ls()); par(pch=20)
library(sna); library(ape)
source('Functions/FunctionsMatVec.R')
source('Functions/FunctionsTree.R')

# Dims
p <- 10; P <- p*(p-1)/2
B <- 1e4

# Parms
beta <- F_Vec2Sym(F_Sym2Vec(matrix(exp(runif(p^2)), p, p)))
prob <- EdgeProba(beta); diag(prob) <- 0
probVec <- F_Sym2Vec(prob)
print(c(2/p, mean(probVec)))

# Function
rSpanTree <- function(prob){
   # Approximate sampling of a spanning according to edge probabilities
   # Frequency of edge selection empricially OK
   # No guaranty about the actual distribution if the whole tree
   p <- nrow(prob)
   g <- F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, prob), p, p)))
   while(!is.connected(g)){g <- F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, prob), p, p)))}
   w <- F_Vec2Sym(F_Sym2Vec(matrix(runif(p^2), p, p)))
   w <- -log(w)
   w[which(g==0)] <- Inf
   return(mst(w))
}

# Tree sampling
edgeFreq <- rep(0, P); treeNum <- rep(0, B); binCode <- 1.1^(0:(P-1))
for(b in 1:B){
   if(b%%round(sqrt(B))==0){cat(b, '')}
   tree <- rSpanTree(prob)
   treeVec <- F_Sym2Vec(tree)
   treeNum[b] <- treeVec %*% binCode
   edgeFreq <- edgeFreq + treeVec
}
sum(edgeFreq)
edgeFreq <- edgeFreq/B
print(c(2/p, mean(probVec), mean(edgeFreq)))
plot(cumsum(sort(table(treeNum), decreasing=TRUE)))

# Plots
par(mfrow=c(3, 1), mex=.6)
plot(cumsum(sort(table(treeNum), decreasing=TRUE)))
plot(sort(edgeFreq))
plot(probVec, edgeFreq-probVec); 
abline(0, 0)
lines(sort(probVec), qbinom(.025, B, sort(probVec))/B-sort(probVec), lty=2)
lines(sort(probVec), qbinom(.975, B, sort(probVec))/B-sort(probVec), lty=2)

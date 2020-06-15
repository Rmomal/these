# Sampling from the tree distribution

rm(list=ls());
library(sna); library(ape)
source('Functions/FunctionsMatVec.R')
source('Functions/FunctionsTree.R')
par(mfcol=c(3, 2), mex=.6, pch=20)

# Dims
p <- 15; P <- p*(p-1)/2; B <- 1e3

# Parms
unifCst <- SumTree(matrix(1, p, p))
beta <- F_Vec2Sym(F_Sym2Vec(matrix(rgamma(p^2, 1, 2), p, p)))
beta <- beta / SumTree(beta)^(1/(p-1))
betaVec <- F_Sym2Vec(beta)
prob <- EdgeProba(beta); diag(prob) <- 0
probVec <- F_Sym2Vec(prob)
print(c(2/p, mean(probVec)))
hist(betaVec, breaks=p)
hist(probVec, breaks=p)
plot(betaVec, probVec); abline(0, 1); abline(0, sqrt(p))

# Function
rSpanTreeV0 <- function(prob){
   # Approximate sampling of a spanning according to edge probabilities
   # Frequency of edge selection empricially OK
   # No guaranty about the actual distribution if the whole tree
   p <- nrow(prob)
   graph <- F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, prob), p, p)))
   while(!is.connected(graph)){graph <- F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, prob), p, p)))}
   w <- F_Vec2Sym(F_Sym2Vec(matrix(runif(p^2), p, p)))
   w <- -log(w)
   w[which(graph==0)] <- Inf
   return(mst(w))
}

# Function
rSpanTreeV1 <- function(beta, prob){
   # Approximate sampling of a spanning according to edge probabilities
   p <- nrow(prob); OK <- FALSE; tries <- 0; pb <- 0
   while(!OK){
      tries <- tries + 1
      
      graph <- F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, prob), p, p)))
      while(!is.connected(graph)){graph <- F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, prob), p, p)))}
      
      w <- F_Vec2Sym(F_Sym2Vec(matrix(runif(p^2), p, p)))
      w <- -log(w)
      w[which(graph==0)] <- Inf
      tree <- mst(w) 
      
      qTree <- p^2 *prod(F_Sym2Vec(prob)[which(F_Sym2Vec(tree)==1)]) / SumTree(graph)
      pTree <- prod(F_Sym2Vec(beta)[which(F_Sym2Vec(tree)==1)])
      if(pTree > qTree){pb <- pb + 1}
      
      if(runif(1) < pTree / qTree){OK <- TRUE}
   }
      
   return(list(tree=tree, qTree=qTree, pTree=pTree, tries=tries, pb=pb))
}

# Tree sampling
edgeFreq <- rep(0, P); binCode <- 1.1^(0:(P-1))
tries <- qTree <- pTree <- treeNum <- rep(0, B); 
for(b in 1:B){
   if(b%%round(sqrt(B))==0){cat(b, '')}
   tree <- rSpanTreeV1(beta, prob)
   treeVec <- F_Sym2Vec(tree$tree)
   edgeFreq <- edgeFreq + treeVec
   treeNum[b] <- treeVec %*% binCode
   qTree[b] <- tree$qTree
   pTree[b] <- tree$pTree
   tries[b] <- tree$tries
}
sum(edgeFreq)
edgeFreq <- edgeFreq/B
print(c(2/p, mean(probVec), mean(edgeFreq)))
summary(qTree/pTree)
summary(tries); summary(pb/tries)


# Plots
plot(cumsum(sort(table(treeNum), decreasing=TRUE)))
plot(sort(edgeFreq))
plot(probVec, edgeFreq-probVec); 
abline(0, 0)
lines(sort(probVec), qbinom(.025, B, sort(probVec))/B-sort(probVec), lty=2)
lines(sort(probVec), qbinom(.975, B, sort(probVec))/B-sort(probVec), lty=2)

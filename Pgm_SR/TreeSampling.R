# Sampling from the tree distribution

rm(list=ls());
library(sna); library(ape)
source('Functions/FunctionsMatVec.R')
source('Functions/FunctionsTree.R')
par(mfcol=c(3, 2), mex=.6, pch=20)

# Dims
p <- 10; P <- p*(p-1)/2; B <- 1e2
seed <- 1; set.seed(seed)

# Parms
unifCst <- SumTree(matrix(1, p, p))
beta <- F_Vec2Sym(F_Sym2Vec(matrix(rgamma(p^2, 1, 2), p, p)))
beta <- beta / SumTree(beta)^(1/(p-1)); betaVec <- F_Sym2Vec(beta)

# Edge prob
prob <- EdgeProba(beta); diag(prob) <- 0; probVec <- F_Sym2Vec(prob)
print(c(2/p, mean(probVec)))
hist(betaVec, breaks=p); hist(probVec, breaks=p)
plot(betaVec, probVec); abline(0, 1); abline(0, sqrt(p)); abline(0, max(prob/beta, na.rm=TRUE))

# Function
# rSpanTreeV0 <- function(prob){
#    # Approximate sampling of a spanning according to edge probabilities
#    # Frequency of edge selection empricially OK
#    # No guaranty about the actual distribution if the whole tree
#    p <- nrow(prob)
#    graph <- F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, prob), p, p)))
#    while(!is.connected(graph)){graph <- F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, prob), p, p)))}
#    w <- F_Vec2Sym(F_Sym2Vec(matrix(runif(p^2), p, p)))
#    w <- -log(w)
#    w[which(graph==0)] <- Inf
#    return(mst(w))
# }

rSpanTreeV1 <- function(beta, prob){
   # Approximate sampling of a spanning according to edge probabilities
   # Rejection sampling approach
   # !!! Beta's must be normalized !!! (constant B = 1)
   # Max ratio between proposal and target heuristicaly set to p^2
   p <- nrow(prob); P <- p*(p-1)/2
   prob <- prob / max(prob) # To enforce connectivity
   probVec <- F_Sym2Vec(prob); betaVec <- F_Sym2Vec(beta)
   # Optimal constant
   w <- log(prob/beta); optTreeVec <- F_Sym2Vec(mst(w)); 
   M <- p^(p-2) / exp(sum(F_Sym2Vec(w)[which(optTreeVec==1)]))
   M <- p^2
   OK <- FALSE; tries <- 0; pTree <- qTree <- rep(0, 1e4)
   while(!OK){
      tries <- tries + 1
      
      # graph = heterogeneous ER connected graph
      graph <- F_Vec2Sym(matrix(rbinom(P, 1, probVec)))
      while(!is.connected(graph)){graph <- F_Vec2Sym(matrix(rbinom(P, 1, probVec)))}
      # gplot(graph, gmode='graph', main=round(SumTree(graph)))
      graphVec <- F_Sym2Vec(graph)
      # qGraph <- prod(dbinom(graphVec, 1, probVec))
      
      # tree = uniform spanning tree
      w <- F_Vec2Sym(F_Sym2Vec(matrix(runif(p^2), p, p)))
      w[which(graph==0)] <- Inf
      tree <- mst(w)
      treeVec <- F_Sym2Vec(tree)

      # Monte-Carlo sample of graphs containing the tree
      biasedProbVec <- probVec; biasedProbVec[which(treeVec==1)] <- 1
      invSpanTreeMean <- 0; G <- 1e3
      for(g in 1:G){invSpanTreeMean <- invSpanTreeMean + 1/SumTree(F_Vec2Sym(matrix(rbinom(P, 1, biasedProbVec))))}
      invSpanTreeMean <- invSpanTreeMean/G
      qTree[tries] <- prod(probVec[which(treeVec==1)]) * invSpanTreeMean
      
      # qTree[tries] <- qGraph / SumTree(graph)
      # Actually qTree should be qGraph / SumTree(graph)
      # qTree[tries] <- prod(probVec[which(treeVec==1)]) / SumTree(graph)
      pTree[tries] <- prod(betaVec[which(treeVec==1)])

      if(pTree[tries] > (M * qTree[tries])){cat('[*] ')}
      if(runif(1) < pTree[tries] / (M * qTree[tries])){OK <- TRUE}
   }
   pTree <- pTree[1:tries]; qTree <- qTree[1:tries]
   cat(tries, '')
   
   return(list(tree=tree, qTree=qTree, pTree=pTree, qTree=qTree))
}

# Tree sampling
edgeFreq <- rep(0, P); binCode <- 1.1^(0:(P-1))
pb <- tries <- maxRatio <- minRatio <- treeNum <- rep(0, B); 
for(b in 1:B){
   tree <- rSpanTreeV1(beta, prob)
   treeVec <- F_Sym2Vec(tree$tree)
   edgeFreq <- edgeFreq + treeVec
   treeNum[b] <- treeVec %*% binCode
   tries[b] <- length(tree$pTree)
   pb[b] <- sum(tree$pTree/tree$qTree > p^2)
   minRatio[b] <- min(tree$pTree/tree$qTree)
   maxRatio[b] <- max(tree$pTree/tree$qTree)
   if(b%%round(sqrt(B))==0){
      cat('\n [', b, '] ', sep='')
      print(rbind(summary(tries[1:b]), summary(minRatio[1:b]), summary(maxRatio[1:b]), summary(pb[1:b])))
   }
}
sum(edgeFreq)
edgeFreq <- edgeFreq/B
print(c(2/p, mean(probVec), mean(edgeFreq)))


# Plots
plot(cumsum(sort(table(treeNum), decreasing=TRUE)))
plot(sort(edgeFreq))
plot(probVec, edgeFreq-probVec); 
abline(0, 0)
lines(sort(probVec), qbinom(.025, B, sort(probVec))/B-sort(probVec), lty=2)
lines(sort(probVec), qbinom(.975, B, sort(probVec))/B-sort(probVec), lty=2)

# Variable clustering from hierachical PCA

###############################################################################
WhichMinMat <- function(A){
   jmin = apply(A, 1, which.min)
   imin = which.min(sapply(1:nrow(A), function(i){A[i, jmin[i]]}))
   jmin = jmin[imin]
   return(c(imin, jmin))
}

###############################################################################
MinMat <- function(A){ijmin = WhichMinMat(A); return(A[ijmin[1], ijmin[2]])}

###############################################################################
F_VarClustPCA <- function(S,traceS=FALSE){
   p = ncol(S)
   
   # Cost matrix
   C = matrix(Inf, p, p); PC = array(dim=c(n, p, p))
   sapply(1:(p-1), function(j){sapply((j+1):p, function(k){
      C[j, k] <<- eigen(S[c(j, k), c(j, k)])$values[2]
   })})
   # image(1:p, 1:p, C)
   
   # Hierarchical clustering
   if(traceS){listS = list()}
   Stmp = S; Ctmp = C; 
   clustPath = matrix(0, (p-1), 8); colnames(clustPath) = c('j', 'k', 'j.num', 'k.num', 'clust.num', 'coef.j', 'coef.k', 'cost')
   clustContent = as.list(1:p)
   varNum = 1:p; step = 0; 
   while(step < p-2){
      step = step + 1; 
      # Pair to merge
      jkmin = sort(WhichMinMat(Ctmp)); 
      clustPath[step, 1:2] = varNum[jkmin]; 
      clustPath[step, 3:4] = jkmin; 
      clustPath[step, 5] = p+step; 
      clustContent[[p+step]] = sort(c(clustContent[[varNum[jkmin][1]]], clustContent[[varNum[jkmin][2]]]))
      # Update covariance
      eigenSjk = eigen(Stmp[jkmin, jkmin]); 
      clustPath[step, 6:7] = eigenSjk$vectors[, 1]
      Spp = Stmp[, jkmin]%*%eigenSjk$vectors[, 1]
      Stmp = rbind(cbind(Stmp, Spp), cbind(t(Spp), eigenSjk$values[1])); 
      # Update distances
      Cpp = sapply(1:(ncol(Stmp)-1), function(j){
         eigen(Stmp[c(j, ncol(Stmp)), c(j, ncol(Stmp))])$values[2]
      }); 
      Ctmp = rbind(cbind(Ctmp, Cpp), rep(Inf, (ncol(Ctmp)+1))); 
      clustPath[step, 8] = eigenSjk$values[2]
      # Removing merges rows and columns
      Stmp = Stmp[-jkmin[2], ]; Stmp = Stmp[, -jkmin[2]]; 
      Stmp = Stmp[-jkmin[1], ]; Stmp = Stmp[, -jkmin[1]]; 
      if(traceS){listS[[step]]=Stmp}
      Ctmp = as.matrix(Ctmp)[-jkmin[2], ]; Ctmp = as.matrix(Ctmp)[, -jkmin[2]]; 
      Ctmp = as.matrix(Ctmp)[-jkmin[1], ]; Ctmp = as.matrix(Ctmp)[, -jkmin[1]]; 
      varNum = varNum[-jkmin]; varNum[p-step] = p+step
      cat(step, ':', clustPath[step, ], '\n')
      # image(Ctmp, main=paste(step, dim(Ctmp)[1]), xlab='', ylab='')
   }
   # Last step
   eigenSjk = eigen(Stmp); 
   clustPath[p-1, ] = c(varNum, c(1, 2), 2*p-1, eigenSjk$vectors[, 1], eigenSjk$values[2])
   clustContent[[2*p-1]] = sort(c(clustContent[[varNum[jkmin][1]]], clustContent[[varNum[jkmin][2]]]))
   clustPath = as.data.frame(clustPath); clustPath$height = cumsum(clustPath$cost)
   
   # Clustering matrix
   clustMatrix = (1:p) %o% rep(1, (p-1))
   sapply(1:(p-1), function(h){
      if(h>1){clustMatrix[, h] <<- clustMatrix[, (h-1)]}
      clustMatrix[clustContent[[p+h]], h] <<- p+h
   })
   
   # Merging path 
   clustMerge = clustPath[, 1:2]
   sapply(1:2, function(c){
      clustMerge[which(clustMerge[, c] <= p), c] <<- -clustMerge[which(clustMerge[, c] <= p), c]
      clustMerge[which(clustMerge[, c] > p), c] <<- clustMerge[which(clustMerge[, c] > p), c] - p
   })
   
   
   return(list(clustPath=clustPath, clustContent=clustContent, clustMatrix=clustMatrix, 
               lastCost=eigenSjk$values[1], clustMerge=clustMerge))
}

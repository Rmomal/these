######################################################################################
# Simul
SimErdos <- function(p, prob){
   G = matrix(rbinom(p^2, 1, prob), p, p)
   G[lower.tri(G, diag=T)] = 0; G = G+t(G)
}
SimTree <- function(p){
   W = matrix(runif(p^2), p, p); W = W + t(W)
   Tree = spantree(W)
   G = matrix(0, p, p)
   invisible(sapply(1:length(Tree$kid), function(i){G[i+1, Tree$kid[i]] <<- 1}))
   G = G + t(G)
   return(G)
}
SimCluster <- function(p, k, d, r){
   # k = nb clusters, d = graph density, r = within/between connection probability
   beta = d/(r/k + (k-1)/k); alpha = r*beta
   while(alpha > 1){r = .9*r; beta = d/(r/k + (k-1)/k); alpha = r*beta}
   Z = t(rmultinom(p, 1, rep(1/k, k)))
   ZZ = Z %*% t(Z); diag(ZZ) = 0; ZZ = F_Sym2Vec(ZZ)
   G = F_Vec2Sym(rbinom(p*(p-1)/2, 1, alpha*ZZ + beta*(1-ZZ)))
   # gplot(G, gmode='graph', label=1:p, vertex.col=Z%*%(1:k))
   return(G)
}



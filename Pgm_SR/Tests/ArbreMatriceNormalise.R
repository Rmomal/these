# Normalisation des poids -> thm arbre-matrice

rm(list=ls())
library(vegan); library(sna);

# Parms
p = 5
Bprob = 1e2; Buv = Bprob
par(pch=20)

# Fonctions
MakeSym <- function(A){
   A[lower.tri(A)] = 0; A = A + t(A); diag(A) = 0
   return(A)
}
SampTree <- function(p){
   U =  MakeSym(matrix(runif(p^2), p, p))
   Span = spantree(U)
   return(Span$kid)
}
AdjTree <- function(K){
   p = length(K)+1; G = matrix(0, p, p); 
   sapply(2:p, function(i){G[i, K[i-1]] <<- 1})
   sapply(2:p, function(j){G[K[j-1], j] <<- 1})
   return(G)
}
ProdTree <- function(K, W){
   # Calcul du produit en prenant le 1er noeud comme racine 
   # => utiliser le thm arbre-matrice avec L^{1,v}, par exemple L^{1,1}
   P = 1; p = length(K)+1
   # Codage enfant : notation de spantree$kid : K[i] = enfant du noeud i, i > 1
   sapply(2:p, function(i){P <<- P * W[i, K[i-1]]})
   # # Codage parent : notation inverse de spantree$kid : K[j] = parent du noeud j, j > 1
   # sapply(2:p, function(j){P <<- P * W[K[j-1], j]})
   return(P)
}

# Poids originaux
# A = MakeSym(matrix(1, p, p))
A = MakeSym(matrix(exp(rnorm(p^2)), p, p))
L = diag(rowSums(A)) - A; 
# sapply(1:p, function(i){sapply(1:p, function(j){(-1)^(i+j)*det(L[-i, -j])})})

# Normalisation de A
d = 1/colSums(A)
# d = exp(rnorm(p))
AA = diag(d) %*% A # normalisation en ligne
# AA = A %*% diag(d) # normalisation en colonne
LL = diag(rowSums(AA)) - AA; 
# sapply(1:p, function(i){sapply(1:p, function(j){(-1)^(i+j)*det(LL[-i, -j])})})

# Test
K = SampTree(p); G = AdjTree(K); print(G); gplot(G, gmode='graph', label=(1:p))
cat(ProdTree(K, A), ProdTree(K, AA), ProdTree(K, AA)*d[1]/prod(d), '\n'); 

# Effet de (u, v)
UV = matrix(0, Buv, 2); ProbRatio = TestRatio = rep(0, Buv)
for (b in 1:Buv){
   uv = sample(1:p, 2, replace=T); u = uv[1]; v = uv[2]
   UV[b, ] = uv
   TestRatio[b] = d[u] / d[1] # normalisation en ligne
   # TestRatio[b] = d[v] / d[1] # normalisation en colonne
   Atot = (-1)^sum(uv)*det(L[-u, -v])
   AAtot = (-1)^sum(uv)*det(LL[-u, -v])
   cat(Atot, AAtot, AAtot*d[v]/prod(d), '\n')
   Prob = matrix(0, Bprob, 2)
   for (bb in 1:Bprob){
      G = SampTree(p)
      Prob[bb, ] = c(ProdTree(G, A)/Atot, ProdTree(G, AA)/AAtot)
      }
   ProbRatio[b] = mean(Prob[, 2]/Prob[, 1])
}
# cbind(UV, ProbRatio)
table(UV[which(abs(ProbRatio-1)<1e-10), 1])
table(UV[which(abs(ProbRatio-1)<1e-10), 2])
plot(TestRatio, ProbRatio, log='xy', col=(UV[, 1])); abline(a=0, b=1, col=2)

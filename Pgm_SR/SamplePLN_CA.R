# Sampling from tre-PLN

rm(list=ls())
library(ape)
library(sna)
library(poilog)

# Parms
p = 10

# Function
F_PLN <- function(Y, Tree, mu, sigma, Rho){
  p = length(Y)
  phi = sapply(1:p, function(j){dpoilog(Y[j], mu=mu[j], sig=sigma[j])})
  PLN = prod(phi)
  sapply(1:(p-1), function(j){sapply((j+1):p, function(k){
    if(Tree[j, k]==1){
      PLN <<- PLN * dbipoilog(Y[j], Y[k], mu1=mu[j], mu2=mu[k], 
                              sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k]) / 
        phi[j] / phi[k]
    }
  })})
  return(PLN)
}

F_Ratio <- function(Y, Tree, mu, sigma, Rho, Lambda){
  Pois.Y = prod(dpois(Y, Lambda))
  PLN.Y = F_PLN(Y, Tree, mu, sigma, Rho)
  return(PLN.Y/Pois.Y)
}

# Tree sampling
W = matrix(runif(p^2), p, p); W = W + t(W)
Tree = mst(W)
gplot(Tree, gmode='graph')

# Parms
mu = rep(1, p)
coef = 1.1; Omega = Tree + coef*diag(p)
while(min(eigen(Omega)$values)<0){coef = 1.1*coef; Omega = Tree + coef*diag(p)}
Sigma = solve(Omega); Rho = cov2cor(Sigma); sigma = sqrt(diag(Sigma))
Lambda = exp(mu )

# Simul
B <- 1e4 
s<-0
ratio =  rep(0, B); 
Y = matrix(0, B, p);
for (b in 1:B){
  if (b %% round(sqrt(B))==0){cat(b, '')}
  Ytmp = rpois(p, Lambda)
  ratio[b] = F_Ratio(Ytmp, Tree, mu, sigma, Rho, Lambda)
  if (runif(1) < ratio[b]){s = s+1; Y[s, ] = Ytmp}
}
s/B
summary(log(ratio))
hist(log(ratio))
Y = Y[1:s, ]

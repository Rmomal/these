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
F_LogRatio <- function(Y, Tree, mu, sigma, Rho, Lambda){
   logPois.Y = sum(dpois(Y, Lambda, log=T))
   logPLN.Y = log(F_PLN(Y, Tree, mu, sigma, Rho))
   return(logPLN.Y - logPois.Y)
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
Lambda = exp(mu + diag(Sigma)/2)

# Simul
B = 1.2e4; 
Y.cur = rpois(p, Lambda); logRatio.cur = F_LogRatio(Y.cur, Tree, mu, sigma, Rho, Lambda)
alpha =  rep(0, B); Y = matrix(0, B, p); 
Y[1, ] = Y.cur; accept = 1
for (b in 2:B){
   if (b %% round(sqrt(B))==0){cat(b, '')}
   Y.prop = rpois(p, Lambda)
   logRatio.prop = F_LogRatio(Y.prop, Tree, mu, sigma, Rho, Lambda)
   alpha[b] = exp(logRatio.prop - logRatio.cur)
   if (runif(1) < alpha[b]){
      accept = accept+1; 
      Y[b, ] = Y.cur = Y.prop
      logRatio.cur = logRatio.prop
      }else{Y[b, ] = Y.cur}
}
accept/B
summary(log(alpha))
hist(log(alpha))
Y = Y[, -(1:round(.2*B))] # remove burn-in (needed ?)
Y = Y[100*(1:round(B/100)), ]

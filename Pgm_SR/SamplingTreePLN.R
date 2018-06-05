# Sampling from tre-PLN

rm(list=ls()); par(mfrow=c(1, 1))
library(ape); library(sna); library(poilog)

# Parms
p = 10; B.path = 1e5; R.subsample = 5e2

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
F_LogRatioPois <- function(Y, Tree, mu, sigma, Rho, Esp.Y){
   logPois.Y = sum(dpois(Y, Esp.Y, log=T))
   logPLN.Y = log(F_PLN(Y, Tree, mu, sigma, Rho))
   return(logPLN.Y - logPois.Y)
}
F_LogRatioNegBin <- function(Y, Tree, mu, sigma, Rho, size.Y, Esp.Y){
   logNegBin.Y = sum(dnbinom(Y, size=size.Y, mu=Esp.Y, log=T))
   logPLN.Y = log(F_PLN(Y, Tree, mu, sigma, Rho))
   return(logPLN.Y - logNegBin.Y)
}

# Tree sampling
W = matrix(runif(p^2), p, p); W = W + t(W)
Tree = mst(W); gplot(Tree, gmode='graph', label=1:p)
degree = colSums(Tree)

# Parms
mu = rep(1, p); coef = 1.1; Omega = Tree + coef*diag(degree)
while(min(eigen(Omega)$values)<0){coef = 1.1*coef; Omega = Tree + coef*diag(p)}
Sigma = solve(Omega); Rho = cov2cor(Sigma); sigma = sqrt(diag(Sigma))
Esp.Y = exp(mu + sigma^2/2)
Var.Y = Esp.Y + Esp.Y^2*(exp(sigma^2)-1)
Cov.Y = (Esp.Y %o% Esp.Y) * (exp(Sigma)-1)
size.Y = Esp.Y^2 / (Var.Y - Esp.Y)

# # Proposal check
# Y.NegBin = matrix(0, B.path, p)
# invisible(sapply(1:B.path, function(b){Y.NegBin[b, ] <<- rnbinom(p, size=size.Y, mu=Esp.Y)}))
# print(cbind(Esp.Y, colMeans(Y.NegBin)))
# print(cbind(Var.Y, apply(Y.NegBin, 2, var)))

# Simul
B.burn = round(.2*B.path); B = B.burn+B.path
# Y.cur = rpois(p, Esp.Y); logRatio.cur = F_LogRatioPois(Y.cur, Tree, mu, sigma, Rho, Esp.Y)
Y.cur = rnbinom(p, size=size.Y, mu=Esp.Y); 
logRatio.cur = F_LogRatioNegBin(Y.cur, Tree, mu, sigma, Rho, size.Y, Esp.Y)
alpha =  rep(0, B); Y.path = matrix(0, B, p); 
Y.path[1, ] = Y.cur; accept = 1
for (b in 2:B){
   if (b %% round(sqrt(B))==0){
      cat(b, '')
      par(mfrow=c(5, 2), mex=.3); sapply(1:p, function(j){plot(log(1+Y.path[1:b, j]), type='l', main='', xlab='', ylab='')})
      }
   # Y.prop = rpois(p, Esp.Y)
   # logRatio.prop = F_LogRatioPois(Y.prop, Tree, mu, sigma, Rho, Esp.Y)
   Y.prop = rnbinom(p, size=size.Y, mu=Esp.Y)
   logRatio.prop = F_LogRatioNegBin(Y.prop, Tree, mu, sigma, Rho, size.Y, Esp.Y)
   alpha[b] = exp(logRatio.prop - logRatio.cur)
   if (runif(1) < alpha[b]){
      accept = accept+1; 
      Y.path[b, ] = Y.cur = Y.prop
      logRatio.cur = logRatio.prop
      }else{Y.path[b, ] = Y.cur}
}
Y.sample = Y.path[-(1:B.burn), ] # remove burn-in (needed ?)
Y.sample = Y.sample[R.subsample*(1:floor(B.path/R.subsample)), ] # sample 1 out of 100
B.sample = nrow(Y.sample)

# Check
par(mfrow=c(5, 2), mex=.3); 
sapply(1:p, function(j){plot(log(1+Y.path[, j]), type='l', main='', xlab='', ylab='')})
par(mfrow=c(4, 3), mex=.5)
for(j in 1:p){plot(log(1+Y.sample[1:(B.sample-1), j]), log(1+Y.sample[2:B.sample, j]), pch=20)}
print(rbind(Esp.Y, colMeans(Y.sample)))
print(rbind(Var.Y, apply(Y.sample, 2, var)))
par(mfrow=c(2, 2))
plot(Esp.Y, colMeans(Y.sample), pch=20, log='xy'); abline(0, 1)
plot(Var.Y, apply(Y.sample, 2, var), pch=20, log='xy'); abline(0, 1)
plot(Cov.Y, cov(Y.sample), pch=20); abline(0, 1)
plot(cov2cor(Cov.Y), cor(Y.sample), pch=20); abline(0, 1)
accept/B
max(diff(which(Y.path[2:B, 1]!=Y.path[1:(B-1), 1])))

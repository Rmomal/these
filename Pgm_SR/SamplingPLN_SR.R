# Sampling from tre-PLN

rm(list = ls())
library(ape)
library(sna)
library(poilog)

# Parms
p = 10

# Function
F_PLN <- function(Y, Tree, mu, sigma, Rho) {
  p = length(Y)
  uni = sapply(1:p, function(j) {
    dpoilog(Y[j], mu = mu[j], sig = sigma[j])
  })
  PLN = prod(uni)
  sapply(1:(p - 1), function(j) {
    sapply((j + 1):p, function(k) {
      if (Tree[j, k] == 1) {
        PLN <<- PLN * dbipoilog(
          Y[j],
          Y[k],
          mu1 = mu[j],
          mu2 = mu[k],
          sig1 = sigma[j],
          sig2 = sigma[k],
          rho = Rho[j, k]
        ) /(uni[j] *uni[k])
      }
    })
  })
  return(PLN)
}

F_Ratio <- function(Ytmp, Tree, mu, sigma, Rho, Lambda) {
  Pois.Y = prod(dpois(Ytmp, Lambda))
  PLN.Y = F_PLN(Ytmp, Tree, mu, sigma, Rho)
  return(PLN.Y / Pois.Y)
}


# Tree sampling
W = matrix(runif(p ^ 2), p, p)
W = W + t(W)
Tree = mst(W)
gplot(Tree, gmode = 'graph',
      vertex.col="deepskyblue",
      vertex.border="white")

# Parms
mu = rep(1, p)
coef = 1.1
Omega = Tree + coef * diag(p)
while (min(eigen(Omega)$values) < 0) {
  coef = 1.1 * coef
  Omega = Tree + coef * diag(p)
}
Sigma = solve(Omega)
Rho = cov2cor(Sigma)
sigma = sqrt(diag(Sigma))
Lambda = exp(mu + diag(Sigma) / 2)

# Simul
B = 1e4

Ytmp = rpois(p, Lambda)
while (F_Ratio(Ytmp, Tree, mu, sigma, Rho, Lambda) < runif(1)  ) {
  Ytmp = rpois(p, Lambda)
}

ratio =  rep(0, B)
Y = matrix(0, B, p)
Y[1,] = Ytmp
s = 1

for (b in 2:B) {
  if (b %% round(sqrt(B)) == 0) {
    cat(b, '')
  }
  Ytmp = rpois(p, Lambda)
  ratio[b] = F_Ratio(Ytmp, Tree, mu, sigma, Rho, Lambda)
  if (runif(1) < ratio[b]) {
    s = s + 1
    Y[b,] = Ytmp
  } else{
    Y[b,] = Y[b - 1,]
  }
}

s / B
summary(log(ratio))
hist(log(ratio))
Y = Y[100 * (1:round(B / 100)),]





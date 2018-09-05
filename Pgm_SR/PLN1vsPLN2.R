# Joint distribution of PLN1 and PLN2 models
# Graphical model 2 -- 1 -- 3

rm(list=ls())
par(mfrow=c(1, 1), pch=20);
library(mvtnorm); library(sna); library(utils); library(poilog)

# Parms
B = 1e6
Omega = matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3, 3)
p = nrow(Omega)
gplot(1*(Omega!=0), gmode='graph', label=(1:p))
Sigma = solve(Omega)
rho = cov2cor(Sigma)
sigma = sqrt(diag(Sigma))
y.max = 4

##################################################################################
p = 1
##################################################################################
# Y table
Y.tab = expand.grid(0:y.max); nb.tab = nrow(Y.tab)

# PLN1
Z.PLN1 = rnorm(B, sd=sigma[1])
Y.PLN1 = rpois(B, exp(Z.PLN1))
PLN1.tab = sapply(1:nb.tab, function(i){sum(Y.PLN1==Y.tab[i, 1])})/B

# PLN2
PLN2.tab = sapply(1:nb.tab, function(i){dpoilog(Y.tab[i, 1], mu=0, sig=sigma[1])})

cbind(Y.tab, PLN1.tab, PLN2.tab)
plot(PLN2.tab, PLN1.tab-PLN2.tab); abline(h=0)
points(PLN2.tab, qbinom(.025, B, PLN2.tab)/B-PLN2.tab, pch='+', col=4)
points(PLN2.tab, qbinom(.975, B, PLN2.tab)/B-PLN2.tab, pch='+', col=4)

invisible(sapply(1:p, function(j){
   cat('Y', j, '\n')
   P1 = as.vector(by(PLN1.tab, Y.tab[, j], sum))
   P2 = as.vector(by(PLN2.tab, Y.tab[, j], sum))
   L2 = qbinom(.025, B, P2)/B
   U2 = qbinom(.975, B, P2)/B
   print(cbind(sort(unique(Y.tab[, j])), P2, L2, P1, U2))
}))

##################################################################################
p = 2
##################################################################################
# Y table
Y.tab = expand.grid(0:y.max, 0:y.max); nb.tab = nrow(Y.tab)

# PLN1
Z.PLN1 = rmvnorm(B, sigma=Sigma[1:2, 1:2])
Y.PLN1 = matrix(rpois(B*p, exp(Z.PLN1)), B, p)
PLN1.tab = sapply(1:nb.tab, function(i){sum((Y.PLN1[, 1]==Y.tab[i, 1]) &
                                               (Y.PLN1[, 2]==Y.tab[i, 2]))})/B

# PLN2
PLN2.tab = sapply(1:nb.tab, function(i){
   dbipoilog(Y.tab[i, 1], Y.tab[i, 2], mu1=0, mu2=0, sig1=sigma[1], sig2=sigma[2], rho=rho[1, 2])
})

cbind(Y.tab, PLN1.tab, PLN2.tab)
plot(PLN2.tab, PLN1.tab-PLN2.tab); abline(h=0)
points(PLN2.tab, qbinom(.025, B, PLN2.tab)/B-PLN2.tab, pch='+', col=4)
points(PLN2.tab, qbinom(.975, B, PLN2.tab)/B-PLN2.tab, pch='+', col=4)

invisible(sapply(1:p, function(j){
   cat('Y', j, '\n')
   P1 = as.vector(by(PLN1.tab, Y.tab[, j], sum))
   P2 = as.vector(by(PLN2.tab, Y.tab[, j], sum))
   L2 = qbinom(.025, B, P2)/B
   U2 = qbinom(.975, B, P2)/B
   print(cbind(sort(unique(Y.tab[, j])), P2, L2, P1, U2))
}))

##################################################################################
p = 3
##################################################################################
# Y table
Y.tab = expand.grid(0:y.max, 0:y.max, 0:y.max); nb.tab = nrow(Y.tab)

# PLN1
Z = rmvnorm(B, sigma=Sigma)
dpois.Y.tab = array(0, dim=c(B, (1+y.max), p))
invisible(sapply(1:p, function(j){
   cat(j, '')
   sapply(0:y.max, function(y){dpois.Y.tab[, (1+y), j] <<- dpois(y, exp(Z[, j]))})
   }))
PLN1.tab = sapply(1:nb.tab, function(i){
   mean(dpois.Y.tab[, (1+Y.tab[i, 1]), 1] * dpois.Y.tab[, (1+Y.tab[i, 2]), 2] *
           dpois.Y.tab[, (1+Y.tab[i, 3]), 3])
})

# PLN2
PLN2.tab = sapply(1:nb.tab, function(i){
   dbipoilog(Y.tab[i, 1], Y.tab[i, 2], mu1=0, mu2=0, sig1=sigma[1], sig2=sigma[2], rho=rho[1, 2]) *
      dbipoilog(Y.tab[i, 1], Y.tab[i, 3], mu1=0, mu2=0, sig1=sigma[1], sig2=sigma[3], rho=rho[1, 3]) /
      dpoilog(Y.tab[i, 1], mu=0, sig=sigma[1])
})

cbind(Y.tab, PLN1.tab, PLN2.tab)
plot(PLN2.tab, PLN1.tab-PLN2.tab); abline(h=0)
points(PLN2.tab, qbinom(.025, B, PLN2.tab)/B-PLN2.tab, pch='+', col=4)
points(PLN2.tab, qbinom(.975, B, PLN2.tab)/B-PLN2.tab, pch='+', col=4)

P1.1 = P2.1 = L2.1 = U2.1 = c()
invisible(sapply(1:p, function(j){
   cat('Y', j, '\n')
   P1 = as.vector(by(PLN1.tab, Y.tab[, j], sum))
   P2 = as.vector(by(PLN2.tab, Y.tab[, j], sum))
   L2 = qbinom(.025, B, P2)/B
   U2 = qbinom(.975, B, P2)/B
   P1.1 <<- c(P1.1, P1); P2.1 <<- c(P2.1, P2); L2.1 <<- c(L2.1, L2); U2.1 <<- c(U2.1, U2)
   print(cbind(sort(unique(Y.tab[, j])), P2, L2, P1, U2))
}))
plot(P2.1, P1.1-P2.1); abline(h=0)
points(P2.1, L2.1-P2.1, pch='+', col=4); 
points(P2.1 ,U2.1-P2.1, pch='+', col=4)

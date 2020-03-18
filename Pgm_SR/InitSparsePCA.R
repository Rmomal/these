# Using sparse PCA (sPCA) for the intialization of EMtree with a missing actor
# Notes on sPCA : 
# scores = data %*% loadings
# approximation od the data = scores %*% t(transform)


rm(list=ls()); par(pch=20)
library(sparsepca); library(mvtnorm); library(PLNmodels)

###############################################################################
# Function
FitSparsePCA <- function(Y, q=1, alphaGrid=10^(seq(-4, 0, by=.1))){
   # Fit sparse PCA on a grid of alpha
   # Needs an estimate o Sigma: empirical variance of the rotated scores 
   #    + diagonal risidual variances (why not?)
   n <- nrow(Y); p <- nrow(Y); alphaNb <- length(alphaGrid); 
   # Fits all sparsePCA
   sPCA <- list()
   for(a in 1:alphaNb){
      sPCA[[a]] <- spca(Y, k=q, alpha=alphaGrid[a])
      sPCA[[a]]$Sigma <- cov(sPCA[[a]]$scores%*%t(sPCA[[a]]$transform))
      resVar <- (n-1)*apply(Y - sPCA[[a]]$scores %*% t(sPCA[[a]]$transform), 2, var)/n
      sPCA[[a]]$Sigma <- sPCA[[a]]$Sigma  + diag(resVar)
      sPCA[[a]]$df <- 1 + sum(sPCA[[a]]$loadings!=0)
      sPCA[[a]]$loglik <- sum(dmvnorm(scale(Y), sigma=sPCA[[a]]$Sigma, log=TRUE))
      sPCA[[a]]$bic <- sPCA[[a]]$loglik - log(n)*sPCA[[a]]$df/2
   }
   # Selects alpha via pseudo-BIC
   loglik <- unlist(lapply(sPCA, function(sPca){sPca$loglik}))
   bic <- unlist(lapply(sPCA, function(sPca){sPca$bic}))
   aOpt <- which.max(bic)
   # Find the cliques
   alphaOpt <- alphaGrid[aOpt]
   sPcaOpt <- sPCA[[aOpt]]
   sPcaOpt$loadings
   cliques <- list(); 
   sapply(1:ncol(sPcaOpt$loadings), function(j){cliques[[j]] <<- which(sPcaOpt$loadings[, j]!=0)})
   return(list(sPcaOpt=sPcaOpt, alphaGrid=alphaGrid, alphaOpt=alphaOpt, loglik=loglik, bic=bic, cliques=cliques))
}

###############################################################################
# Simulated data
n <- 50; q0 <- 3; r <- 3; p <- q0*r
U <- 3*matrix(rnorm(n*q0), n, q0)
Y <- sapply(1:p, function(j){U[, 1+floor((j-1)/r)] +  rnorm(n)})
# Y <- cbind(matrix(rnorm(n*r), n, r)+U[, 1]%o%rep(1, r), matrix(rnorm(n*(p-r)), n, (p-r))); q0 <- 1
Y <- scale(Y)
image(1:p, 1:p, cov(Y))
Sigma <- cov(Y)
# Finds cliques
sPCA <- FitSparsePCA(Y)
plot(sPCA$alphaGrid, sPCA$loglik, ylim=c(min(sPCA$bic), max(sPCA$loglik)), type='b', log='x')
points(sPCA$alphaGrid, sPCA$bic, type='b', col=2)
abline(v=sPCA$alphaOpt, col=2, lty=2)
sPCA$cliques

###############################################################################
# Simulated data
# Fits PLN + Creates fake Gaussian latents
load('../Data_SR/BarentsFish_Group.Rdata')
# load('../Data_SR/BarentsFish.Rdata')
PLNfit <- PLN(Data$count ~ 1)
Y <- PLNfit$var_par$M; S <- PLNfit$var_par$S
sapply(1:nrow(Y), function(i){Y[i, ] <<- Y[i, ] + rmvnorm(1, sigma=diag(S[i, ]))})
# Y <- PLNfit$var_par$M
Y <- scale(Y)

sPCA <- FitSparsePCA(Y, q=3)
plot(sPCA$alphaGrid, sPCA$loglik, ylim=c(min(sPCA$bic), max(sPCA$loglik)), type='b', log='x')
points(sPCA$alphaGrid, sPCA$bic, type='b', col=2)
abline(v=sPCA$alphaOpt, col=2, lty=2)
sPCA$cliques


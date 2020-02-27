n <- 50; p <- 10; 
cliqueList <- list(); cliqueList[[1]] <- (1:3); cliqueList[[2]] <- (6:7)
cliqueNb <- length(cliqueList)

# Simul de Sigma
sigma <- rgamma(p, 2, 1)
Sigma <- exp(-as.matrix(dist(matrix(rnorm(2*p), p, 2), diag=TRUE)))
Sigma <- diag(sigma)%*%Sigma%*%diag(sigma)

# Init hidden node
# principal components associates with each cliques are not orthogonal ...
Corr <- cov2cor(Sigma); sigma <- sqrt(diag(Sigma))
coef <- matrix(0, p, cliqueNb); lambda <- rep(0, cliqueNb)
sapply(1:cliqueNb, function(c){
   pca <- eigen(cov2cor(Corr[cliqueList[[c]], cliqueList[[c]]]))
   coef[, c] <<- rep(0, p); 
   coef[cliqueList[[c]], c] <<- pca$vectors[, 1]
   lambda[c] <<- pca$values[1]
}) 

# Recontructing Sigma
CorrFull <- rbind(cbind(Corr, Corr%*%coef), cbind(t(coef)%*%Corr, 1.1*t(coef)%*%Corr%*%coef))
isSymmetric(CorrFull); eigen(CorrFull)$values; plot(CorrFull[1:p, 1:p], Corr); abline(0, 1)
sigmaFull <- c(sigma, rep(1, cliqueNb)) 
SigmaFull <- diag(sigmaFull) %*% CorrFull %*% diag(sigmaFull)
isSymmetric(SigmaFull); eigen(SigmaFull)$values; plot(SigmaFull[1:p, 1:p], Sigma); abline(0, 1)

# Initialising Omega
OmegaFull <- solve(SigmaFull)
isSymmetric(OmegaFull); eigen(OmegaFull)$values
# sapply(1:cliqueNb, function(c){
#    OmegaFull[(1:p)[-cliqueList[[c]]], (p+c)] <<- OmegaFull[(p+c), (1:p)[-cliqueList[[c]]]] <<- 0
# })
# isSymmetric(OmegaFull); eigen(OmegaFull)$values
coefDiag <- c(rep(1, p), 1/sqrt(diag(OmegaFull[p+(1:cliqueNb), p+(1:cliqueNb)])))
OmegaFull <- diag(coefDiag) %*% OmegaFull %*% diag(coefDiag)

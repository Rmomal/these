# Using sparse PCA (sPCA) for the intialization of EMtree with a missing actor
# Notes on sPCA : 
# scores = data %*% loadings
# approximation od the data = scores %*% t(transform)


rm(list=ls()); par(pch=20)
library(sparsepca); library(mvtnorm); library(PLNmodels)

###############################################################################
# Function
FitSparsePCA <- function(Y, q=1, alphaGrid=10^(seq(-4, 0, by=.1)) ){
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
    sPCA[[a]]$loglik <- sum(dmvnorm((Y), sigma=sPCA[[a]]$Sigma, log=TRUE))
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
  loadings=do.call(cbind,lapply(sPCA, function(x){x$loadings}))
  sapply(1:ncol(sPcaOpt$loadings), function(j){cliques[[j]] <<- which(sPcaOpt$loadings[, j]!=0)})
  return(list(loadings=loadings,sPcaOpt=sPcaOpt, alphaGrid=alphaGrid, alphaOpt=alphaOpt, loglik=loglik, bic=bic, cliques=cliques))
}
verdictSpca<-function(sPCA, trueClique){
  plot(sPCA$alphaGrid, sPCA$loglik, ylim=c(min(sPCA$bic), max(sPCA$loglik)), pch=20, log='x')
  points(sPCA$alphaGrid, sPCA$bic, pch=20, col=2)
  abline(v=sPCA$alphaOpt, col=2, lty=2)
  sPCA$cliques
  N=setdiff(1:p,trueClique)
  FP=sum(sPCA$cliques[[1]]%in%N)/length(N)
  FN=sum(setdiff(1:p,sPCA$cliques[[1]])%in%trueClique)/length(trueClique)
  c(FN, FP)
  c(length(sPCA$cliques[[1]]), length(trueClique))
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
#######
set.seed(10)
n=200
p=14
r=1
type="scale-free"
data=data_from_scratch(type = type,p = p+r,n = n,signed = FALSE,prob = 5/p,v = 0)
omega=data$omega
hidden=which(diag(omega)%in%sort(diag(omega), decreasing = TRUE)[1:r])[1:r] # on cache les r plus gros
trueClique=which(omega[hidden,-hidden]!=0)
Kh  <- omega[hidden,hidden]
Ko  <- omega[-hidden,-hidden]
Koh <- omega[-hidden,hidden]
Km  <- Ko - Koh %*%solve(Kh)%*% t(Koh)
sigmaO=solve(Km)
# counts generation
counts=generator_PLN(sigmaO,covariates = NULL,n=n)
PLNfit = PLN(counts~1, control = list(trace=FALSE))
Y <- PLNfit$var_par$M; S <- PLNfit$var_par$S
sapply(1:nrow(Y), function(i){Y[i, ] <<- Y[i, ] + rmvnorm(1, sigma=diag(S[i, ]))})


sPCA <- FitSparsePCA((Y))
varySize<-colSums(1*(sPCA$loadings!=0))
plot(sPCA$alphaGrid,varySize, pch=20, log='x')
abline(v=sPCA$alphaOpt, col=2, lty=2)

plot(sPCA$alphaGrid, sPCA$loglik, ylim=c(min(sPCA$bic), max(sPCA$loglik)), pch=20, log='x')
points(sPCA$alphaGrid, sPCA$bic, pch=20, col=2)
abline(v=sPCA$alphaOpt, col=2, lty=2)
sPCA$cliques
max(sPCA$bic)
N=setdiff(1:p,trueClique)
FP=sum(sPCA$cliques[[1]]%in%N)/length(N)
FN=sum(setdiff(1:p,sPCA$cliques[[1]])%in%trueClique)/length(trueClique)
c(FN, FP)
c(length(sPCA$cliques[[1]]), length(trueClique))

clique=which(sPCA$loadings[,30]!=0)
FP=sum(clique%in%N)/length(N)
FN=sum(setdiff(1:p,clique)%in%trueClique)/length(trueClique)
c(FN, FP)
###############################################################################
# Simulated data
# Fits PLN + Creates fake Gaussian latents
load('../Data_SR/BarentsFish_Group.Rdata')
load("/Users/raphaellemomal/these/Data_SR/BarentsFish.Rdata")
# load('../Data_SR/BarentsFish.Rdata')
PLNfit <- PLN(Data$count ~ 1)
Y <- PLNfit$var_par$M; S <- PLNfit$var_par$S
sapply(1:nrow(Y), function(i){Y[i, ] <<- Y[i, ] + rmvnorm(1, sigma=diag(S[i, ]))})
# Y <- PLNfit$var_par$M
Y <- scale(Y)

q <- 1
sPCA <- FitSparsePCA(Y, q=q)
plot(sPCA$alphaGrid, sPCA$loglik, ylim=c(min(sPCA$bic), max(sPCA$loglik)), type='b', log='x')
points(sPCA$alphaGrid, sPCA$bic, type='b', col=2)
abline(v=sPCA$alphaOpt, col=2, lty=2)
sPCA$cliques
cor(Data$covariates)
cor(cbind(sPCA$sPcaOpt$scores, Data$covariates))[1:q, q+(1:ncol(Data$covariates))]

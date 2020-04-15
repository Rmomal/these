# Model selection for PLN tree-shaped models

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
library(PLNmodels); library(EMtree); library(LITree); library(mvtnorm)
source("/home/robin/RECHERCHE/RESEAUX/R-Momal/these/R/codes/missingActor/fonctions-missing.R")

# Barents data
dataDir <- '../Data_SR/'; dataName <- 'BarentsFish_Group'
load(paste0(dataDir, dataName, '.Rdata'))
Y <- Data$count; n <- nrow(Y); p <- ncol(Y); d <- ncol(Data$covariates) 
logYfact <- sum(lgamma(Y+1)); log2pi <- log(2*acos(-1))

# Function
F_Crit <- function(Y, X, theta, Sigma, Omega=solve(Sigma), M, S, beta=NULL, betat=NULL, Pt=EdgeProba(betat)){
   # Y=Y; X=X; theta=t(pln$model_par$Theta); Sigma=pln$model_par$Sigma; M=pln$var_par$M; S=pln$var_par$S
   n <- nrow(Y); p <- ncol(Y); q <- ncol(M); r <- q-p; d <- ncol(X)
   # Y
   EqLogpY <- -sum(exp(X%*%theta + M[, 1:p] + S[, 1:p]/2)) + sum(Y*(X%*%theta + M[, 1:p])) - logYfact
   dfY <- prod(dim(theta))
   if(is.null(beta)){ # PLN
      # T
      EqLogpT <- HqT <- dfT <- 0
      # Z
      EqLogpZ <- -n*log(det(Sigma))/2 - sum((t(M)%*%M+diag(colSums(S))) * Omega)/2 - n*q*log2pi/2
      dfZ <- p*(p+1)/2
      HqZ <- sum(log(S))/2 + n*q*(1+log2pi)/2
   }else{ # EM tree or miss
      # T
      EqLogpT <- sum(Pt*log(beta + (beta==0)))/2 - log(SumTree(beta))
      dfT <- q*(q-1)/2 - 1
      HqT <- -sum(Pt*log(betat + (betat==0)))/2 + log(SumTree(betat))
      # Z
      logPsi <- log(1 - cov2cor(Omega)^2 + (cov2cor(Omega)==1))/2
      EqLogpZ <- n*sum(log(diag(Omega)))/2 + n*sum(Pt*logPsi)/2 - 
         sum(Pt * (t(M)%*%M+diag(colSums(S))) * Omega)/2 - n*q*log2pi/2
      dfZ <- q*(q+1)/2
      HqZ <- sum(log(S))/2 + n*q*(1+log2pi)/2
   }
   EqLogpTZY <- EqLogpT + EqLogpZ + EqLogpY
   JZ <- EqLogpZ+HqZ
   HqTZ <- HqT + HqZ
   JY <- EqLogpTZY + HqTZ
   pBIC <- (dfT+dfZ+dfY)*log(n)/2
   return(list(d=d, r=r, EqLogpT=EqLogpT, EqLogpZ=EqLogpZ, EqLogpY=EqLogpY, dfT=dfT, dfZ=dfZ, dfY=dfY, HqT=HqT, HqZ=HqZ, 
               EqLogpTZY=EqLogpTZY, JZ=JZ, HqTZ=HqTZ, JY=JY, 
               pBIC=pBIC, BIC=JY-pBIC, ICLT=EqLogpTZY-HqT, ICLZ=EqLogpTZY-HqZ, ICLTZ=EqLogpTZY-HqTZ))
}

# PLN fit
if(!file.exists(paste0(dataDir, dataName, '-plnList.Rdata'))){
   XList <- sapply(0:d, function(h){X <- matrix(1, nrow(Y), 1); if(h>0){X <- cbind(X, Data$covariates[, 1:h])}; return(X)})
   plnList <- sapply(0:d, function(h){PLN(Y ~ -1+XList[[1+h]])})
   save(XList, plnList, file=paste0(dataDir, dataName, '-plnList.Rdata'))
}
load(paste0(dataDir, dataName, '-plnList.Rdata'))

# PLN crit
critPln <- t(sapply(0:d, function(h){ # h <- 2
   pln <- plnList[[1+h]]; X <- XList[[1+h]]
   F_Crit(Y=Y, X=X, theta=t(pln$model_par$Theta), Sigma=pln$model_par$Sigma, M=pln$var_par$M, S=pln$var_par$S)
   }))

# EMtree fit
if(!file.exists(paste0(dataDir, dataName, '-emTreeList.Rdata'))){
   emTreeList <- list()
   sapply(0:d, function(h){emTreeList[[1+h]] <<- EMtree(plnList[[1+h]])})
   save(emTreeList, file=paste0(dataDir, dataName, '-emTreeList.Rdata'))
}
load(paste0(dataDir, dataName, '-emTreeList.Rdata'))

# EMtree crit
critEmTree <- t(sapply(0:d, function(h){ # h <- 2
   pln <- plnList[[1+h]]; emTree <- emTreeList[[1+h]]; X <- XList[[1+h]]
   Sigma <- pln$model_par$Sigma; beta <- emTree$edges_weight
   logPsi <- -log(1 - cov2cor(Sigma)^2 + (cov2cor(Sigma)==1))/2
   betat <- beta * exp(n*logPsi)
   F_Crit(Y=Y, X=X, theta=t(pln$model_par$Theta), Sigma=pln$model_par$Sigma, M=pln$var_par$M, S=pln$var_par$S, 
          beta=beta, betat=betat)
}))

# EMmiss fit
rMax <- 2; 
if(!file.exists(paste0(dataDir, dataName, '-emMissList.Rdata'))){
   library(tidyverse); library(useful); library(MASS); library(reshape2); library(parallel); library(sparsepca)
   eps <- 1e-3; alpha <- 1; cores <- 3; plot <- FALSE; maxIter <- 100; B <- 1
   emMissList <- list()
   for(h in 0:d){ # h <- 1
      emMissList[[1+h]] <- list()
      MO <-plnList[[1+h]]$var_par$M ; SO <- plnList[[1+h]]$var_par$S; 
      theta <- plnList[[1+h]]$model_par$Theta; SigmaO <- plnList[[1+h]]$model_par$Sigma
      #initialize 
      cat('\n init \n')
      init0 = initVEM(Y, initviasigma=NULL, SigmaO, r=0)
      Wginit= init0$Wginit; Winit= init0$Winit; Omegainit=init0$omegainit 
      # vem r 0
      cat('\n', c(h, 0), '\n')
      emMissList[[1+h]]$VEM[[1]]
      VEM0 <-VEMtree(Y, MO, SO, MH=NULL, ome_init=Omegainit, W_init=Winit, eps=eps,
                     Wg_init=Wginit, plot=plot, maxIter=maxIter, print.hist=FALSE,
                     vraiOm=NULL, alpha=alpha, verbatim=FALSE, filterPg=TRUE, 
                     filterWg=TRUE)
      emMissList[[1+h]]$VEM0 <- VEM0
      #  vary r
      VEMr<-lapply(1:rMax, function(r){
         cat('\n', c(h, r), '\n')
         cat(paste0(r," missing actors: "))
         cliques_spca <- boot_FitSparsePCA(scale(Y),B,r=r)
         ListVEM <- List.VEM(cliqueList=cliques_spca, Y, SigmaO, MO, SO, alpha, r=r, cores=cores,
                             eps=eps, maxIter)
         return(ListVEM)
      })
      emMissList[[1+h]]$VEMr <- VEMr
   }
   save(emMissList, file=paste0(dataDir, dataName, '-emMissList.Rdata'))
}
load(paste0(dataDir, dataName, '-emMissList.Rdata'))

# Reshape EmMiss fit
for(h in 0:d){ # h <- 1
   emMissList[[1+h]]$VEM <- list()
   emMissList[[1+h]]$VEM[[1]] <- emMissList[[1+h]]$VEM0
   for(r in 1:rMax){cat(r, ''); emMissList[[1+h]]$VEM[[1+r]] <- emMissList[[1+h]]$VEMr[[r]][[1]]}
}

# EMmiss crit
critEmMiss <- c()
invisible(t(sapply(0:d, function(h){ # h <- 2
   pln <- plnList[[1+h]]; X <- XList[[1+h]]; emMiss <- emMissList[[1+h]]
   sapply(0:rMax, function(r){ # r <- 1
      vem <- emMiss$VEM[[1+r]]; 
      crit <- unlist(F_Crit(Y=Y, X=X, theta=t(pln$model_par$Theta), Sigma=NULL, Omega=vem$omega, M=vem$M, S=vem$S, 
                            beta=vem$W, betat=vem$Wg, Pt=EdgeProba(vem$W)))
      crit <- unlist(c(crit, vem$lowbound[nrow(vem$lowbound), ]))
      critEmMiss <<- rbind(critEmMiss, crit)
   })
})))

# critEmMiss <- critOrg
critEmMiss <- as.data.frame(critEmMiss)
critEmMiss$klT <- critEmMiss$EqLogpT + critEmMiss$HqT 
print(critEmMiss[, c('d', 'r', 'EqLogpT', 'HqT', 'T2', 'klT', 'JY', 'J')])
critEmMiss$JZ <- critEmMiss$EqLogpZ + critEmMiss$HqZ 
critEmMiss$T3 <- critEmMiss$J - critEmMiss$T1 - critEmMiss$T2 
critEmMiss$T13 <- critEmMiss$T1 + critEmMiss$T3 
print(critEmMiss[, c('d', 'r', 'EqLogpZ', 'HqZ', 'JZ', 'T1', 'T3', 'T13')])
critEmMiss$JTZ <- critEmMiss$EqLogpT + critEmMiss$HqT + critEmMiss$EqLogpZ + critEmMiss$HqZ 
print(critEmMiss[, c('d', 'r', 'EqLogpT', 'EqLogpZ', 'EqLogpY', 'HqT', 'HqZ', 'JTZ', 'J')])
# print(critEmMiss[1:6, c('d', 'r', 'EqLogpT', 'EqLogpZ.T', 'EqLogpY.Z', 'HqT', 'HqZ')])
# modeSelEmMiss <- critEmMiss[, c('d', 'r', 'EqLogpTZY', 'HqTZ', 'J', 'pBIC', 'BIC', 'ICL')]
# print(modeSelEmMiss)

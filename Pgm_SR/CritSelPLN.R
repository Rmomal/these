# Model selection for PLN tree-shaped models

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
library(PLNmodels); library(EMtree); library(LITree); library(mvtnorm)
# source("/home/robin/RECHERCHE/RESEAUX/R-Momal/these/R/codes/missingActor/fonctions-missing.R")
source("/home/robin/RECHERCHE/RESEAUX/R-Momal/these/R/codes/missingActor/fonctions-missing-SR.R")
# source("/home/robin/RECHERCHE/RESEAUX/R-Momal/these/R/codes/missingActor/fonctions-missing-14-04-20.R")

################################################################################
# Functions
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
      EqLogpZ <- -n*log(det(Sigma))/2 - sum((t(M)%*%M+diag(colSums(S))) * Omega)/2 - n*p*log2pi/2
      dfZ <- p*(p+1)/2
      HqZ <- sum(log(S))/2 + n*p*(1+log2pi)/2
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
   JT <- EqLogpT + HqT
   JZ <- EqLogpZ + HqZ
   HqTZ <- HqT + HqZ
   JTZY <- EqLogpTZY + HqTZ
   pBIC <- (dfT+dfZ+dfY)*log(n)/2
   # pBIC <- dfZ*n/q*log(n)*2 + dfY*log(n)/2
   return(list(d=d, r=r, 
               EqLogpT=EqLogpT, EqLogpZ=EqLogpZ, EqLogpY=EqLogpY, 
               dfT=dfT, dfZ=dfZ, dfY=dfY, HqT=HqT, HqZ=HqZ, 
               JT=JT, JZ=JZ, EqLogpTZY=EqLogpTZY, HqTZ=HqTZ, JTZ=JT+JZ, JTZY=JTZY, 
               pBIC=pBIC, BIC=JTZY-pBIC, ICLT=EqLogpTZY-HqT, ICLZ=EqLogpTZY-HqZ, ICLTZ=EqLogpTZY-HqTZ))
}
F_PlotSlopeCrit <- function(pen, crit, r, d, title=''){
   # pen=critEmMiss$pBIC; crit=critEmMiss$JTZY; r=critEmMiss$r; d=critEmMiss$d; title=''
   reg <- t(sapply(1:max(d), function(dd){lm(crit[which(d==dd)] ~ pen[which(d==dd)])$coef}))
   plot(pen, crit, col=d, pch=1+r, main=title)
   sapply(1:max(d), function(dd){
      abline(reg[dd, 1], reg[dd, 2], col=dd, lwd=2)
      # yLeg <- min(crit)+(max(d):1)/(max(d)+1)*(max(crit)-min(crit))
      # text(min(pen[which(d==dd)]), yLeg[dd], round(reg[dd, 2], 1), col=dd, cex=1, pos=3)
      xLeg <- min(pen)+(1:max(d))/(max(d)+1)*(max(pen)-min(pen))
      text(xLeg[dd], max(crit), round(reg[dd, 2], 1), col=dd, cex=1, pos=1)
   })
}

################################################################################
# Data
# Barents data
dataDir <- '../Data_SR/'
#dataDir <- '/Users/raphaellemomal/these/Data_SR/'; dataName <- 'BarentsFish_Group'
load(paste0(dataDir, dataName, '.Rdata'))
Y <- Data$count; 
# Y <- Data$count[, 1:round(ncol(Y)/2)]; dataName <- paste0(dataName, '_1stHalf')
# Y <- Data$count[, (round(ncol(Y)/2)+1):ncol(Y)]; dataName <- paste0(dataName, '_2ndHalf')

# # Fake Barents
# seed <- 1; set.seed(seed); fakeName <- paste0('fakeSeed', seed, '-', dataName)
# if(!file.exists(paste0(dataDir, fakeName, '.Rdata'))){
#    load(paste0(dataDir, dataName, '-plnList.Rdata'))
#    pln <- plnList[[1]]; X <- XList[[1]]
#    Z <- rmvnorm(nrow(X), sigma=diag(diag(pln$model_par$Sigma)))
#    mu <- X%*%t(pln$model_par$Theta) + Z
#    Y <- matrix(rpois(prod(dim(Z)), exp(mu)), nrow(Z), ncol(Z))
#    Data <- list(count=Y, covariates=Data$covariates)
#    save(Data, file=paste0(dataDir, fakeName, '.Rdata'))
# }
# dataName <- fakeName 
# load(paste0(dataDir, dataName, '.Rdata')); Y <- Data$count; 
# plot(mu, log(1+Y)); abline(0, 1)

# Dimensions
n <- nrow(Y); p <- ncol(Y); dMax <- 1+ncol(Data$covariates) 
rMax <- 3; logYfact <- sum(lgamma(Y+1)); log2pi <- log(2*acos(-1))

################################################################################
# PLN fit
if(!file.exists(paste0(dataDir, dataName, '-plnList.Rdata'))){
   XList <- sapply(1:dMax, function(d){X <- matrix(1, nrow(Y), 1); if(d>1){X <- cbind(X, Data$covariates[, 1:(d-1)])}; return(X)})
   plnList <- sapply(1:dMax, function(d){PLN(Y ~ -1+XList[[d]])})
   save(XList, plnList, file=paste0(dataDir, dataName, '-plnList.Rdata'))
}
load(paste0(dataDir, dataName, '-plnList.Rdata'))

# PLN crit
critPln <- t(sapply(1:dMax, function(d){ # d <- 1
   pln <- plnList[[d]]; X <- XList[[d]]
   F_Crit(Y=Y, X=X, theta=t(pln$model_par$Theta), Sigma=pln$model_par$Sigma, M=pln$var_par$M, S=pln$var_par$S)
}))

################################################################################
# EMtree fit
if(!file.exists(paste0(dataDir, dataName, '-emTreeList-r', rMax, '.Rdata'))){
   emTreeList <- list()
   sapply(1:dMax, function(d){emTreeList[[d]] <<- EMtree(plnList[[d]])})
   save(emTreeList, file=paste0(dataDir, dataName, '-emTreeList-r', rMax, '.Rdata'))
}
load(paste0(dataDir, dataName, '-emTreeList-r', rMax, '.Rdata'))

# EMtree crit
critEmTree <- t(sapply(1:dMax, function(d){ # d <- 2
   pln <- plnList[[d]]; emTree <- emTreeList[[d]]; X <- XList[[d]]
   Sigma <- pln$model_par$Sigma; beta <- emTree$edges_weight
   logPsi <- -log(1 - cov2cor(Sigma)^2 + (cov2cor(Sigma)==1))/2
   betat <- beta * exp(n*logPsi)
   F_Crit(Y=Y, X=X, theta=t(pln$model_par$Theta), Sigma=pln$model_par$Sigma, M=pln$var_par$M, S=pln$var_par$S, 
          beta=beta, betat=betat)
}))

################################################################################
# EMmiss fit
noBeta = FALSE

if(!file.exists(paste0(dataDir, dataName, '-emMissList-r', rMax, '.Rdata'))){
   library(tidyverse); library(useful); library(MASS); library(reshape2); library(parallel); library(sparsepca)
   eps <- 1e-3; alpha <- 1; cores <- 3; plot <- FALSE; maxIter <- 100; B <- 1
   emMissList <- list()
   for(d in 1:dMax){ # d <- 1
      emMissList[[d]] <- list()
      MO <-plnList[[d]]$var_par$M ; SO <- plnList[[d]]$var_par$S; 
      theta <- plnList[[d]]$model_par$Theta; SigmaO <- plnList[[d]]$model_par$Sigma
      #initialize 
      cat('\n init \n')
      init0 = initVEM(Y, initviasigma=NULL, SigmaO, r=0)
      Wginit= init0$Wginit; Winit= init0$Winit; Omegainit=init0$omegainit 
      # vem r 0
      cat('\n', c(d, 0), '\n')
      emMissList[[d]]$VEM[[1]]
      q=p+r
      D=.Machine$double.xmax
      alpha = (1/n)*((1/(q-1))*log(D) - log(q))
      VEM0 <-VEMtree(Y, MO, SO, MH=NULL, ome_init=Omegainit, W_init=Winit, eps=eps,
                     Wg_init=Wginit, plot=TRUE, maxIter=maxIter, print.hist=FALSE,
                     alpha=alpha, verbatim=FALSE,   nobeta=noBeta)
      emMissList[[d]]$VEM0 <- VEM0
      #  vary r
      VEMr<-lapply(1:rMax, function(r){
         cat('\n', c(d, r), '\n')
         cat(paste0(r," missing actors: "))
         cliques_spca <- boot_FitSparsePCA(scale(Y),B,r=r)
         ListVEM <- List.VEM(cliquesObj=cliques_spca, Y, SigmaO, MO, SO,  r=r, cores=1,
                             eps=eps, maxIter=maxIter,  nobeta=noBeta)
         return(ListVEM)
      })
      emMissList[[d]]$VEMr <- VEMr
   }
   save(emMissList, file=paste0(dataDir, dataName, '-emMissList-r', rMax, '.Rdata'))
}
load(paste0(dataDir, dataName, '-emMissList-r', rMax, '.Rdata'))

# Reshape EmMiss fit
for(d in 1:dMax){ # d <- 1
   emMissList[[d]]$VEM <- list()
   emMissList[[d]]$VEM[[1]] <- emMissList[[d]]$VEM0
   for(r in 1:rMax){cat(r, ''); emMissList[[d]]$VEM[[1+r]] <- emMissList[[d]]$VEMr[[r]][[1]]}
}

# EMmiss crit
critEmMiss <- c()
invisible(t(sapply(1:dMax, function(d){ # d <- 2
   pln <- plnList[[d]]; X <- XList[[d]]; emMiss <- emMissList[[d]]
   sapply(0:rMax, function(r){ # r <- 1
      vem <- emMiss$VEM[[1+r]]; 
      crit <- unlist(F_Crit(Y=Y, X=X, theta=t(pln$model_par$Theta), Sigma=NULL, Omega=vem$omega, M=vem$M, S=vem$S, 
                            beta=vem$W, betat=vem$Wg))
      crit <- unlist(c(crit, vem$lowbound[nrow(vem$lowbound), ]))
      critEmMiss <<- rbind(critEmMiss, crit)
   })
})))

################################################################################
# Results
# Comp parms EmTree / EmMiss
par(mfrow=c(dMax, 5), mex=.6)
# par(mfrow=c(3, 2), mex=.6)
t(sapply(1:dMax, function(d){ # d <- 1
   pln <- plnList[[d]]; X <- XList[[d]]; emTree <- emTreeList[[d]]; emMiss <- emMissList[[d]]$VEM0
   Sigma <- pln$model_par$Sigma; logPsi <- -log(1 - cov2cor(Sigma)^2 + (cov2cor(Sigma)==1))/2;
   # beta
   betaTree <- emTree$edges_weight; betaMiss <- emMiss$W
   betatTree <- betaTree * exp(n*logPsi)
   betatMiss <- emMiss$Wg
   logBetaTree <- log(betaTree) - log(SumTree(betaTree))/(p-1) + log((p-1)) 
   logBetatTree <- log(betatTree) - log(SumTree(betatTree))/(p-1) + log((p-1)) 
   logBetaMiss <- log(betaMiss) - log(SumTree(betaMiss))/(p-1) + log((p-1)) 
   logBetatMiss <- log(betatMiss) - log(SumTree(betatMiss))/(p-1) + log((p-1)) 
   c(SumTree(exp(logBetaTree)), SumTree(exp(logBetatTree)), SumTree(exp(logBetaMiss)), SumTree(exp(logBetatMiss)))
   plot(logBetatTree, logBetatMiss, col=2, xlab='tree', ylab='miss', main=paste0('beta: d=', d)); abline(a=0, b=1, v=0)
   points(logBetaTree, logBetaMiss, col=1)
   points(logBetatTree-logBetaTree, logBetatMiss-logBetaMiss, col=4)
   # M
   plot(pln$var_par$M, emMiss$M, xlab='tree', ylab='miss', main=paste0('M: d=', d)); abline(a=0, b=1, v=0)
   # S
   plot(pln$var_par$S, emMiss$S, xlab='tree', ylab='miss', main=paste0('S: d=', d)); abline(a=0, b=1, v=0)
   # Omega
   colMat <- matrix(1, p, p); diag(colMat) <- 2
   plot(solve(Sigma), emMiss$omega, col=colMat, xlab='tree', ylab='miss', main=paste0('Omega: d=', d)); abline(a=0, b=1, h=0, v=0)
   # Z.hat * Omega
   colMat <- matrix(1, p, p); diag(colMat) <- 2
   SShatTree <- (t(pln$var_par$M)%*%pln$var_par$M+diag(colSums(pln$var_par$S)))
   SShatMiss <- (t(emMiss$M)%*%emMiss$M+diag(colSums(emMiss$S)))
   plot(SShatTree*solve(Sigma), SShatMiss*emMiss$omega, col=colMat, xlab='tree', ylab='miss', main=paste0('SS*Omega: d=', d)); abline(a=0, b=1, h=0, v=0)
   # # Sigma
   # colMat <- matrix(1, p, p); diag(colMat) <- 2
   # plot((Sigma), (solve(emMiss$omega)), col=colMat, xlab='tree', ylab='miss', main=paste0('Sigma: d=', d)); abline(a=0, b=1, h=0, v=0)
}))

# Check RM / SR
critEmMiss <- as.data.frame(critEmMiss)
print(critEmMiss[, c('d', 'r', 'EqLogpT', 'HqT', 'T2', 'JT', 'JTZY', 'J')])
critEmMiss$T3 <- critEmMiss$J - critEmMiss$T1 - critEmMiss$T2 
critEmMiss$T13 <- critEmMiss$T1 + critEmMiss$T3 
print(critEmMiss[, c('d', 'r', 'EqLogpZ', 'HqZ', 'JZ', 'T1', 'T3', 'T13')])
print(critEmMiss[, c('d', 'r', 'EqLogpT', 'EqLogpZ', 'EqLogpY', 'HqT', 'HqZ', 'JTZ', 'J')])

critEmMiss$ICLT <- critEmMiss$BIC - critEmMiss$HqT
critEmMiss$ICLZ <- critEmMiss$BIC - critEmMiss$HqZ
critEmMiss$ICLTZ <- critEmMiss$BIC - critEmMiss$HqTZ
critEmMiss$q <- p + critEmMiss$r
critEmMiss$cardT <- (critEmMiss$q-2) * log(critEmMiss$q)

varList <- c('d', 'r', 'EqLogpT', 'HqT', 'JT', 'EqLogpZ', 'HqZ', 'JZ', 'EqLogpY', 'JTZY', 'JTZ') #, 'pBIC', 'BIC', 'ICLTZ')
critEmMiss[, varList]
critEmMiss[which(critEmMiss$r==0), varList]
critEmTree[, varList]
critPln[, varList]

par(mfrow=c(2, 2))
F_PlotSlopeCrit(n*critEmMiss$q*log(critEmMiss$q), critEmMiss$JZ, critEmMiss$r, critEmMiss$d)



par(mfrow=c(3, 2), mex=.6)
F_PlotSlopeCrit(critEmMiss$pBIC, critEmMiss$JTZY, critEmMiss$r, critEmMiss$d)
F_PlotSlopeCrit(critEmMiss$pBIC, critEmMiss$JTZ, critEmMiss$r, critEmMiss$d)
F_PlotSlopeCrit(critEmMiss$pBIC, critEmMiss$JZ, critEmMiss$r, critEmMiss$d)
F_PlotSlopeCrit(critEmMiss$pBIC+critEmMiss$cardT, critEmMiss$JTZY, critEmMiss$r, critEmMiss$d)
F_PlotSlopeCrit(critEmMiss$pBIC+critEmMiss$cardT, critEmMiss$JTZ, critEmMiss$r, critEmMiss$d)
F_PlotSlopeCrit(critEmMiss$pBIC+critEmMiss$cardT, critEmMiss$JZ, critEmMiss$r, critEmMiss$d)

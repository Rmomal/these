# F.mean1 <- function(x){1 - mean(P / (x+M))}
# F.logsum0 <- function(x){-sum(log(-x+P)-log(M))}
#########################################################################
F_NegLikelihood <- function(beta.vec, log.phi, P){
   return(- sum(F_Sym2Vec(P)*(log(beta.vec)+F_Sym2Vec(log.phi))) + log(SumTree(F_Vec2Sym(beta.vec))))
}
F_NegGradient <- function(beta.vec, log.phi, P){
   M = Kirshner(F_Vec2Sym(beta.vec)*exp(log.phi))$Q
   lambda = SetLambda(P, M)
   return(- F_Sym2Vec(P)/beta.vec + F_Sym2Vec(M) + rep(lambda, length(beta.vec)))
}
#########################################################################
SetLambda <- function(P, M, eps = 1e-6){
   # F.x has to be increasing. The target value is 0
   F.x <- function(x){1 - sum(P / (x+M))}
   suite=TRUE
   if(F.x(1e-16) >0){
     F.x <- function(x){0.99 - sum(P / (x+M))}
     if(F.x(1e-16) >0){
       suite=FALSE
       x=NA
     }
   }
   if(suite){
       x.min = 1e-4;
    # print( F.x(x.min))
    # if (F.x(x.min)>0) browser()
     while(F.x(x.min)>0){x.min = x.min / 2}
     x.max = 10; while(F.x(x.max)<0){x.max = x.max * 2}
     x = (x.max+x.min)/2;
     f.min = F.x(x.min); f.max = F.x(x.max); f = F.x(x)
     # x.list = exp(seq(log(x.min), log(x.max), length.out=50))
     # plot(x.list, sapply(x.list, function(l){F.x(l)}), type='l', log='xy'); abline(h=0)
     # points(c(x.min, x, x.max), c(f.min, f, f.max), col=c(1, 2, 1))
     while(abs(x.max-x.min) > eps){
        if(f > 0){x.max = x; f.max = f}else{x.min = x; f.min = f}
        x = (x.max+x.min)/2;
        f = F.x(x)
        # points(c(x.min, x, x.max), c(f.min, f, f.max), col=c(1, 2, 1))
        # cat(x.min, x, x.max, '/', f.min, f, f.max, '\n')
     }
   }
   return(x)
}

#########################################################################
FitBetaStatic <- function(beta.init, phi, iterMax = 20, eps = 1e-4,print){
   # beta.init = beta.unif; iterMax = 1e3; eps = 1e-6; log.phi = log(phi)
   beta.tol = 1e-4
   beta.min = 1e-30
   beta.old = beta.init / sum(beta.init)
   log.phi = log(phi)
   iter = 0
   logpY = rep(0, iterMax)
   diff = 2*eps
   while ((diff > eps) & (iter < iterMax)){
      iter = iter+1
     # P = Kirshner(beta.old*phi)$P # sum(P)

      P = EdgeProba(beta.old*phi)

      beta = F_Vec2Sym(optim(F_Sym2Vec(beta.old), F_NegLikelihood, gr=F_NegGradient,
                             method='BFGS', log.phi, P)$par)
      beta[which(beta < beta.min)] = beta.min; diag(beta) = 0
      logpY[iter] = - F_NegLikelihood(F_Sym2Vec(beta), log.phi, P)
      # Test
      diff = max(abs(beta.old-beta))
      d<-ncol(P)
      beta.old = beta

     if(print) cat(" max(P) =",max(P)," sum(P)/2 =", sum(P)/2,'\n','sum(beta)= ', sum(beta),"\n Diff =",diff)

   }
   #cat(" max(P) =",max(P)," sum(P)/2 =", sum(P)/2,'\n','sum(beta)= ', sum(beta),"\n")
   cat("iter: ",iter,", diff= ",diff)
    logpY = logpY[1:iter]
   # plot(logpY)
   return(list(beta=beta, logpY=logpY))
}

FitBeta1step <- function(beta.init, phi, iterMax = 1e3, eps = 1e-6){
  # beta.init = beta.unif; iterMax = 1e3; eps = 1e-6; log.phi = log(phi)
  beta.tol = 1e-4
  beta.min = 1e-30
  beta.old = beta.init / sum(beta.init)
  log.phi = log(phi)
  iter = 0
  logpY = 0
  diff = 2*eps

    P = EdgeProba(beta.old*phi)

    beta = F_Vec2Sym(optim(F_Sym2Vec(beta.old), F_NegLikelihood, gr=F_NegGradient,
                           method='BFGS', log.phi, P)$par)
    beta[which(beta < beta.min)] = beta.min; diag(beta) = 0
    logpY= - F_NegLikelihood(F_Sym2Vec(beta), log.phi, P)
    # Test
    diff = max(abs(beta.old-beta))
    #cat( ':', min(P), max(P), sum(P),'/', sum(beta)/2, '/', logpY, diff, '\n')
    beta.old = beta


  # plot(logpY)
  return(list(beta=beta, logpY=logpY))
}


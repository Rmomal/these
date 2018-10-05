# F.mean1 <- function(x){1 - mean(P / (x+M))}
# F.logsum0 <- function(x){-sum(log(-x+P)-log(M))}
#########################################################################
F_NegLikelihood <- function(beta.vec, log.psi, P){
  return(- sum(F_Sym2Vec(P)*(log(beta.vec)+F_Sym2Vec(log.psi))) + log(SumTree(F_Vec2Sym(beta.vec))))
}
F_NegGradient <- function(beta.vec, log.psi, P){
  M = Kirshner(F_Vec2Sym(beta.vec)*exp(log.psi))$Q
  # options(error = recover)
  lambda = SetLambda(P, M)
  return(- F_Sym2Vec(P)/beta.vec + F_Sym2Vec(M) + rep(lambda, length(beta.vec)))
}
#########################################################################
SetLambda <- function(P, M, eps = 1e-6){
  # F.x has to be increasing. The target value is 0
  F.x <- function(x){
    if(x!=0){
      1 - sum(P / (x+M))
    }else{
      1 - (2*sum(P[upper.tri(P)] / M[upper.tri(M)]))
    }
  }
  suite=TRUE
  # if(F.x(1e-16) >0){
  #   F.x <- function(x){0.99 - sum(P / (x+M))}
  #   if(F.x(1e-16) >0){
  #     suite=FALSE
  #     x=NA
  #     browser()
  #     cat("ECHEC GRADIENT DESCENT : F.X(1e-16)=",F.x(1e-16)," \n")
  #   }
  # }
  if(suite){
    x.min = ifelse(F.x(1e-16) >0,-20,1e-4);
    # if (F.x(x.min)>0) browser()
    while(F.x(x.min)>0){x.min = x.min -x.min/2}
    x.max = 10; while(F.x(x.max)<0){x.max = x.max * 2}
    x = (x.max+x.min)/2;

    f.min = F.x(x.min); f.max = F.x(x.max); f = F.x(x)
    # x.list = exp(seq(log(x.min), log(x.max), length.out=50))
    # plot(x.list, sapply(x.list, function(l){F.x(l)}), type='l', log='xy'); abline(h=0)
    # points(c(x.min, x, x.max), c(f.min, f, f.max), col=c(1, 2, 1))
    while(abs(x.max-x.min) > eps){
      if(f > 0) {
        x.max = x
        f.max = f
      } else{
        x.min = x
        f.min = f
      }
      x = (x.max+x.min)/2;
      f = F.x(x)

      # points(c(x.min, x, x.max), c(f.min, f, f.max), col=c(1, 2, 1))
      # cat(x.min, x, x.max, '/', f.min, f, f.max, '\n')
    }
   # cat("lambda : ", x,"F.x(lambda) : ",f,"\n")
  }
  return(x)
}

#########################################################################
# Choice of alpha for alpha * n
F_AlphaN <- function(CorY, n, cond.tol=1e-10){
   # Grid on alpha
   alpha.grid = (1:n)/n; alpha.nb = length(alpha.grid); cond = rep(0, alpha.nb)
   for (a in 1: alpha.nb){
      psi.vec = F_Sym2Vec(-alpha.grid[a]*n*log(1 - CorY^2)/2); 
      psi.vec = psi.vec - mean(psi.vec)
      psi = F_Vec2Sym(exp(psi.vec))
      # image(psi, main=alpha.grid[a])
      lambda = svd(psi)$d
      cond[a] = (min(abs(lambda))/max(abs(lambda)))
      # plot(alpha.grid, cond, log='y')
   }
   alpha = alpha.grid[max(which(cond>cond.tol))]
   psi.vec = F_Sym2Vec(-alpha*n*log(1 - CorY^2)/2); 
   psi.vec = psi.vec - mean(psi.vec)
   psi = F_Vec2Sym(exp(psi.vec))
   return(list(psi=psi, alpha=alpha))
}

#########################################################################
FitBetaStatic <- function(beta.init, psi, maxIter, eps = 1e-6, print=FALSE){

  # beta.init = beta.unif; maxIter = 1e3; eps = 1e-6; log.psi = log(psi)
  beta.tol = 1e-4
  beta.min = 1e-30
  beta.old = beta.init / sum(beta.init)
  log.psi = log(psi)
  iter = 0
  logpY = rep(0, maxIter)
  diff = 2*eps

  while ((diff > eps) & (iter < maxIter)){
    #print(iter)
    iter = iter+1
    # P = Kirshner(beta.old*psi)$P # sum(P)
    # options(error = recover)
    P = EdgeProba(beta.old*psi)

    optimise = optim(F_Sym2Vec(beta.old), F_NegLikelihood, gr=F_NegGradient,
                     method='BFGS', log.psi, P)
    beta = F_Vec2Sym(optimise$par)
    beta[which(beta < beta.min)] = beta.min; diag(beta) = 0
    #cat("beta = ", beta)
    logpY[iter] = - optimise$value
    # Test
    beta.diff = max(abs(beta.old-beta))
    # beta.diff = max(abs(log(beta.old/beta)), na.rm=T)
    if(iter > 1){diff = max(beta.diff, abs(logpY[iter]-logpY[iter-1]))}else{diff=beta.diff}
    d<-ncol(P)
    cat(iter, logpY[iter], beta.diff, diff ,'\n')
    beta.old = beta

    #if(print) cat(" max(P) =",max(P)," sum(P)/2 =", sum(P)/2,'\n','sum(beta)= ', sum(beta),"\n Diff =",diff)

  }
  if(iter==1){
    eps=1e-6
    while ((diff > eps) & (iter < maxIter)){
      #print(iter)
      iter = iter+1
      # P = Kirshner(beta.old*psi)$P # sum(P)
      # options(error = recover)
      P = EdgeProba(beta.old*psi)
      
      beta = F_Vec2Sym(optim(F_Sym2Vec(beta.old), F_NegLikelihood, gr=F_NegGradient,
                             method='BFGS', log.psi, P)$par)
      beta[which(beta < beta.min)] = beta.min; diag(beta) = 0
      #cat("beta = ", beta)
      logpY[iter] = - F_NegLikelihood(F_Sym2Vec(beta), log.psi, P)
      # Test
      diff = max(abs(beta.old-beta))
      d<-ncol(P)
      beta.old = beta
      
      #if(print) cat(" max(P) =",max(P)," sum(P)/2 =", sum(P)/2,'\n','sum(beta)= ', sum(beta),"\n Diff =",diff)
      
    }
  }
  #cat(" max(P) =",max(P)," sum(P)/2 =", sum(P)/2,'\n','sum(beta)= ', sum(beta),"\n")
  cat("\n Fin while final iter: ",iter,", diff= ",diff)
  logpY = logpY[1:iter]
  # plot(logpY)
  P = EdgeProba(beta.old*psi)
  return(list(beta=beta, logpY=logpY, P=P,maxIter=iter))
}

FitBeta1step <- function(beta.init, psi, maxIter = 1e3, eps = 1e-6){
  # beta.init = beta.unif; maxIter = 1e3; eps = 1e-6; log.psi = log(psi)
  beta.tol = 1e-4
  beta.min = 1e-30
  beta.old = beta.init / sum(beta.init)
  log.psi = log(psi)
  iter = 0
  logpY = 0
  diff = 2*eps
  # options(error = recover)
  P = EdgeProba(beta.old*psi)

  beta = F_Vec2Sym(optim(F_Sym2Vec(beta.old), F_NegLikelihood, gr=F_NegGradient,
                         method='BFGS', log.psi, P)$par)
  beta[which(beta < beta.min)] = beta.min; diag(beta) = 0
  logpY= - F_NegLikelihood(F_Sym2Vec(beta), log.psi, P)
  # Test
  diff = max(abs(beta.old-beta))
  #cat( ':', min(P), max(P), sum(P),'/', sum(beta)/2, '/', logpY, diff, '\n')
  beta.old = beta


  # plot(logpY)
  return(list(beta=beta, logpY=logpY))
}

#########################################################################
TreeGGM<-function(CorY, n, step, print, maxIter, cond.tol=1e-10){
   p = ncol(CorY);
   alpha.psi = F_AlphaN(CorY, n, cond.tol=cond.tol)
   psi = alpha.psi$psi
   cat('alpha =', alpha.psi$alpha, '\n')
   
   beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)
   
   FitEM = switch(step,"FALSE"=FitBetaStatic(beta.init=beta.unif, psi=psi, print=print, maxIter = maxIter),
                  "TRUE"=FitBeta1step(beta.init=beta.unif, psi=psi))
   
   return(list(P=Kirshner(FitEM$beta)$P, L=FitEM$logpY, probaCond=FitEM$P, maxIter=FitEM$maxIter, alpha=alpha.psi$alpha))
}

#################################################################################
# Resampling
F_ResampleTreePLN <- function(Y, X, O, v=0.8, B=1e2, maxIter, cond.tol=1e-10){
   # v = 0.8; B = 1e2
   # Resampling edge probability
   # Y, X, O: same as for PLN
   # v = (subsample size) / (total sample size)
   # B = nb resamples
   # Out = Pmat = B x p(p-1)/2 matrix with edge probability for each resample
   # Y = Data$count; X = matrix(1, n, 1); O = matrix(0, n, p)
   n = nrow(Y); p = ncol(Y); P = p*(p-1)/2; V = round(v*n); Pmat = matrix(0, B, P); 
   alpha = rep(0, B)
   for (b in 1:B){
      cat('\n', b, '')
      sample = sample(1:n, V, replace=F)
      Y.sample = Y[sample, ]; X.sample = X[sample, ]; O.sample = O[sample, ]
      # if(ncol(Y)>1){Y.sample = Y[sample, ]}else{Y.sample = matrix(Y[sample], V, 1)}; 
      # if(ncol(X)>1){X.sample = X[sample, ]}else{X.sample = matrix(X[sample], V, 1)}; 
      # if(ncol(O)>1){O.sample = O[sample, ]}else{O.sample = matrix(O[sample], V, 1)}; 
      PLN.sample = PLN(Y.sample ~ -1 + X.sample + offset(O.sample))
      Sigma.sample = PLN.sample$model_par$Sigma
      TGGM = TreeGGM(cov2cor(Sigma.sample), n=length(sample), step='FALSE', maxIter=maxIter, cond.tol=cond.tol)
      Pmat[b, ] = F_Sym2Vec(TGGM$probaCond)
      alpha[B] = TGGM$alpha
   }
   return(Pmat, alpha)
}

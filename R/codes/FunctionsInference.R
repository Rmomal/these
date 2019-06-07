# F.mean1 <- function(x){1 - mean(P / (x+M))}
# F.logsum0 <- function(x){-sum(log(-x+P)-log(M))}
#########################################################################
F_NegLikelihood <- function(beta.vec, log.psi, P){
  return(- sum(F_Sym2Vec(P)*(log(beta.vec)+F_Sym2Vec(log.psi))) + 
           log(SumTree(F_Vec2Sym(beta.vec))))
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
  alpha.grid = (1:n)/n; alpha.nb = length(alpha.grid); 
  # cond = rep(0, alpha.nb)
  # for (a in 1:alpha.nb){
  #   psi.vec = F_Sym2Vec(-alpha.grid[a]*n*log(1 - CorY^2)/2);
  #   psi.vec = psi.vec - mean(psi.vec)
  #   psi = F_Vec2Sym(exp(psi.vec))
  #   # image(psi, main=alpha.grid[a])
  #   lambda = svd(psi)$d
  #   cond[a] = (min(abs(lambda))/max(abs(lambda)))
  #   # plot(alpha.grid, cond, log='y')
  # }
  # alpha = alpha.grid[max(which(cond>cond.tol))]
  cond = Inf; a = 0
  while(cond > cond.tol && a<length(alpha.grid)){
    a = a+1
    psi.vec = F_Sym2Vec(-alpha.grid[a]*n*log(1 - CorY^2)/2);
    psi.vec = psi.vec - mean(psi.vec)
    psi = F_Vec2Sym(exp(psi.vec))
    lambda = svd(psi)$d
    cond = min(abs(lambda))/max(abs(lambda))
  }
  alpha = alpha.grid[a-1]
  psi.vec = F_Sym2Vec(-alpha*n*log(1 - CorY^2)/2); 
  psi.vec = psi.vec - mean(psi.vec)
  psi = F_Vec2Sym(exp(psi.vec))
  return(list(psi=psi, alpha=alpha))
}

#########################################################################
FitBetaStatic <- function(beta.init, psi, maxIter, eps1 = 1e-6,eps2=1e-4, print=FALSE,
                          indices_save=indices){
  beta.tol = 1e-4;  beta.min = 1e-30; beta.old = beta.init / sum(beta.init)
  log.psi = log(psi)
  iter = 0;  logpY = rep(0, maxIter);  times = c(); listP=vector("list", length(indices_save))
  beta.diff = diff.loglik=2*eps2
  
  T1<-Sys.time()
  while (((beta.diff > eps1) || (diff.loglik>eps2) ) && iter < maxIter ){
    iter = iter+1
    P = EdgeProba(beta.old*psi)
    beta = F_Vec2Sym(optim(F_Sym2Vec(beta.old), F_NegLikelihood, gr=F_NegGradient,method='BFGS', log.psi, P)$par)
    beta[which(beta< beta.min)] = beta.min
    diag(beta) = 0
    logpY[iter] = -F_NegLikelihood(F_Sym2Vec(beta),log.psi,P) 
    
    # Test
    beta.diff = max(abs(beta.old-beta))
    beta.old = beta
    
    if(iter > 1){diff.loglik =  abs(logpY[iter]-logpY[iter-1])}else{diff.loglik=1}
    # if(iter %in% indices_save){
    #   listP[[which(indices_save==iter)]]<-P
    #   times<-c(times,difftime(Sys.time(),T1))
    # }
    
  }
  times<-c(times,difftime(Sys.time(),T1))
  cat("\n Fin while final iter: ",iter,", diff= ",diff.loglik ,'\n')
  logpY = logpY[1:iter]
  # if(iter<maxIter){
  #   listP<-lapply(listP, function(x){
  #     if(is.null(x)){
  #       x<-P  
  #       times<<-c(times,difftime(Sys.time(),T1))
  #     }else{
  #       x<-x
  #     }
  #   })
  # }
  # plot(logpY)
  P = EdgeProba(beta.old*psi)
  return(list(beta=beta, logpY=logpY, listP=listP,P=P,maxIter=iter, times=times))
}

FitBeta1step <- function(beta.init, psi, maxIter = 1e3, eps = 1e-6){
  beta.min = 1e-30
  beta.old = beta.init / sum(beta.init)
  log.psi = log(psi)
  iter = 0
  logpY = 0
  
  P = EdgeProba(beta.old*psi)
  beta = F_Vec2Sym(optim(F_Sym2Vec(beta.old), F_NegLikelihood, gr=F_NegGradient,
                         method='BFGS', log.psi, P)$par)
  beta[which(beta < beta.min)] = beta.min; diag(beta) = 0
  logpY= - F_NegLikelihood(F_Sym2Vec(beta), log.psi, P)
  beta.old = beta
  print("iter=1 !")
  return(list(beta=beta, logpY=logpY, P=P))
}

#########################################################################
TreeGGM<-function(CorY, n, print, maxIter, cond.tol=1e-10,indices=0){
  p = ncol(CorY);
  alpha.psi = F_AlphaN(CorY, n, cond.tol=cond.tol)
  psi = alpha.psi$psi
  cat('alpha =', alpha.psi$alpha, '\n')
  
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)
  
  FitEM = FitBetaStatic(beta.init=beta.unif, psi=psi, print=print, maxIter = maxIter,
                        indices_save = indices)
  
  return(list( L=FitEM$logpY, probaCond=FitEM$P, maxIter=FitEM$maxIter, alpha=alpha.psi$alpha, 
               listP=FitEM$listP, times=FitEM$times))
}

#################################################################################
# Resampling
F_ResampleTreePLN <- function(Y, X, O, v=0.8, B=1e2, maxIter, cond.tol=1e-14,cores=3){
  # v = 0.8; B = 1e2
  # Resampling edge probability
  # Y, X, O: same as for PLN
  # v = (subsample size) / (total sample size)
  # B = nb resamples
  # Out = Pmat = B x p(p-1)/2 matrix with edge probability for each resample
  # Y = Data$count; X = matrix(1, n, 1); O = matrix(0, n, p)
  n = nrow(Y); p = ncol(Y); P = p*(p-1)/2; V = round(v*n); Pmat = matrix(0, B, P); 
  # time<-rep(0,B)
  
  obj<-mclapply(1:B,function(b){
    #  cat('\n', b, '')
    set.seed(b)
    sample = sample(1:n, V, replace = F)
    Y.sample = Y[sample,]
    X.sample = X[sample,]
    O.sample = O[sample,]
    
    PLN.sample = PLN(Y.sample ~ -1 + X.sample + offset(O.sample))
    Sigma.sample = PLN.sample$model_par$Sigma
    inf<-TreeGGM(cov2cor(Sigma.sample), n=n, maxIter=maxIter, 
                 cond.tol=cond.tol)[c("times")]
    
    return(inf)
  },mc.cores=cores)
  #browser()
  #  Pmat<-do.call(rbind,lapply(obj,function(x){F_Sym2Vec(x$probaCond)}))
  # ListAlpha = do.call(c,lapply(obj,function(x){x$alpha}))
  #  summaryiter = do.call(c,lapply(obj,function(x){x$maxIter}))
  times<-do.call(c,lapply(obj,function(x){x$time}))
  #list(Pmat,maxIter=summaryiter,times=times)
  return(times)
}


save_iterations <- function(Y, X, O,maxit ,B=1e2, cond.tol=1e-14){
  
  res<-mclapply(1:maxit,function(x){
    cat('\n', x, ' max iteration')
    F_ResampleTreePLN_oak(maxIter=x,Y=Y, X=X, O=O, B=B, cond.tol=cond.tol)
  },mc.cores=1 )
  
  return(res)
}


F_ResampleTreePLN_oak<-function(Y, X, O, maxIter,v=0.8, B=1e2, cond.tol=1e-14, indices){
  n = nrow(Y); P = p*(p-1)/2; V = round(v*n); Pmat = matrix(0, B, P); 
  
  obj<-lapply(1:B,function(b){
    #  cat('\n', b, '')
    set.seed(b)
    sample = sample(1:n, V, replace = F)
    Y.sample = Y[sample,]
    X.sample = X[sample,]
    O.sample = O[sample,]
    
    PLN.sample = PLN(Y.sample ~ -1 + X.sample + offset(O.sample))
    Sigma.sample = PLN.sample$model_par$Sigma
    inf<-TreeGGM(cov2cor(Sigma.sample), n=n, maxIter=maxIter, cond.tol=cond.tol,
                 indices=indices)[c("listP","alpha","L","maxIter","times")]
    
    return(inf)
  })
  #browser()
  
  listPmat<-lapply(seq_along(indices),function(i){
    do.call(rbind,lapply(obj,function(x){F_Sym2Vec(x$listP[[i]])}))
  })
  times<-do.call(c,lapply(obj,function(x){x$times}))
  L<-do.call(c,lapply(obj,function(x){x$L}))
  alpha<-do.call(c,lapply(obj,function(x){x$alpha}))
  maxIter<-do.call(c,lapply(obj,function(x){x$maxIter}))
  
  return(list(listPmat,times,L,alpha,maxIter))
}







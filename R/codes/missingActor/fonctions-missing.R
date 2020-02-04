#############
# initialization
#Find the clique created by the missing actor

init.mclust<-function(S,nb.missing=1, n.noise=50,plot=TRUE, title="",trueClique=NULL){
  # browser()
  X<-data.frame(t(t(eigen(S)$vectors[,1:2])*sqrt(eigen(S)$values[1:2])))
  
  p=nrow(X)+nb.missing
  b=apply(X, 2, range)
  poissonNoise<-apply(b, 2, function(x,n=n.noise){
    runif(n,min=x[1]-0.1, max=x[2]+0.1)
  })
  data=rbind(X, poissonNoise)
  noiseInit<-sample(c(T,F), size=nrow(X)+n.noise, replace=T, prob=c(3, 1))
  
  datapolar=cart2pol(x=data[,1],y=data[,2])[,1:2]
  datapolar=datapolar %>% mutate(theta2=ifelse(theta>pi,theta-pi,theta)) %>% dplyr::select(r,theta2)
  newdata=pol2cart(datapolar$r,datapolar$theta2)[,1:2]
  
  clust=Mclust(data=newdata,initialization = list(noise=noiseInit),G=nb.missing)
  groups<-mclust::map(clust$z)
  res<-which(groups==1)[which(groups==1)<=nrow(X)] 
  if(plot){
    if(!is.null(trueClique)){
      #False positives rate
      N=setdiff(1:p,trueClique)
      FP=sum(res%in%N)/length(N)
      # False negatives rate
      FN=sum(setdiff(1:p,res)%in%trueClique)/length(trueClique)
      title=paste0(title,":"," FN=", round(FN,2),",FP=",round(FP,2))
    }
    
    g= ggplot(X,aes(X1,X2, label=rownames(X),color=rownames(X)%in%res))+geom_point(size=0.1)+
      theme_light()+geom_text()+labs(x="eig vect 1",y="eig vect 2", title=title)+
      guides(color=FALSE)+scale_color_brewer(palette="Dark2")+
      geom_hline(yintercept=0, color="gray50")+geom_vline(xintercept=0, color="gray50")
    print(g)
  }
  return(res)
}

#################
# For OPTIM

# Min for VE
argminKL <- function(gamma, Pg, M,S,omega,W,p,lambda ){
  
  r=ncol(omega)-p
  O = 1:p
  H=(p+1):(p+r)
  omegaH=omega[H,H]
   
  EhZoZo = t(M[,O])%*%M[,O]+ diag(colSums(S[,O]))
  
  KL <- -n*0.5*sum(log(diag(omega)))+0.5*omegaH*( t(M[,H])%*%M[,H]+sum(S[,H]) ) +   
    sum(diag(Pg[O,H]*omega[O,H] %*% (t(M[,H])%*%M[,O]) )) + 
    0.5*sum(diag( (Pg[O,O]*omega[O,O]) %*% EhZoZo)) -
    n*0.5*2*sum(F_Sym2Vec(Pg)*log(F_Sym2Vec(CorOmegaMatrix(omega))))+
    2*sum(F_Sym2Vec(Pg)*(gamma-log(F_Sym2Vec(W))))-log(SumTree(F_Vec2Sym(exp(gamma))))+
    log(SumTree(W))-lambda*(sum(exp(gamma)) - 0.5) -0.5*sum(log(S[,H]))
  
 
  return(KL)
}

Grad_KL_Wg <- function(gamma, Cg, Pg, M,S,omega,W){
  Mei = Meila(F_Vec2Sym(exp(gamma)))
  lambda = SetLambda(Pg, Mei)
  return( F_Sym2Vec(Pg) - exp(gamma)*(F_Sym2Vec(Mei) + lambda)) # inversion du signe du lambda
}

# Max for Mstep

argmaxJ<-function(gamma,Pg,omega,sigmaTilde,lambda,n, trim=TRUE, verbatim=TRUE){
  
  
  
  W=F_Vec2Sym(exp(gamma))
  
  #browser()
  logphi=log(CorOmegaMatrix(omega))
  
  maxJ <- sum(Pg*F_Vec2Sym(F_Sym2Vec(log(W)) + n*0.5*F_Sym2Vec(logphi))) + 
    n*0.5*sum(log(diag(omega))) - log(SumTree(W)) - 
    0.5*sum( Pg * omega * sigmaTilde * n) - n*0.5*sum(diag(omega*sigmaTilde)) 
  + lambda*(sum(W)-0.5)
  
  return(maxJ)
}



wrap_optimDiag <- function(omega_ii, i) {
  return(optimDiag(i, omega_ii, omega, SigmaTilde, Pg))
}

optimDiag <- function(omega_ii,i,  omega, SigmaTilde, Pg) {
  # browser()
  q <- nrow(omega)
  phi = CorOmegaMatrix(omega)
  newphi=(1- phi)/phi
  diag(newphi)=0
  quotient = diag(Pg %*% newphi)[i]
  
  grad=(1/omega_ii)*(sum(quotient)+1) - SigmaTilde[i, i]
  return( grad)
}

dichotomie <- function(a, b, f, epsilon){
  # ---------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   find the zero of monotonous function f by binary search.
  # INPUT
  #   a        : minimum argument of f
  #   b        : maximum argument of f
  #   epsilon  : tolerance (length of last interval)
  # OUTPUT
  #   c        : approximate zero of f with epsilon tolerance
  # ---------------------------------------------------------------------------------------------------------
  
  min <- a
  max <- b
  sgn <- sign(f(max)-f(min))
  c <- (max+min)/2
  # browser()
  print(f(c))
  while(abs(f(c))>1e-5){
    c <- (max+min)/2
    if(sgn*f(c)>0){
      max <- c
    } else{
      min <- c
    }
  }
  return(c)
}

computeOffDiag<-function(omegaDiag,SigmaTilde){
  q=length(omegaDiag)
  omega=matrix(0,q,q)
  sapply(1:(q-1),
         function(j){
           sapply((j+1):q,
                  function(k){
                    omega[k, j] <<- (1 - sqrt( 1+4*SigmaTilde[j,k]^2*omegaDiag[j]*omegaDiag[k]))/(2*SigmaTilde[j,k])
                    omega[j, k] <<- omega[k, j]
                  }
           )
         }
  )
  diag(omega)=omegaDiag
  return(omega)
}

EdgeProba <- function(W, verbatim=FALSE, p.min=1e-10){
  it=-1
  Wcum = SumTree(W)
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
  # while(!is.finite(Wcum)){
  #   #handles numerical issues with matrix tree theorem
  #   it=it+1
  #   borne=30-it
  #   if(verbatim) message(cat("W corrected, bound=",borne))
  #   
  #   W.log=log(F_Sym2Vec(W))
  #   W.center=W.log-mean(W.log)
  #   W.center[which(W.center<(-borne))]=-borne
  #   W=F_Vec2Sym(exp(W.center))
  #   Wcum = SumTree(W)
  # }
  
  p = nrow(W); P = matrix(0, p, p)
  #core of computation
  sapply(1:(p-1),
         function(j){
           sapply((j+1):p,
                  function(k){
                    W_jk = W; W_jk[j, k] = W_jk[k, j] = 0 #kills kj edge in W_kj
                    P[k, j] <<- 1 - SumTree(W_jk) / Wcum
                    P[j, k] <<- P[k, j]
                  }
           )
         }
  )
  P[which(P<p.min)]=p.min
  diag(P)=0
  return(P)
}

SetLambda2<- function(P, M,n, eps = 1e-6){
  # F.x has to be increasing. The target value is 0
  #browser()
  p=ncol(P)
  F.x <- function(x){
    if(x!=0){
      sum(log(P / (x+M)))/(p*(p-1))
    }else{
      (2*sum(log(P[upper.tri(P)] / M[upper.tri(M)])))/(p*(p-1))
    }
  }
  x.min = ifelse(F.x(0) >0,-20,1e-4);
  while(F.x(x.min)>0){x.min = x.min -x.min/2}
  x.max = 10
  while(F.x(x.max)<0){x.max = x.max * 2}
  x = (x.max+x.min)/2
  f.min = F.x(x.min)
  f.max = F.x(x.max)
  f = F.x(x)
  
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
  }
  if(abs( 1 - sum(P / (x+M)))>2) browser()
  return(x)
}

#####################
# For hidden proba

HiddenEdgeProba <- function(W,r=1, verbatim=FALSE){
  #computes the probability for two links to a hidden covariates to NOT be there
  #Coded for 1 hidden covariate
  it=-1
  Wcum = SumTree(W)
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
  # while(!is.finite(Wcum)){
  #   #handles numerical issues with matrix tree theorem
  #   it=it+1
  #   borne=30-it
  #   if(verbatim)cat("W corrected, bound=",borne)
  #   
  #   W.log=log(F_Sym2Vec(W))
  #   W.center=W.log-mean(W.log)
  #   W.center[which(W.center<(-borne))]=-borne
  #   W=F_Vec2Sym(exp(W.center))
  #   Wcum = SumTree(W)
  # }
  
  q = nrow(W); P = matrix(0, q, q)
  p=q-1 #one h
  h=q # the hidden covariate is stored in last row/column
  #core of computation
  sapply(1:p-1,
         function(j){
           sapply((j+1):p, #visits all combinations of observed nodes
                  function(k){
                    W_jk = W 
                    W_jk[h, k] = W_jk[k, h] = 0 #kills kh edge in W_kj
                    W_jk[h, j] = W_jk[j, h] = 0 #kills jh edge in W_kj
                    P[k, j] <<- 1 - SumTree(W_jk) / Wcum
                    P[j, k] <<- P[k, j]
                  }
           )
         }
  )
  P[which(P<1e-10)]=1e-10 # why not 0 ?
  diag(P)=0
  return(P)
}

CgMatrix<-function(Pg,Pghkl,omega,p){ # code for 1 hidden covariate. Otherwise, must be a sum on h
  h=ncol(omega)
  omegaH=omega[h,h]
  Cg = matrix(1/omegaH,p,p)
  
  sapply(1:p-1,
         function(j){
           sapply((j+1):p, #visits all combinations of observed nodes
                  function(k){
                    Cg[k, j] <<- omega[h,k]*omega[h,j] * (Pg[h,k]+Pg[h,j]+Pghkl[k,j] -1)
                    Cg[j, k] <<- Cg[k, j]
                  }
           )
         })
  return(Cg)
}
CorOmegaMatrix<-function(omega){ # code for 1 hidden covariate. Otherwise, must be a sum on h
  q=ncol(omega)
  CorOmega = matrix(0,q,q)
  
  sapply(1:q-1,
         function(j){
           sapply((j+1):q, #visits all combinations of observed nodes
                  function(k){
                    CorOmega[k, j] <<- 1-(omega[k,j]^2/(omega[k,k]*omega[j,j]))
                    CorOmega[j, k] <<- CorOmega[k, j]
                  }
           )
         })
  return(CorOmega)
}

computeWtree<-function(omega, W, Wg, MH, MO, SO,phi, alpha=0,trim=TRUE){
  n=nrow(MH) ; r=ncol(MH) ; p=ncol(MO)
  O = 1:p ; H=(p+1):(p+r)
  # browser()
  sigmaO<-(t(MO)%*%MO + diag(colSums(SO)))/n
  
  A<-log(Wg*phi^(-n*0.5)/W)
  #  browser()
  B<-(0.5*omega[O,O]*t(MO)%*%MO*alpha)
  C<-(omega[O,H]*t(MH)%*%MO)# - 2*n*omega[H,O]*diag(omega[O,O]%*%sigmaO)/omega[H,H])
  # browser()
  logWtree<-matrix(0,p+r,p+r)
  logWtree[O,O]<-A[O,O]+B
  logWtree[O,H]<-A[O,H]+C
  logWtree[H,O]<-t(logWtree[O,H])
  logWtree[H,H]<- A[H,H]
  diag(logWtree) = 0
  if(trim){
    
    logWtree=logWtree-mean(logWtree)
    logWtree[which(logWtree<(-30))]=-30
    
  }
  
  return(logWtree)
}
adjustW<-function(W, mult=50){
  #    
  gamma=log(W)
  while(!is.finite(SumTree(exp(gamma)))){
    gamma=gamma-mean(gamma)
    gamma[which(gamma<(-30))]=-30
    browser()
    q95=quantile(F_Sym2Vec(gamma),probs = 0.95)
    if(max(gamma)>3*q95) gamma[which(gamma>3*q95)]=gamma[which(gamma>3*q95)]-(max(gamma)-3*q95)
  }
  # while(abs(mean(gamma))>0.2){
  #   gamma=gamma-mean(gamma)
  #   gamma[which(gamma<(-30))]=-30
  # }
  return(exp(gamma))
}
trimW<-function(W, trim=TRUE,verbatim=FALSE){
  #browser()
  gamma=log(W)
  if(trim){
    gamma=gamma-mean(gamma)
    gamma[which(gamma<(-30))]=-30
  } 
  Wcum = SumTree(exp(gamma))
  W=exp(gamma)
  it=0
  while(!is.finite(Wcum)){
    #handles numerical issues with matrix tree theorem
    it=it+1
    borne=30-it
    if(verbatim) message(cat("W corrected, bound=",borne))
    
    W.log=log(F_Sym2Vec(W))
    W.center=W.log-mean(W.log)
    W.center[which(W.center<(-borne))]=-borne
    W=F_Vec2Sym(exp(W.center))
    Wcum = SumTree(W)
  }
  
  return(W)
}

########################################
# GENERAL

generator_PLN<-function(Sigma,covariates=NULL, n=50){
  p<-ncol(Sigma)
  if(!is.null(covariates)){
    n<-nrow(covariates)
    c<-ncol(covariates)
    
    string<-paste0("~", paste(colnames(covariates), collapse=" + "))
    formula<-as.formula(string)
    m<- model.matrix(formula,covariates)[,-1]
    
    mc<-ncol(m)
    beta<-matrix(runif(p*mc),mc,p)
    prod=m %*% beta
  }else{
    prod=0
  }
  
  Z<- rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
  while(sum((Z+prod)>21)>0){
    #  browser()
    badlines=which((Z+prod)>21, arr.ind = TRUE)[,1]
    Zbad<-Z[badlines,]
    Zgood=Z[-badlines,]
    Znew= rmvnorm(nrow(Zbad), rep(0,nrow(Sigma)), Sigma)
    Z=rbind(Zgood,Znew)
  }
  Y = matrix(rpois(n*p, exp(Z+prod )), n, p)
  return(Y)
}

B_alpha <- function(B, n, cond.tol=1e-10){
  #TODO
  # Grid on alpha
  # alpha.grid = (1:n)/n; alpha.nb = length(alpha.grid);
  # cond = Inf; a = 0
  # while(cond > cond.tol && a<length(alpha.grid)){
  #   a = a+1
  #   psi.vec = F_Sym2Vec(-alpha.grid[a]*n*log(1 - CorY^2)/2);
  #   psi.vec = psi.vec - mean(psi.vec)
  #   psi = F_Vec2Sym(exp(psi.vec))
  #   lambda = svd(psi)$d
  #   cond = min(abs(lambda))/max(abs(lambda))
  # }
  # alpha = alpha.grid[a-1]
  # psi.vec = F_Sym2Vec(-alpha*n*log(1 - CorY^2)/2);
  # psi.vec = psi.vec - mean(psi.vec)
  # psi = F_Vec2Sym(exp(psi.vec))
  # return(list(psi=psi, alpha=alpha))
}

auc<-function(pred,obs){
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  ROC_auc <- round(performance(prediction,"auc")@y.values[[1]],digits=2)
  return(ROC_auc)
}

accmax<-function(obs,label){
  pred <- prediction(obs,label)
  perf <- performance(pred, "acc")
  
  accmax= max(perf@y.values[[1]])
  return(accmax)
}

ggimage<-function(data){
  melted_data <- melt(data)
  ggplot(melted_data, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +guides(fill=FALSE)+ theme(plot.title = element_text(size=10, hjust=0.5))
}


accppvtpr<-function(probs,ome,h, seuil=0.5){
  Acc=round(mean(1*(probs>seuil)==ome),2) #(TP+TN)/(TP+FP+TN+FN)
  AccH=round(mean(1*(probs[h,]>seuil)==ome[h,]),2)
  AccO=round(mean(1*(probs[-h,-h]>seuil)==ome[-h,-h]),2)
  PPV=round(sum((ome==1)*(probs>seuil))/(sum((ome==1)*(probs>seuil))+
                                         sum((ome==0)*(probs>seuil))),2)#TP/(TP+FP)
  PPVH=round(sum((ome[h,]==1)*(probs[h,]>seuil))/(sum((ome[h,]==1)*(probs[h,]>seuil))+
                                                  sum((ome[h,]==0)*(probs[h,]>seuil))),2)
  PPVO=round(sum((ome[-h,-h]==1)*(probs[-h,-h]>seuil))/(sum((ome[-h,-h]==1)*(probs[-h,-h]>seuil))+
                                                        sum((ome[-h,-h]==0)*(probs[-h,-h]>seuil))),2)
  
  TPR=round(sum((ome==1)*(probs>seuil))/sum(ome==1), 2)
  TPRH=round(sum((ome[h,]==1)*(probs[h,]>seuil))/sum(ome[h,]==1), 2)
  TPRO=round(sum((ome[-h,-h]==1)*(probs[-h,-h]>seuil))/sum(ome[-h,-h]==1), 2)
  
  return(c(Acc,AccH,AccO,PPV,PPVH,PPVO,TPR,TPRH,TPRO))
}

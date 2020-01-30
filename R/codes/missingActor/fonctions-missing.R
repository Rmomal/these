# fonctions for generate data with missing actor
gener.data<-function(p=20,n=200,counts=FALSE, type="starerdos", v=0.001,plot=TRUE,prop=10/p){
  #  data=data_from_scratch(type = type,p = p,n = n,signed = FALSE,prob = 5/p,v = v)
  
  
  #hidden=which.max(diag(omega)) #on cache le noeud de plus fort degré
  star.graph <- graphModel$new(type = type,size=p, p.or.m = prop)
  star.model <- GGMmodel$new(graph=star.graph,nb.missing.var= 1)
  adjmat=star.model$getAdjmat()
  hidden=star.model$missing.var.list
  
  #omega=covarianceFromGraph(adjmat,missing.var.list = hidden,
  #                         method="raph",prop.positive.cor = 0 ,
  #                         alpha.hidden = 1, alpha.observed = 1.5)
  omega=star.model$K
  # omega=data$omega
  sigma=star.model$Sigma
  Kh  <- omega[hidden,hidden]
  Ko  <- omega[-hidden,-hidden]
  Koh <- omega[-hidden,hidden]
  Km  <- Ko - Koh %*%solve(Kh)%*% t(Koh)
  cat(paste0("\ndistance Km-Ko = ",mean((Km-Ko)^2)))
  sigma_obs=solve(Km)
  #  K  <- params
  # les voisins du noeuds qu'on cache
  trueClique=which(omega[hidden,-hidden]!=0)
  
  counts_obs=generator_PLN(sigma_obs,covariates = NULL,n=n)
  # counts_obs=counts[,-hidden]
  sigma_obs=PLN(counts_obs~1)$model_par$Sigma # ajout d'un bruit d'estimation
  # les résultats $get... de LITree ne suffisent pas pour PLN
  
  if(plot){
    G=draw_network(adjmat,groupes=1*(diag(omega)==diag(omega)[hidden]), 
                   layout="nicely",curv=0,nb=2,pal="black",nodes_label = 1:p)$G
    print(G)
  }
  #sigma_obs est ici l'estimation de SigmaO, inverse de Omega marginal
  return(list(sigmaO=sigma_obs, omegaO=omega_obs,omega=omega,clique=trueClique, hidden=hidden))
}


init.mclust<-function(S,nb.missing=1, n.noise=50,plot=TRUE, title="",trueClique=NULL){
  X<-data.frame(t(t(eigen(S)$vectors[,1:2])*sqrt(eigen(S)$values[1:2])))
  
  p=nrow(X)+nb.missing
  b=apply(X, 2, range)
  poissonNoise<-apply(b, 2, function(x,n=n.noise){
    runif(n,min=x[1]-0.1, max=x[2]+0.1)
  })
  data=rbind(X, poissonNoise)
  noiseInit<-sample(c(T,F), size=nrow(X)+n.noise, replace=T, prob=c(3, 1))
  
  datapolar=cart2pol(x=data[,1],y=data[,2])[,1:2]
  datapolar=datapolar %>% mutate(theta2=ifelse(theta>pi,theta-pi,theta)) %>% select(r,theta2)
  newdata=pol2cart(datapolar$r,datapolar$theta2)[,1:2]
  
  clust=Mclust(data=newdata,initialization = list(noise=noiseInit),G=nb.missing)
  groups<-map(clust$z)
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
argminKL <- function(gamma, Cg, Pg, M,S,omega,W,lambda,trim=TRUE){
  p=ncol(Cg)
  r=ncol(omega)-p
  O = 1:p
  H=(p+1):(p+r)
  omegaH=omega[H,H]
  # gamma = log(beta), here for Wg
  if(trim){
    gamma=gamma-mean(gamma)
    gamma[which(gamma<(-30))]=-30
  }
  #browser()
  #Mei=Meila( F_Vec2Sym(exp(gamma)))
  # lambda = SetLambda(Pg, Mei )
  EhZoZo = t(M[,O])%*%M[,O]+ diag(colSums(S[,O]))
  
  # terme 1 : T1=-Egh[log p(ZH | ZO)]
  T1 <- -(n*0.5**log(omegaH)-0.5*omegaH*( t(M[,H])%*%M[,H]+sum(S[,H]) ) + 0.5*sum(diag(Cg%*%EhZoZo)) -
            sum(diag(Pg[O,H]*omega[O,H] %*% (t(M[,H])%*%M[,O]))) )
  
  # terme 2 : T2=-Eg[log p(T)] - entr(g(T))
  # browser()
  T2 <- -(2*sum(F_Sym2Vec(Pg)*(log(F_Sym2Vec(W)) - gamma) - log(SumTree(W)) + #2* parce que les vecteurs sont la moitié des matrices
                  log(SumTree(F_Vec2Sym(exp(gamma)))) + # - car minimisation
                  lambda*(sum(exp(gamma)) - 0.5)))  # 0.5 pour somme des beta /2
  
  # terme 3 : T3=Eh[log h(ZH | ZO)]
  T3 <- -0.5*sum(log(S[,H])) - n*0.5*(1+log(2*pi))
  
  # terme 4 : T4=EghO* [log p(ZO | T)]
  T4 <- -(n*0.5*(2*sum(F_Sym2Vec(Pg)*log(F_Sym2Vec(CorOmegaMatrix(omega)))) + sum(log(diag(omega[O,O]))))-
            0.5*sum(diag( (Pg[O,O]*omega[O,O] - Cg) %*% EhZoZo)) )
  
  KL = T1+T2+T3+T4
  if(is.nan(KL)){
    cat(max(gamma),": higher bound ")
    gamma[which(gamma>(30))]=30
    Mei = Meila(F_Vec2Sym(exp(gamma)))
    lambda = SetLambda(Pg, Mei )
    T2 <- -(2*sum(F_Sym2Vec(Pg)*(log(F_Sym2Vec(W)) - gamma) - log(SumTree(W)) + #2* parce que les vecteurs sont la moitié des matrices
                    log(SumTree(F_Vec2Sym(exp(gamma)))) + # - car minimisation
                    lambda*(sum(exp(gamma)) - 0.5)))  # 0.5 pour somme des beta /2
    
    KL = T1+T2+T3+T4
    if(is.nan(KL))   cat("\nbeta tilde optimization failed\n")
  }
  return(KL)
}

Grad_KL_Wg <- function(gamma, Cg, Pg, M,S,omega,W){
  Mei = Meila(F_Vec2Sym(exp(gamma)))
  lambda = SetLambda(Pg, Mei)
  return( F_Sym2Vec(Pg) - exp(gamma)*(F_Sym2Vec(Mei) + lambda)) # inversion du signe du lambda
}

# Max for Mstep

argmaxJ<-function(gamma,Pg,omega,sigmaTilde,n){
  W=F_Vec2Sym(log(gamma))
  Mei = Meila(F_Vec2Sym(exp(gamma)))
  lambda = SetLambda(Pg, Mei )
  
  maxJ <- sum(Pg*(log(W) + n*0.5*log(CorOmegaMatrix(omega)))) + n*0.5*sum(log(diag(omega))) - log(SumTree(W)) - 
    0.5*sum( Pg * omega * sigmaTilde * n) + lambda(sum(W)-0.5)
  return(maxJ)
}



wrap_optimDiag <- function(omega_ii, i) {
  return(optimDiag(i, omega_ii, omega, SigmaTilde, Pg))
}

optimDiag <- function(i, omega_ii, omega, SigmaTilde, Pg) {
  
  q <- nrow(omega)
  phi = CorOmegaMatrix(omega)
  quotient = diag(Pg %*% ((1- phi)/phi))[i]
  
  grad=(1/omega_ii)*(sum(quotient)+1) - Sigma[i, i]
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

EdgeProba <- function(W, verbatim=FALSE){
  it=-1
  Wcum = SumTree(W)
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
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
  P[which(P<1e-10)]=1e-10
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
  while(!is.finite(Wcum)){
    #handles numerical issues with matrix tree theorem
    it=it+1
    borne=30-it
    if(verbatim)cat("W corrected, bound=",borne)
    
    W.log=log(F_Sym2Vec(W))
    W.center=W.log-mean(W.log)
    W.center[which(W.center<(-borne))]=-borne
    W=F_Vec2Sym(exp(W.center))
    Wcum = SumTree(W)
  }
  
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

computeWtree<-function(omega, W, Wg, MH, MO, SO){
  n=nrow(MH) ; r=ncol(MH) ; p=ncol(MO)
  O = 1:p ; H=(p+1):(p+r)
  # browser()
  phi <- CorOmegaMatrix(omega)
  
  sigmaO<-(t(MO)%*%MO + diag(colSums(SO)))/n
  
  A<-log(Wg*phi^(-n*0.5)/W)
  B<-(n*0.5*omega[O,O]*sigmaO)
  C<-(omega[O,H]*t(MH)%*%MO - 2*n*omega[H,O]*diag(omega[O,O]%*%sigmaO)/omega[H,H])
  
  logWtree<-matrix(0,p+r,p+r)
  logWtree[O,O]<-A[O,O]+B
  logWtree[O,H]<-A[O,H]+C
  logWtree[H,O]<-t(logWtree[O,H])
  logWtree[H,H]<- A[H,H]
  return(logWtree)
}











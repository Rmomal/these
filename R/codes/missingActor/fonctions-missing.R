#############
# initialization
initEM <- function(Sigma = NULL,n=1e6,cst=1.1,
                   #                   S = NULL,
                   pca=TRUE,cliqueList) {
  # -----------------------------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   Initialize EM by selecting a possible clique and taking the principal component of this clique
  #   as initial value for the hidden variable
  # INPUT
  #   X      : data matrix (nxp)
  #   clique : vector containing the indices of the nodes in the clique. Default is NULL and the clique
  #            is determined with hierarchical clustering.
  # OUTPUT
  #   Sigma0    : initial value of complete covariance matrix ((p+1)x(p+1) matrix)
  #   K0        : initial value of complete precision matrix ((p+1)x(p+1) matrix)
  #   clique    : vector containing the indices of the nodes in the clique
  # -----------------------------------------------------------------------------------------------------------------------------
  
  cliqueNb=length(cliqueList) ; p=ncol(Sigma)
  # code sans simulation de données
  Corr <- cov2cor(Sigma); sigma <- sqrt(diag(Sigma))
  coef <- matrix(0, p, cliqueNb); lambda <- rep(0, cliqueNb)
  sapply(1:cliqueNb, function(c){
    pca <- eigen(cov2cor(Corr[cliqueList[[c]], cliqueList[[c]]]))
    coef[, c] <<- rep(0, p); 
    coef[cliqueList[[c]], c] <<- pca$vectors[, 1]
    lambda[c] <<- pca$values[1]
  }) 
  
  # Recontructing Sigma
  CorrFull <- rbind(cbind(Corr, Corr%*%coef), cbind(t(coef)%*%Corr, cst*t(coef)%*%Corr%*%coef))
 # isSymmetric(CorrFull); eigen(CorrFull)$values; plot(CorrFull[1:p, 1:p], Corr); abline(0, 1)
  sigmaFull <- c(sigma, rep(1, cliqueNb)) 
  SigmaFull <- diag(sigmaFull) %*% CorrFull %*% diag(sigmaFull)
 # isSymmetric(SigmaFull); eigen(SigmaFull)$values; plot(SigmaFull[1:p, 1:p], Sigma); abline(0, 1)
  
  # Initialising Omega
  OmegaFull <- solve(SigmaFull)
  #isSymmetric(OmegaFull); eigen(OmegaFull)$values
  # sapply(1:cliqueNb, function(c){
  #    OmegaFull[(1:p)[-cliqueList[[c]]], (p+c)] <<- OmegaFull[(p+c), (1:p)[-cliqueList[[c]]]] <<- 0
  # })
  # isSymmetric(OmegaFull); eigen(OmegaFull)$values
  coefDiag <- c(rep(1, p), 1/sqrt(diag(OmegaFull)[p+(1:cliqueNb)]))
  OmegaFull <- diag(coefDiag) %*% OmegaFull %*% diag(coefDiag)
  
  # ajout Rmomal 
##################  
  # X=rmvnorm(n,rep(0,nrow(covX)), covX)
  # #  n<-nrow(X)
  # nbCliques <- length(cliquelist)
  # if(pca){
  #   pr <- lapply(1:nbCliques, function(k)
  #     prcomp(t(X[, cliquelist[[k]]]), center = TRUE, scale = TRUE))
  #   mat <- lapply(1:nbCliques, function(k)
  #     scale(pr[[k]]$rotation[, 1], center = TRUE, scale = TRUE))
  #   newX <- scale(X, center = TRUE, scale = TRUE)
  #   newX <- scale(cbind(newX, do.call(cbind, mat)), center = TRUE, scale = TRUE)
  # }else{
  #   
  #   # mat <- lapply(1:nbCliques, function(k)
  #   #   rowMeans(X[, cliquelist[[k]]]) )
  #   mat <- lapply(1:nbCliques, function(k){
  #     apply(X[, cliquelist[[k]]], 1, function(x) median(x))
  #   }
  #   )
  #   
  #   mat=do.call(cbind, mat)
  #   newX <- scale(cbind(X, mat), center = TRUE, scale = TRUE) # pourquoi scale
  # }
  # 
  # 
  # # pour chacune des cliques identifiées dans X, on résume les variables liées par une acp.
  # # en récupérant le vecteur principal de chaque acp, qui vont jouer le rôle de nouvelles variables
  # # qui étaient manquantes et corrélées à toutes les variables de la clique.
  # mu <- matrix(0, ncol(newX), 1)
  # newX <- newX + mvrnorm(n, mu, 0.01 * diag(ncol(newX)))
  # Sigma0 <- 1 / n * t(newX) %*% (newX)
  # K0 <- solve(Sigma0)
  ############
  return(structure(list(
    Sigma0 = SigmaFull,
    K0 = OmegaFull,
    cliquelist = cliqueList
  )))
}
#Find the clique created by the missing actor
findCliques<-function(covX,k){
  # -----------------------------------------------------------------------------------------------------------------------------
  # FUNCTION
  # Find k groups of nodes which are 'clique likes"
  # INPUT
  #   X      : data matrix (nxp)
  #   k      : number of groups
  #            is determined with hierarchical clustering.
  # OUTPUT
  #   cliquelist    : a list of groups of indices
  # -----------------------------------------------------------------------------------------------------------------------------
  
  Y = 1- abs(cov2cor(covX))
  res.kmeans<-kmeans(Y,k)
  # Compute the mean distance within each cluster (clique should have the smallest mean distances)
  scores<-res.kmeans$withinss / ((res.kmeans$size)*(res.kmeans$size-1)/2)
  # Renumbering of the cluster according the score (cluster one will have the smallest score)
  reordered.cluster<-order(scores)[res.kmeans$cluster]
  return(cliquelist = split(1:length(reordered.cluster),reordered.cluster))
}
init.mclust<-function(S,nb.missing=1, n.noise=50,plot=TRUE, title="",trueClique=NULL){
  #browser()
  #X <-data.frame(t(t(eigen(S)$vectors[,1:2])*sqrt(eigen(S)$values[1:2])))
  Scomp=prcomp(S,scale. = TRUE)
  data=data.frame(Scomp$rotation[,1:2]%*%diag(Scomp$sdev[1:2]))
  
  datapolar=cart2pol(x=data[,1],y=data[,2])[,1:2] # recup r et theta
  
  
  datapolar_half=datapolar %>% mutate(theta2=ifelse(theta>pi,theta-pi,theta)) %>%
    dplyr::select(r,theta2)
  
  # noise in polar coords
  r <- sqrt(runif(n.noise))
  theta2 <- runif(n.noise, 0, pi)
  
  datapolarall=rbind(datapolar_half,cbind(r,theta2)) 
  
  newdata=pol2cart(datapolarall$r,datapolarall$theta2)[,1:2]
  #plot(newdata,col=col, xlim=c(-1,1),ylim=c(-1,1), pch=20)
  
  #  noiseInit<-sample(c(T,F), size=nrow(X)+n.noise, replace=T, prob=c(3, 1))
  noiseInit<-c(rep(F,ncol(S)),rep(T,n.noise))
  clust=Mclust(data=newdata,
               initialization = list(noise=noiseInit),
               G=nb.missing)
  groups<-mclust::map(clust$z)
  res<-which(groups==1)[which(groups==1)<=nrow(S)] 
  # if(length(res)==0)browser()
  
  #   browser()
  if(!is.null(trueClique)){
    #False positives rate
    p=ncol(S)
    N=setdiff(1:p,trueClique)
    FP=sum(res%in%N)/length(N)
    # False negatives rate
    FN=sum(setdiff(1:p,res)%in%trueClique)/length(trueClique)
    title=paste0(title,":"," FN=", round(FN,2),",FP=",round(FP,2))
  }
  if(plot){   
    g= ggplot(data,aes(X1,X2, label=1:p,color=(1:p)%in%res))+geom_point(size=0.1)+
      theme_light()+geom_text()+labs(x="eig vect 1",y="eig vect 2", title=title)+
      guides(color=FALSE)+scale_color_brewer(palette="Dark2")+
      geom_hline(yintercept=0, color="gray50")+geom_vline(xintercept=0, color="gray50")
    print(g)
  }
  return(list(init=res, FPN=c(FN,FP)))
}

#################
# For OPTIM

#=====
# VE step
VE<-function(MO,SO,SH,omega,W,Wg,MH,Pg,maxIter,minIter,eps, 
             alpha,form,beta.min=1e-10, plot=FALSE,verbatim=FALSE){
  t1=Sys.time()
  M=cbind(MO,MH) ;  S<-cbind(SO, matrix(rep(SH,n),n,1))
  LB=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p)$J
  cat(paste0("\n entree VE: ",LB,"\n"))
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO); H=(p+1):ncol(omega);r=length(H); iter=0 ; 
  omegaH=omega[H,H]
  cat(paste0("\n omegaH = ",omegaH,"\n"))
  diffKL=-1 ; diff=c(1000); Wdiff=1; diffW=c(0.1); diffMH=1; diffM=c(0.1)
  
  phi=CorOmegaMatrix(omega)
  SH <- 1/omegaH
  S<-cbind(SO, matrix(rep(SH,n),n,1))# all SHi have same solution, depending only on Mstep
  
  LB=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p)$J
  cat(paste0("\n S: ",LB,"\n"))
  
  KL<-numeric(maxIter)
  
  #  while( iter < maxIter ){#((diffW > eps) && (iter < maxIter)) || (iter < minIter)
  iter=iter+1
  #   cat(iter)
  #-- Probabilities estimates
  # Pg = EdgeProba(Wg)
  
  # if(sum(is.nan(Pg))!=0){
  #   browser()
  #   cat( ": adjust ", summary(c(logWtree)), "\n")
  #   logWtree=log(adjustW(exp(logWtree)))
  #   cat("to: ",summary(c(logWtree)), "\n")
  #   Pg = EdgeProba(exp(logWtree))
  # } 
  #-- Updates
  
  #---- Wg
  Wg.new= computeWg(phi,omega,W,MH,MO,alpha)  
  diag(Wg.new)=1
  Wg.new[which(Wg.new< beta.min)] = beta.min
  Pg = EdgeProba(Wg.new)
  
  
  Wdiff=max(abs(F_Sym2Vec(Wg.new)-F_Sym2Vec(Wg)))
  if(is.nan(Wdiff)) browser()
  
  #---- MH 
  MH.new<- (-MO) %*% (Pg[O,H] * omega[O,H]) #omegaH=1
  diffMH<-max(abs(MH-MH.new))
  M=cbind(MO,MH.new)
  
  #---- Wg
  Wg.new= computeWg(phi,omega,W,MH,MO,alpha )  
  diag(Wg.new)=1
  Wg.new[which(Wg.new< beta.min)] = beta.min
  Pg = EdgeProba(Wg.new)
  
  
  Wdiff=max(abs(F_Sym2Vec(Wg.new)-F_Sym2Vec(Wg)))
  if(is.nan(Wdiff)) browser()
  
  KL[iter]<-argminKL(F_Sym2Vec(log(Wg.new)), Pg, M,S,omega,phi,W,p)
  if(iter>1) diffKL =(KL[iter] - KL[iter-1])
#  browser()
  MH=MH.new
  M=cbind(MO,MH)
  Wg=Wg.new
  #-- end
  if(iter>1){
    diff = c(diff,diffKL)
    diffW=c(diffW,Wdiff)
    diffM=c(diffM,diffMH)
  } 
  # }
  KL=KL[1:iter]
  t2=Sys.time(); time=t2-t1
  if(verbatim) cat(paste0("\nVE step converged in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal W difference: ",round(diffW[iter],5),
                          "\nFinal KL difference: ",round(diff[iter],4)))
  if(plot){
    g=data.frame(Diff.W=diffW,  diff.KL=diff,diff.MH=diffM, PartofKL=KL) %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+ facet_wrap(~key, scales="free")+ labs(x="iter",y="") +
      mytheme.dark
    print(g)
  }
  
  res=list(Gprobs=Pg,Gweights=Wg,Hmeans=M,Hvar=S, KL=KL[iter], diff=diffM, diffW=diffW )
  return(res)
}

argminKL <- function(gamma, Pg, M,S,omega,phi,W,p ){
  r=ncol(omega)-p ; O = 1:p ; H=(p+1):(p+r)
  omegaH=omega[H,H]
  EhZoZo = t(M[,O])%*%M[,O]+ diag(colSums(S[,O]))
  
  KL <- 0.5 * sum(diag(omegaH * (t(M[, H]) %*% M[, H] + sum(S[, H])))) +
    sum(diag(Pg[O, H] * omega[O, H] %*% (t(M[, H]) %*% M[, O]))) +
    0.5 * sum(diag((Pg[O, O] * omega[O, O]) %*% EhZoZo)) -
    n * 0.5 * sum(log(diag(omega))) -
    n * 0.5 * 2 * sum(F_Sym2Vec(Pg) * log(F_Sym2Vec(phi))) +
    2 * sum(F_Sym2Vec(Pg) * (gamma - log(F_Sym2Vec(W)))) -
    log(SumTree(F_Vec2Sym(exp(gamma)))) +
    log(SumTree(W)) - 
    0.5 * sum(log(S[, H]))#-lambda*(sum(exp(gamma)) - 0.5)
  
  return(KL)
}
computeWg<-function(phi,omega,W,MH,MO, alpha ){

  p=ncol(MO) ; r=ncol(MH) ; O = 1:p ; H = (p+1):(p+r)
  logWg<-matrix(0,(p+r),(p+r)) ; Wg<-matrix(0,(p+r),(p+r))
  psi=phi^(n*alpha*0.5)
  
  
  logWg[O,O]<-log(W[O,O])+log(psi[O,O])-0.5*alpha*omega[O,O]*(t(MO)%*%MO)
  logWg[O,H]<-log(W[O,H])+log(psi[O,H])-alpha*omega[O,H]*(t(MH)%*%MO)
  
  
  logWg[H,O]<-t(logWg[O,H])
  diag(logWg) = 0
  # shrinking and centering
  gamma=logWg[O,H]
  gamma[which(gamma<(-30))]=-30
  gamma=gamma-mean(gamma)
  gamma[which(gamma>20)]=20
  Wg[O,H]=exp(gamma)
  Wg[H,O]<-t(Wg[O,H])
  gamma=F_Sym2Vec(logWg[O,O])
  gamma[which(gamma<(-30))]=-30
  gamma=gamma-mean(gamma)
  gamma[which(gamma>max(log(Wg[O,H])))]=max(log(Wg[O,H]))
  Wg[O,O]=exp(F_Vec2Sym(gamma))
  
  diag(Wg) = 1
  
  return(Wg)
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

#=====
# M step
Mstep<-function(M,S,Pg, omega,W,maxIter, beta.min, trim=TRUE,plot=FALSE,eps, verbatim=FALSE,Wg,p){
  t1=Sys.time()
  diffJ=(1); diff.J=c(1);  diff.W=c(0.05);diffW=1
  maxJ=c()
  n=nrow(S) 
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  iter=0
  
  LB=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p)$J
  cat(paste0("\n entree M: ", LB,"\n"))
  
  omegaDiag<-optimDiag(Pg, omega, SigmaTilde)
  
  omega.new=computeOffDiag(omegaDiag,SigmaTilde)
  #  browser()
  LB=LowerBound(Pg = Pg, omega=omega.new, M=M, S=S,W=W, Wg=Wg,p)$J
  cat(paste0("\n omega : ", LB,"\n")) 
  
  
  phi=CorOmegaMatrix(omega.new) # de omega ou omega.new ?
  
  while((diffW>eps && iter <= maxIter) ){
    iter=iter+1
    if(verbatim)  cat(paste0("\nIter n°:",iter))
    # browser()
    Mei=Meila(W)  # c'est Meila qui justifie de faire un while
    W.new= Pg/(Mei) 
    
    diag(W.new)=1
    W.new[which(W.new< beta.min)] = beta.min
    diffW=max(abs(F_Sym2Vec(W.new)-F_Sym2Vec(W)))
    
    maxJ[iter]<-argmaxJ(log(W.new),Pg,omega.new,SigmaTilde,phi,n)
    if(is.nan(maxJ[iter])) browser()
    if(iter>1){
      diffJ = (maxJ[iter] - maxJ[iter-1])
      diff.J = c(diff.J,diffJ)
      diff.W = c(diff.W,diffW)
    } 
    W=W.new
  }
  t2=Sys.time(); time=t2-t1
  if(verbatim) cat(paste0("\nM step converged in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal maxJ difference: ",round(diff.J[iter],4)))
  
  if(plot){
    g=data.frame(Diff.W=diff.W,Diff.J=diff.J, Jbound=maxJ) %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+ mytheme.dark+
      facet_wrap(~key, scales="free")+labs(x="iter",y="") 
    print(g)
  }
  omega=omega.new
  res=list(W=W, omega=omega, diff=diff, diffW=diffW, finalJ=maxJ[iter])
  
  LB=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p)$J
  cat(paste0("\n W: ",LB,"\n"))
  return(res)
}

argmaxJ<-function(gamma,Pg,omega,sigmaTilde,phi,n, trim=TRUE, verbatim=TRUE){
  
  maxJ <- sum(Pg*F_Vec2Sym(F_Sym2Vec(gamma) + n*0.5*F_Sym2Vec(log(phi)))) + 
    n*0.5*sum(log(diag(omega))) - log(SumTree(exp(gamma))) - 
    0.5*n*sum(diag( Pg * omega %*% sigmaTilde))# - n*0.5*sum(diag(omega*sigmaTilde)) #termes diagonaux
  # + lambda*(sum(W)-0.5)
  
  return(maxJ)
}

optimDiag<-function(Pg,omega,SigmaTilde){
  
  quantity<-Pg*omega*SigmaTilde
  diag(quantity)=0
  vectorSuml<-colSums(quantity)
  omegaDiag = (1-vectorSuml)/diag(SigmaTilde)
   
  
  return(omegaDiag)
}
optimDiag2 <- function(omega_ii,i,  omega, SigmaTilde, Pg) {
  # browser()
  q <- nrow(omega)
  phi = CorOmegaMatrix(omega)
  newphi=(1- phi)/phi
  diag(newphi)=0
  quotient = diag(Pg %*% newphi)[i]
  
  grad=(1/omega_ii)*(quotient+1) - SigmaTilde[i, i]
  return( grad)
}
dichotomie <- function(a, b, F.x, eps, bool=FALSE){
  # if(bool) browser()
  x.min = ifelse(F.x(0) >0,b,1e-4);
  while(F.x(x.min)<0){x.min = x.min -x.min/2}
  
  x.max = a
  #browser()
  while(F.x(x.max)>0){x.max = x.max * 2}
  #  cat(paste0("\n [xmin ; xmax] = [",round(x.min,3),",",round(x.max,3),"]"))
  x = (x.max+x.min)/2
  
  f.min = F.x(x.min)
  f.max = F.x(x.max)
  f = F.x(x)
  iter=0
  while(abs(x.max-x.min) > eps){
    iter=iter+1
    
    # if(iter%%10 == 0) cat(paste0(iter,", "))
    # if(iter > 30) browser()
    if(f < 0) {
      x.max = x
      f.max = f
    }else{
      x.min = x
      f.min = f
    }
    x = (x.max+x.min)/2;
    f = F.x(x)
  }
  #  cat(paste0("f=",round(f,8),", om_ii=",round(x,3),"\n"))
  return(x)
}

computeOffDiag<-function(omegaDiag,SigmaTilde){
  q=length(omegaDiag)
  omega=matrix(0,q,q)
  sapply(1:(q-1),
         function(j){
           sapply((j+1):q,
                  function(k){
                    omega[k, j] <<- (1 - ( 1+4*SigmaTilde[k,j]^2*omegaDiag[j]*omegaDiag[k])^0.5)/(2*SigmaTilde[k,j])
                    omega[j, k] <<- omega[k, j]
                  }
           )
         }
  )
  diag(omega)=omegaDiag
  return(omega)
}


#=====
#Lower bound

LowerBound<-function(Pg ,omega, M, S, W, Wg,p){
  n=nrow(M)
  q=nrow(omega)
  r=q-p
  O=1:p
  H=(p+1):q
  psi=CorOmegaMatrix(omega)
  #Egh lop (Z |T)
  # browser()
  T1<- 2*sum(F_Sym2Vec(n*0.5* Pg * log (psi ))) + n*0.5* sum(log(diag(omega))) - 0.5* sum( (Pg*omega) * (t(M)%*%M + diag(colSums(S)) )) - q*n* 0.5 * log(2*pi)
  
  # Eg log(p) - Eg log(g)
  
  T2<-2*sum(F_Sym2Vec(Pg) * (log(F_Sym2Vec(W)) - log(F_Sym2Vec(Wg)) )) - log(SumTree(W))+ log(SumTree(Wg))
  
  #Eh log h(ZH)
  T3<- 0.5*sum(log(S[,H])) - n*r*0.5*(1+log(2*pi))
  
  
  J=T1+T2+T3
  if(is.nan(J)) browser()
  return(data.frame(J=J, "Egh(log pZ|T)"=T1, "Eg(log p/g)"=T2, "Eh(log hZH)"=T3))
}



#####################
# For double hidden proba

HiddenEdgeProba <- function(W,r=1, verbatim=FALSE){
  #computes the probability for two links to a hidden covariates to NOT be there
  #Coded for 1 hidden covariate
  it=-1
  Wcum = SumTree(W)
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
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
  P[which(P<1e-10)]=1e-10  
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


########################################
# GENERAL
F_Vec2Sym <- function(A.vec){
  # Makes a symmetric matrix from the vector made of its lower tirangular part
  n = (1+sqrt(1+8*length(A.vec)))/2
  A.mat = matrix(0, n, n)
  A.mat[lower.tri(A.mat)] = A.vec
  A.mat = A.mat + t(A.mat)
  return(A.mat)
}
F_Sym2Vec <- function(A.mat){
  # Makes a vector from the lower triangular par of a symmetric matrix
  return(A.mat[lower.tri(A.mat)])
}
EMtree_corZ<-function(CovZ,n,  maxIter=30, cond.tol=1e-10, verbatim=TRUE, plot=FALSE){
  CorZ=cov2cor(CovZ)
  p = ncol(CorZ)
  alpha.psi = Psi_alpha(CorZ, n, cond.tol=cond.tol)
  psi = alpha.psi$psi
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)
  FitEM = FitBetaStatic(beta.init=beta.unif, psi=psi, maxIter = maxIter,
                        verbatim=verbatim, plot=plot)
  return(FitEM)
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

generator_PLN<-function(Sigma,covariates=NULL, n=50){
  # ajout d'une constante, pa rapport à EMtree::generator_PLN
  p<-ncol(Sigma)
  if(!is.null(covariates)){
    n<-nrow(covariates)
    c<-ncol(covariates)
    string<-paste0("~", paste(colnames(covariates), collapse=" + "))
    formula<-as.formula(string)
    m<- model.matrix(formula,covariates)[,-1]
    mc<-ncol(m)
    beta<-matrix(runif(p*mc),mc,p)
    prod=m %*% beta+2
  }else{
    prod=2
  }
  Z<- rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
  Y = matrix(rpois(n*p, exp(Z+prod )), n, p)
  return(Y)
}
generator_param<-function(G,signed=FALSE,v=0){
  lambda = 1;  p=ncol(G);  sumlignes=rowSums(matrix(G,p,p));  D=diag(sumlignes+v)
  if(signed){
    Gsign = F_Vec2Sym(F_Sym2Vec(G * matrix(2*rbinom(p^2, 1, .3)-1, p, p)))
    omega = lambda*D - Gsign
    while(min(eigen(omega)$values) < 1e-10 & lambda<1e3){
      lambda = 1.1*lambda
      omega = lambda*D - Gsign
    }
  }else{
    omega = lambda*D + G
    while (min(eigen(omega)$values) < 1e-10){
      lambda = 1.1*lambda
      omega =lambda*D + G
    }
  }
  sigma = cov2cor(solve(omega))
  sim=list(sigma=sigma,omega=omega,cste=lambda)
  return(sim)
}
#rewrite data from scratch
library(huge)
library(Matrix)
data_from_scratch<-function(type, p=20,n=50, r=5, covariates=NULL,
                            prob=log(p)/p,dens=log(p)/p,
                            signed=FALSE,v=0,draw=FALSE){
  # make graph
  graph<- generator_graph(graph=type,p=p,prob=prob,dens=dens,r=r)
  param<-generator_param(G=as.matrix(graph),signed=signed,v=v)
  data<-generator_PLN(param$sigma,covariates,n)
  if(draw){
    g=as_tbl_graph(as.matrix(graph)) %>%
      ggraph(layout="nicely")+
      geom_edge_link()+
      geom_node_point(size=2, color="blue")
    print(g)
  }
  
  return(list(data=data,omega= param$omega))
}

################
# diagnostic functions
auc<-function(pred,label){ #require(ROCR)
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  ROC_auc <- round(performance(prediction,"auc")@y.values[[1]],digits=2)
  return(ROC_auc)
}
accppvtpr<-function(probs,omega,h, seuil=0.5){
  Acc=round(mean(1*(probs>seuil)==omega),2) #(TP+TN)/(TP+FP+TN+FN)
  AccH=round(mean(1*(probs[h,]>seuil)==omega[h,]),2)
  AccO=round(mean(1*(probs[-h,-h]>seuil)==omega[-h,-h]),2)
  PPV=round(sum((omega!=0)*(probs>seuil))/(sum((omega!=0)*(probs>seuil))+ sum((omega==0)*(probs>seuil))),2)#TP/(TP+FP)
  PPVH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/(sum((omega[h,]!=0)*(probs[h,]>seuil))+ sum((omega[h,]==0)*(probs[h,]>seuil))),2)
  PPVO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))+  sum((omega[-h,-h]==0)*(probs[-h,-h]>seuil))),2)
  
  TPR=round(sum((omega!=0)*(probs>seuil))/sum(omega!=0), 2)
  TPRH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/sum(omega[h,]!=0), 2)
  TPRO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/sum(omega[-h,-h]!=0), 2)
  
  return(c(Acc,AccH,AccO,PPV,PPVH,PPVO,TPR,TPRH,TPRO))
}
ggimage<-function(data){
  melted_data <- melt(data)
  ggplot(melted_data, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +guides(fill=FALSE)+ theme(plot.title = element_text(size=10, hjust=0.5))
}
plotVEM<-function(probs,omega,r,seuil){
  # plots heatmaps for the chosen threshold and print verdicts as title
  q=ncol(omega)
  h=(q-r):q
  performance=accppvtpr(probs,omega,h,seuil)
  Acc=performance[1] ;AccH=performance[2] ;AccO=performance[3] 
  PPV=performance[4] ;PPVH=performance[5] ; PPVO=performance[6]
  TPR=performance[7] ;TPRH=performance[8] ;TPRO=performance[9] 
  p1<-ggimage(probs)+labs(title=paste0("G hat"))
  p2<-ggimage(omega)+labs(title="True G")
  auc<-round(auc(pred = probs, label = omega),3)
  grid.arrange(p1,p2,ncol=2, top=paste0("Tpr=",TPR," (TprO=",TPRO," , TprH=",TPRH,
                                        ")\n Ppv=",PPV," (PpvO=",PPVO," , PpvH=",PPVH,")",
                                        "\n AUC=",auc))
  
}

courbes_seuil<-function(probs,omega,h,seq_seuil){
  
  tmp=sapply(seq_seuil,function(x)  accppvtpr(seuil=x,probs=probs,omega=omega,h=h))
  res=data.frame(cbind(t(tmp),seq_seuil))
  colnames(res)=c("Acc","AccH","AccO","PPV","PPVH","PPVO","TPR","TPRH","TPRO","seuil")
  return(as_tibble(res))
  
}
plotVerdict<-function(values,colonne){
  # plots verdict curves along threshold from VEM result
  colonne<-enquo(colonne)
  values %>%as_tibble() %>% gather(key,value, - !!colonne)%>%
    mutate(status=ifelse(substr(key,4,4)=="","all",substr(key,4,4))) %>% 
    mutate(key=substr(key,1,3)) %>% spread(key,value)%>% 
    mutate(PPV=ifelse(PPV<0,NA,PPV)) %>% 
    ggplot(aes(TPR,PPV,color=status))+
    geom_rect(aes(xmin=0.5, xmax=1, ymin=0.5, ymax=1), fill="gray90", color="gray90",alpha=0.1)+
    geom_point()+  geom_line()+
    facet_wrap(~status)+mytheme.dark+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


VEMtree<-function(counts,MO,SO,sigma_obs,ome_init,W_init,Wg_init, verbatim=TRUE,maxIter=20, plot=TRUE, eps=1e-2, alpha){
  # MH = matrix(100,n,r);
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO); H=(p+1):ncol(omega);r=length(H);
  pr=prcomp(t(counts),scale. = FALSE)
  MH = matrix(pr$rotation[,1]*pr$sdev[1],nrow=n,ncol=r)
  omega=omega_init;  W=W_init;  Wg=Wg_init
  iter=0 ; lowbound=list()
  KL=c() ; J=c();diffW=c();diffOm=c()
  t1=Sys.time()
  Pg=matrix(0.5, ncol(W),ncol(W))
  SH=matrix(1,nrow=n,ncol=r)
  while(iter < maxIter){ #(diffW[iter] > eps && diffOm[iter] > eps && iter < maxIter) || iter<1
    
    iter=iter+1 
    cat(paste0("\n Iter n°", iter))
    #VE
    resVe<-VE(MO,SO,SH, sigma_obs,omega,W,Wg,MH=MH,Pg=Pg,maxIter=1,minIter=1,eps=1e-3, plot=FALSE, 
              form="theory",alpha=alpha, verbatim=FALSE)
    KL[iter]=resVe$KL
    M=resVe$Hmeans ; MH=matrix(M[,H],n,r)
    S=resVe$Hvar ; SH=matrix(S[,H],n,r)
    Pg=resVe$Gprobs
    Wg=resVe$Gweights
    #M
    resM<-Mstep(M,S,Pg, omega,W,maxIter=5, beta.min=1e-6, eps=1e-3 ,plot=FALSE, verbatim=FALSE,
                Wg=Wg, p=p)
    W.new=resM$W
    diffW[iter]=abs(max(W.new-W))
    W=W.new
    omega.new=resM$omega
    diffOm[iter]=abs(max(omega.new-omega))
    omega=omega.new
    J[iter]=resM$finalJ
    
    
    lowbound[[iter]] = LowerBound(Pg ,omega, M, S, W, Wg,p)
  }
  lowbound=do.call(rbind,lowbound)
  features<-data.frame(Jbound=J, KL=KL, diffW=diffW, diffOm=diffOm)
  t2=Sys.time()
  t2=Sys.time(); time=t2-t1
  if(verbatim) cat(paste0("\nVEMtree ran in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal Jbound difference: ",round(J[iter]-J[iter-1],5)))
  if(plot){
    g1<-features %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line() + facet_wrap(~key, scales="free")+
      labs(x="",y="", title="Stoping criteria")+ mytheme.dark
    
    
    g2<- lowbound %>% rowid_to_column() %>% gather(key,value,-rowid) %>% 
      ggplot(aes(rowid,value, color=key))+geom_point()+geom_line()+
      facet_wrap(~key, scales="free")+
      labs(x="iteration",y="", title="Lower bound and components")+mytheme
    
    grid.arrange(g1, g2, ncol=1)
  }
  
  
  return(list(M=M,S=S,Pg=Pg,Wg=Wg,W=W,omega=omega, lowbound=lowbound, features=features))
}


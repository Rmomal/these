### fonction krishner
source("/home/momal/Git/these/pack1/R/fonctions.R")
source("/home/momal/Git/these/pack1/R/psi.R")


############
#FONCTIONS
laplacien<-function(matrix){
  lapla<--(matrix)
  diag(lapla)<-rowSums(matrix)
  return(lapla)
}

#laplacien(matrix(1,nrow=5,ncol=5))

# /!\ se rappeler qu'on enleve le coin haut gauche
inv_lap<-function(matrix){
  L<-laplacien(matrix)[-1,-1]
  if (determinant(L,logarithm=FALSE )$modulus[1]>2e-16){
    inv<-solve(L,tol = 1e-17)
    return(inv)
  }else{
      message("not inversible")
    }

}
# beta<-Beta
#Kirshner retourne la matrice des probas
Kirshner<-function(beta){
  #browser()
  if(length(inv_lap(beta))==(ncol(beta)-1)^2){
      Q<-inv_lap(beta)
  colQ<-matrix(diag(Q),nrow=length(diag(Q)),ncol=length(diag(Q)),byrow=TRUE)
  rowQ<-matrix(diag(Q),nrow=length(diag(Q)),ncol=length(diag(Q)),byrow=FALSE)
  matrixM<-(colQ+rowQ-2*Q)
  matrixM<-rbind(diag(Q),matrixM)
  matrixM<-cbind(c(0,diag(Q)),matrixM)
  matrixK<-beta*matrixM
  # matrix<-beta*Meila(beta)
  return(list(matrixK,matrixM))
  }else{
    message("Kirshner not run")
}

}

Kirshner(matrix(-1,nrow=10,ncol=10))

MTT<-function(matrix){
  logdet<-determinant( laplacien(matrix)[-1,-1],logarithm=TRUE)$modulus
  return(logdet)
}

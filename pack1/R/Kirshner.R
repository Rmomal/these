### fonction krishner
source("/home/momal/Documents/codes/fonctions.R")
source("/home/momal/Documents/codes/psi.R")

# simu beta
n<-ncol(data)
Beta<-simMatrix(n)#beta de taille n*n


laplacien<-function(matrix){
  lapla<--(matrix)
  diag(lapla)<-rowSums(matrix)
return(lapla)
}


# /!\ se rappeler qu'on enleve le coin haut gauche
inv_lap<-function(matrix){
  inv<-(laplacien(matrix)[-1,-1])^(-1)
  return(inv)
}

#Kirshner retourne la matrice des probas
Kirshner<-function(beta){
  #browser()
  Q<-inv_lap(beta)
  colQ<-matrix(diag(Q),nrow=length(diag(Q)),ncol=length(diag(Q)),byrow=TRUE)
  rowQ<-matrix(diag(Q),nrow=length(diag(Q)),ncol=length(diag(Q)),byrow=FALSE)
  matrix<-beta[-1,-1]*(colQ+rowQ-2*Q)
  matrix<-rbind(beta[-1,1]*rowQ[,1],matrix)
  matrix<-cbind(c(0,beta[1,-1]*colQ[1,]),matrix)
  # matrix<-beta*Meila(beta)
return(matrix)
}
#retourne la matrice M
Meila<-function(beta){
  #browser
  Q<-inv_lap(beta)
  colQ<-matrix(diag(Q),nrow=length(diag(Q)),ncol=length(diag(Q)),byrow=TRUE)
  rowQ<-matrix(diag(Q),nrow=length(diag(Q)),ncol=length(diag(Q)),byrow=FALSE)
  matrix<-(colQ+rowQ-2*Q)
  matrix<-rbind(rowQ[,1],matrix)
  matrix<-cbind(c(0,colQ[1,]),matrix)
  return(matrix)
}
psi<-calcul_psi(as.matrix(data))
Beta<-simMatrix(11)*exp(psi)
prob_cond<-Kirshner(Beta)

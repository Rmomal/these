source("/home/momal/Git/these/pack1/R/fonctions.R")
source("/home/momal/Git/these/pack1/R/psi.R")
source("/home/momal/Git/these/pack1/R/Kirshner.R")
#Beta<-simMatrix(n)#beta matrice n*n symétrique de probas
library(readxl)

##
#fonctions de traitement
trtmt1<-function(matrix, tronc){
  index<- which(abs(matrix)>tronc,arr.ind=TRUE)
  if(length(index!=0)){
    matrix[index]<-sign(matrix[index])*tronc
    matrix<-(matrix-mean(matrix)+2)/sd(matrix)
    matrix<-matrix-min(matrix)
  }
  return(matrix)
}

trtmt2<-function(matrix, tronc){
  index<- which(abs(matrix)>tronc,arr.ind=TRUE)
  if(length(index!=0)){
  matrix[index]<-sign(matrix[index])*tronc
  matrix<-matrix-min(matrix)+1
  matrix<-log(matrix)
  }
  return(matrix)
}

#####
## initialisation
Y<-as.matrix(read_excel("~/Documents/codes/Data/Data Files/1. cd3cd28.xls"))[1:50,]
n<-ncol(Y)
#Beta<-matrix(1,nrow=n,ncol=n)
rho<-matrix(nrow=n,ncol=n)
for(i in 1:n){
  rho[i,]<-unlist(apply(Y,2,function(x) cor(x,Y[,i])))
}
Beta<-rho-min(rho)

logpsi<-calcul_psi(Y)
#logpsi<-logpsi-max(logpsi)
beta_psi<-(exp(logpsi)^(1/100))*Beta
#beta_psi[ which(beta_psi>(1e+15),arr.ind=TRUE)]<-(1e+15)
heatmap(Kirshner(beta_psi)[[1]])
criterion<-c()
likelihood<-c()
tronc<-1e10/10^(ceiling(log10(max(exp(logpsi)))))
crit<-FALSE

tr<-trtmt2 #choix du traitement

#Beta<-tr(Beta,tronc)

#####
# Algorithme EM
while(!crit){
  beta_psi<-exp(logpsi)*Beta
  
  likelihood<-c(likelihood,MTT(beta_psi)/MTT(Beta))
  #E
  proba_cond<-Kirshner(beta_psi)[[1]]
  length(proba_cond[proba_cond<0])==0
  #M
  Beta_it<-proba_cond/Kirshner(Beta)[[2]]
  diag(Beta_it)<-0
  
  #Beta_it<-tr(Beta_it,tronc)
  criterion<-c(criterion,mean((Beta-Beta_it)^2))
  crit<-mean((Beta-Beta_it)^2)<1e-3|length(criterion)>600
  
  Beta<-Beta_it
}

#####
## Diagnostique
  hist(Beta)
  hist((beta_psi))
  plot(log(criterion))
  plot(likelihood)




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

trtmt3<-function(matrix, tronc){
 max<-max(log(matrix))-100/ncol(matrix)
 matrix<-log(matrix)-max
  return(list(matrix,-max))
}
trtmt3(beta_psi)
#####
## initialisation
Y<-as.matrix(read_excel("~/Documents/codes/Data/Data Files/1. cd3cd28.xls"))[1:50,]

#Y<-X
n<-ncol(Y)
#Beta<-matrix(1,nrow=n,ncol=n)
rho<-matrix(nrow=n,ncol=n)
for(i in 1:n){
  rho[i,]<-unlist(apply(Y,2,function(x) cor(x,Y[,i])))
}
Beta<-rho-min(rho)

logpsi<-calcul_psi(Y)
#logpsi<-logpsi-max(logpsi)
alpha<-1/100
beta_psi<-(exp(logpsi)^(alpha))*Beta
#beta_psi[ which(beta_psi>(1e+15),arr.ind=TRUE)]<-(1e+15)
heatmap(Kirshner(beta_psi)[[1]])
criterion<-c()
likelihood<-c()
tronc<-1e10/10^(ceiling(log10(max(exp(logpsi)))))
crit<-FALSE

tr<-trtmt2 #choix du traitement

Beta<-tr(Beta,tronc)

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
  crit<-mean((Beta-Beta_it)^2)<1e-5|length(criterion)>600

  Beta<-Beta_it
}

#####
## Diagnostique
  hist(Beta)
  hist((beta_psi))
  plot(log(criterion))
  plot(likelihood)





# Simu erdos
library(LITree)

erdos.graph <- graphModel$new(type = "erdos",size=5, p.or.m = 0.5)
plot(erdos.graph,edge.arrow.size=0, vertex.color="gold",
     vertex.size=15, vertex.frame.color="gray",
     vertex.label.color="black",
     vertex.label.cex=0.8,
     vertex.label.dist=0)
model <- GGMmodel$new(graph=erdos.graph)
model$randomSample(n=20)
X=model$getX()
K=model$K
Sigma=model$Sigma

save(X,K,Sigma,file="Erdos20ind5var.Rdata")

load("Erdos20ind5var.Rdata")

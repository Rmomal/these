### mise à jour beta

source("/home/momal/Documents/codes/fonctions.R")
source("/home/momal/Documents/codes/psi.R")
source("/home/momal/Documents/codes/Kirshner.R")


#on normalise en ligne
standard<-function(data){
  res<-t(apply(data,1,function(x) (x-mean(x))/sd(x)))
  return(res)
}

library(readxl)
data <- read_excel("~/Documents/codes/Data/Data Files/1. cd3cd28.xls")
data<-standard(data)
d<-ncol(data)
psi<-calcul_psi(as.matrix(data))
Beta<-simMatrix(d)#beta de taille d*d
Beta_psi<-Beta*exp(psi)

prob_cond<-Kirshner(Beta_psi)
# problème dans les résultats, hyp : mauvais calcul du psi
constante<-determinant.matrix(laplacien(Beta),logarithm = FALSE)$modulus

# while(convergence = FALSE){
#  Beta<-prob_cond/Meila(Beta)
# }


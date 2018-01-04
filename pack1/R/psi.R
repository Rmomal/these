# calcul des psi_kl
library(readxl)
data <- read_excel("~/Documents/codes/Data/Data Files/1. cd3cd28.xls")


###########
# @renvoi une matrice de terme général la somme sur les échantillons du log des psi_kl

#' Title
#'
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
calcul_psi<-function(Y){

  n<-nrow(Y)
  d<-ncol(Y)

  rho<-matrix(nrow=d,ncol=d)
  for(i in 1:d){
    rho[i,]<-unlist(apply(Y,2,function(x) cor(x,Y[,i])))
  }
  prody<-array(0,dim=c(d,d,n))
  for(i in 1:d){
    for(j in 1:d){
      prody[i,j,]<- Y[,j]*Y[,i]
    }
  }
  sumy<-array(0,dim=c(d,d,n))
  for(i in 1:d){
    for(j in 1:d){
      sumy[i,j,]<-(Y[,j]^2+Y[,i]^2)
    }
  }
  tmp<-array(0,dim=c(d,d,n))
  for(i in 1:n){
    tmp[,,i]<-log((matrix(1,nrow=d,ncol=d)-rho^2)^(-1/2))+(rho/(matrix(1,nrow=d,ncol=d)-rho^2))*(prody[,,i]-(rho/2)*sumy[,,i])
  }
  #browser()
  psi<-apply(tmp, MARGIN=c(1, 2), function(x) sum(x,na.rm=TRUE))
  return(psi)
}

psi<-calcul_psi(as.matrix(data))


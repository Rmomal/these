# calcul des psi_kl
###########
# OBSOLETE
###########
# calcul_psi<-function(Y){
#
#   n<-nrow(Y)
#   d<-ncol(Y)
#
#   #calcul de la matrice des rho
#   rho<-matrix(nrow=d,ncol=d)
#   for(i in 1:d){
#     rho[i,]<-unlist(apply(Y,2,function(x) cor(x,Y[,i])))
#   }
#   #matrice 3D des produits de colonnes deux à deux de Y
#   prody<-array(0,dim=c(d,d,n))
#   for(i in 1:d){
#     for(j in 1:d){
#       prody[i,j,]<- Y[,j]*Y[,i]
#     }
#   }
#   #matrice 3D des sommes de carré de colonnes 2 à 2 de Y
#   sumy<-array(0,dim=c(d,d,n))
#   for(i in 1:d){
#     for(j in 1:d){
#       sumy[i,j,]<-(Y[,j]^2+Y[,i]^2)
#     }
#   }
#   #matrice 3D application formule log(psi_kl)
#   tmp<-array(0,dim=c(d,d,n))
#   for(i in 1:n){
#     tmp[,,i]<-log((matrix(1,nrow=d,ncol=d)-rho^2)^(-1/2))+(rho/(matrix(1,nrow=d,ncol=d)-rho^2))*(prody[,,i]-(rho/2)*sumy[,,i])
#   }
#   #browser()
#
#   return(exp(tmp))
# }
#beta_psi<-structure(apply(psi, , function(x) x*Beta), dim=dim(psi))

###########

#' Title
#'
#' @param Y
#' @return matrice de terme général la somme sur les échantillons des log de psi_kl
#' @export
#' @examples

calcul_psi<-function(Y){
  n<-nrow(Y)
  d<-ncol(Y)
  #calcul de la matrice des rho
  rho<-matrix(nrow=d,ncol=d)
  for(i in 1:d){
    rho[i,]<-unlist(apply(Y,2,function(x) cor(x,Y[,i])))
  }
  #matrice 3D application formule log(psi_kl)
  psi<-(-n/2)*log(matrix(1,nrow=d,ncol=d)-rho^2)
  diag(psi)<-0
  #psi<-psi-max(psi) #astuce numérique
  return(psi)
}


source("/home/momal/Git/these/pack1/R/fonctions.R")
source("/home/momal/Git/these/pack1/R/psi.R")
source("/home/momal/Git/these/pack1/R/Kirshner.R")

library(readxl)
library(gridExtra,grid,ggplot2)
library(gplots)
library(lattice)
library(LITree)
library(tidyr,dplyr)
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

trtmt3<-function(matrix, seuil){
 max<-max(matrix)-seuil/ncol(matrix)
 matrix<-matrix-max
 index<-which(matrix<=0,arr.ind=TRUE)
  return(list(index,-max))
}


#####
## Divers

initi_beta<-function(data){
  d<-ncol(data)
  rho<-matrix(nrow=d,ncol=d)
  for(i in 1:d){
    rho[i,]<-unlist(apply(data,2,function(x) cor(x,data[,i])))
  }
  beta<-rho-min(rho)
  return(beta)
}

pal <- colorRampPalette(c("linen", "violetred1", "red"))(n = 49)

draw<-function(frame,Beta,alpha,seuil,tronc,init,msg){
  var1<-frame[,2]
  var2<-frame[,3]
  p <- ggplot(frame,aes(iterations,log(var1)))+
    geom_point()+
    labs(xlab="log(criterion)",title="Convergence criterion")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  p2 <-  ggplot(frame,aes(iterations,var2))+
    geom_point()+
    labs(xlab="log(likelihood)",title="Likelihood")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  r <-ggplot(gather(data.frame(Beta)), aes(value)) +
    geom_histogram(fill="lightseagreen",bins=length(gather(data.frame(Beta))$value)%/%5)+
    labs(title="Edges weights")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))


  t <- textGrob(paste0(" alpha = ",alpha,",\n threshold = ",seuil,",\n troncature = ",tronc,", \n iterations = ",
                       length(var1),",\n min criterion = ",round(min(var1),digits=6),
                       ",\n Beta initialised : ",init,
                       ",\n message : ",msg),just="center")
  return(grid.arrange(t, p, p2, r, ncol=2))
}

decalage<-function(data){
  data2<-data[is.finite(data)]
  if((max(data2)-7*sd(data2))>min(data2)){
   data[which(data2>(max(data2)-7*sd(data2)))]<-max(data2[which(data2<(max(data2)-7*sd(data2)))])
  }
  return(data)
}

decalage(logbeta_psi)
######################

EM_mixTree<-function(Y,alpha,seuil,tronc,init){

## INIT
  if (init){
     Beta<-initi_beta(Y)
  }else{
    n<-ncol(Y)
    Beta<-matrix(1,nrow=n,ncol=n)
  }
  logpsi<-calcul_psi(Y)
  cst<-max(logpsi)
  logpsi<-logpsi-cst
  beta_psi<-(exp(logpsi)^(alpha))*Beta

  criterion<-c()
  loglikelihood<-c()
  msg<-"no problem"
  crit<-FALSE
  crit2<-FALSE
  crit3<-FALSE
  crit4<-FALSE
  tr<-trtmt3 #choix du traitement

## Algorithme EM
  while(!crit){
        #-- ajustement donnees beta_psi
    #browser()
    beta_psi<-(exp(logpsi))*Beta
    beta_psi<-decalage(beta_psi)

    logbeta_psi<-log(Beta)+logpsi
    max<-max(logbeta_psi)
    logbeta_psi<-logbeta_psi-max
    beta_psi<-exp(logbeta_psi)

    logbeta_psi<-log(beta_psi)
    index<-tr(logbeta_psi,300)[[1]]



    beta_psi[index]<-1e-5
    beta_psi2<-beta_psi*exp(cst)# recorrection pour la vraisemblance
    beta_psi2[index]<-tronc
        #-- contrôle du déterminant (limite fixée à précision machine)
   #browser()
    if(determinant(laplacien(beta_psi)[-1,-1],logarithm=FALSE)$modulus[1]>1e-16){
       proba_cond<-Kirshner(beta_psi/max(beta_psi))[[1]]

    ## E step
       #-- contrôle de la positivité des probas (plus stringent, possible avec un det à 1e-6)
      if(length(which(proba_cond<0))!=0){
        crit2<-TRUE
        message("negative probability")
        msg<-"hard inversion"
      }else{

      ## M step
        #browser()
        Beta_it<-proba_cond/Kirshner(Beta)[[2]]
        diag(Beta_it)<-0
          #-- critere de convergence = MSE
        loglikelihood<-c(loglikelihood,MTT(beta_psi2)-MTT(Beta)) #log=TRUE dans det de MTT
        print(loglikelihood[length(loglikelihood)])
        criterion<-c(criterion,mean((Beta-Beta_it)^2))
        crit3<-mean((Beta-Beta_it)^2)<1e-5
        Beta<-Beta_it*10^(ceiling(log10(min(Beta_it)+10)))

      }
    }else{
         message("beta_psi/max not invertible")
         msg<-"not invertible"
      crit4<-TRUE
    }
    crit<- crit3|crit4|length(criterion)>500|crit2
  }

  if(length(criterion)>0){
      res<-data.frame(iterations=seq(1:length(criterion)),criterion=criterion,
                  loglik=loglikelihood)
      draw(res,Beta,alpha,seuil,tronc,init,msg)
  }
return(Beta)
}


#####
## Diagnostique

# deux jeux de donnees : réel et simulé de LITree::erdos
Y<-as.matrix(read_excel("~/Documents/codes/Data/Data Files/1. cd3cd28.xls"))[1:100,]
load("Erdos20ind5var.Rdata")

A<-EM_mixTree(Y,1/100,300,1e-1,TRUE)
heatmap.2(A, Rowv=NA,Colv=NA, density.info="none", trace="none", dendrogram="none",
          symbreaks=F, scale="none",breaks=50,col=pal)


S<-EM_mixTree(Y,1,100,0.001,TRUE)


B<-EM_mixTree(X,1,300,1e-1,TRUE)
heatmap.2(S, Rowv=NA,Colv=NA, density.info="none", trace="none", dendrogram="none",
          symbreaks=F, scale="none",breaks=50,col=pal)

heatmap(K, Rowv=NA,Colv=NA)
# il ajoute une arête

#####
# faire varier alpha
# faire varier seuil
# faire le graph de beta




######################
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



removeDiag<-function(M){
  diag(M)<-NA
  return(matrix(M[!is.na(M)],nrow(M),ncol(M)-1))
}


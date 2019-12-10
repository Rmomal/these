
buildGraph <- function(graph.type=c("sparse",
                                    "GGMselect",
                                    "doublestar",
                                    "hmm",
                                    "regular",
                                    "erdos",
                                    "starerdos",
                                    "complete",
                                    "tree"),
                       p, eta=0.2, extraeta=eta/5,
                       p.or.m=0.02, erdos.type="gnp"){
  
  level_m=p
  switch(graph.type,
         sparse={
           adjmat <- matrix(0,nrow=p,ncol=p)
           for (i in 1:(p-1)){
             adjmat[i,i+1] <- 0.5
             adjmat[i+1,i] <- 0.5
           }
           print("sparse graph generated")
         },
         dense={
           adjmat       <- matrix(1/2,nrow=p,ncol=p)
           diag(adjmat) <- 0
           print("dense graph generated")
         },
         GGMselect={
           if (mode(p) != "numeric" || length(p) != 1 ||
               p <= 1 || (p%%1) !=   0)
             stop("bad value of p")
           if (mode(eta) != "numeric" || length(eta) != 1 ||
               eta < 0 || eta > 1)
             stop("bad value of eta")
           if (mode(extraeta) != "numeric" || length(extraeta) != 1 ||
               extraeta < 0 || extraeta > 1)
             stop("bad value of extraeta")
           # Constants
           EPS      <-  1/10
           PCORMIN  <-  1/1000
           PCORMAX  <- 5/1000
           N        <- p*(p-1)/2
           
           
           Nab                       <- rbinom(1,N,extraeta)
           temp                      <- rep(0,N)
           nonzeros                  <- sample(1:N,Nab)
           temp[nonzeros]            <- 1
           dimlab                    <- as.character(1:p)
           adjmat                    <- matrix(0,
                                               nrow=p,
                                               ncol=p,
                                               dimnames=list(dimlab, dimlab))
           adjmat[lower.tri(adjmat)] <- temp
           
           d    <- floor(p/3)
           D    <- d*(d-1)/2
           temp <- vector("integer", length=D)
           
           for (i in 1:3) {
             Nab             <- rbinom(1,D,eta)
             temp[]          <- 0
             nonzeros        <- sample(1:D,Nab)
             temp[nonzeros]  <- 1
             G               <- matrix(0,nrow=d,ncol=d)
             G[lower.tri(G)] <- temp
             ind            <- ((1+(i-1)*d):(i*d))
             adjmat[ind,ind] <- G
           }
           print("GGMselect graph generated")
         },
         # simone={
         #   alpha    <- c(eta, eta, eta)
         #   pi       <- matrix(extraeta, 3, 3)
         #   diag(pi) <- eta
         #   net      <- rNetwork(p, pi, alpha, directed=FALSE, signed=TRUE)
         #   adjmat   <- net$A/2
         #   print("simone graph generated")
         # },
         star={
           M              <-  p+1
           adjmat         <- matrix(0,M,M)
           adjmat[1:p,M]  <-  1
           print("star tree built")
         },
         doublestar={
           M                                <-  p+2
           adjmat                           <-  matrix(0,M,M)
           adjmat[1:ceiling(p/2),M-1]       <- 1
           adjmat[(ceiling(p/2)+1):(M-1),M] <- 1
           print("double star tree built")
         },
         hmm ={
           M                            <- 2*p-2;
           adjmat                       <- matrix(0,M,M)
           adjmat[1,p+1]                <- 1
           adjmat[p,M]                  <- 1
           adjmat[2:(p-1),(p+1):M]      <- diag(p-2)
           adjmat[(p+1):(M-1),(p+2):M]  <- diag(p-3)
           print("hmm tree built")
         },
         regular={
           adjmat    <- matrix(0,p,p)
           num_nodes <- p
           
           while(num_nodes>2){
             num_p      <- floor(num_nodes/3)
             new_adjmat <- kronecker(diag(num_p), t(c(1,1,1)))
             res_node   <- num_nodes-3*num_p
             if(res_node==1){
               vec        <- matrix(0,num_p,1)
               vec[num_p] <- 1
               new_adjmat <- cbind(new_adjmat, vec)
             } else if(res_node==2){
               vec            <- matrix(0,num_p,2)
               vec[num_p,1:2] <- 1
               new_adjmat     <- cbind(new_adjmat, vec)
             }
             adjmat    <- rbind(adjmat, cbind(matrix(0,num_p,ncol(adjmat)-ncol(new_adjmat)), new_adjmat))
             vec       <- matrix(0,nrow(adjmat),num_p)
             adjmat    <- cbind(adjmat, vec)
             num_nodes <- num_p
             level_m   <- rbind(level_m,nrow(adjmat))
           }
           print('regular tree built')
         },
         erdos={
           g <- erdos.renyi.game(n=p, p.or.m=p.or.m, type=erdos.type, directed=FALSE)
           adjmat <- 0.5*get.adjacency(g, type="both", sparse=FALSE)
           print('Erdös-Renyi graph built')
         },
         starerdos={
           g <- erdos.renyi.game(n=p, p.or.m=p.or.m, type=erdos.type, directed=FALSE)
           adjmat <- 0.5*get.adjacency(g, type="both", sparse=FALSE)
           adjmat[1,2:p]<-0.5
           adjmat[2:p,1]<-0.5
           print('Starred Erdös-Renyi graph built')
         },
         
         complete={
           adjmat <- matrix(0.5, p, p)
           print('complete graph built')
         },
         tree={
           adjmat <- matrix(1, p, p)
           Sigma  <- sampleGraph(adjmat)$S
           adjmat <- 0.5*ChowLiu2(Sigma=Sigma)
           print('Random Chow-Liu tree built')
         },
         stop("Graph not generated !")
  )
  adjmat <- adjmat+t(adjmat)
  return(adjmat)
} # end of buildGraph


# -------------------------------------------------------------------------------------------------------------------
# FUNCTION
#   Build graph from precision matrix for gaussian models by extracting the support
## INPUT
##  K: precision matrix
## OUTPUT
##  adjmat :  adjacency matrix
# -------------------------------------------------------------------------------------------------------------------
graphFromPCor <- function(K){
  adjmat       <- (K!=0)*1
  diag(adjmat) <- 0
  return(adjmat)
} #


covarianceFromGraph <- function(adjmat, prop.positive.cor=1,
                                missing.var.list=NULL,
                                alpha.hidden=1.25,
                                alpha.observed=1.05,
                                method=NULL){
  # ---------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   Construct SPD covariance matrix whose inverse has specified support
  # INPUT
  #   adjmat : adjacency matrix
  # OUTPUT
  #   K  : precision matrix
  # ---------------------------------------------------------------------------------------------------------
  p      <- nrow(adjmat)
  params <- runif(p*p, min=0.8, max=1) * (2*rbinom(p*p,size=1,prob=1-prop.positive.cor)-1) # Define the proportion of negative versus positive interaction
  
  params <- matrix(params, nrow=p, byrow=TRUE)
  params <- (params+t(params))/2
  
 
  #####  default ######
  if(method=="default"){
    params <- params*adjmat
    if (!is.null(missing.var.list)) {
      #KOO
      params[-missing.var.list,-missing.var.list] <- alpha.observed*params[-missing.var.list,-missing.var.list]
      #KHO
      params[missing.var.list,-missing.var.list] <- alpha.hidden*params[missing.var.list,-missing.var.list]
      #KOH
      params[-missing.var.list,missing.var.list] <- alpha.hidden*params[-missing.var.list,missing.var.list]
    }
    
 #   diag(params) <- 1.1 * rowSums(abs(params))
    mineig <- 0.1+0.8*runif(1,min=0,max=1)
    value=mineig- min(eigen(params, symmetric=TRUE)$values)
    if(value>0){
      K  <- params+(value)*diag(p)
    }else{
      K  <- params+(-value+mineig)*diag(p)
    }
    
  #  K  <- params
    
  }
  if(method=="raph"){
    params <- params*adjmat
    if (!is.null(missing.var.list)) {
      #KHH
      params[missing.var.list,missing.var.list] <- alpha.hidden*params[missing.var.list,missing.var.list]
      #KHO
      params[missing.var.list,-missing.var.list] <- alpha.observed*params[missing.var.list,-missing.var.list]
      #KOH
      params[-missing.var.list,missing.var.list] <- t(  params[missing.var.list,-missing.var.list] )
    }
    
    #   diag(params) <- 1.1 * rowSums(abs(params))
    mineig <- 0.1+0.8*runif(1,min=0,max=1)
    value=mineig- min(eigen(params, symmetric=TRUE)$values)
    if(value>0){
      K  <- params+(value)*diag(p)
    }else{
      K  <- params+(-value+mineig)*diag(p)
    }
    Kh  <- K[missing.var.list,missing.var.list]
    Ko  <- K[-missing.var.list,-missing.var.list]
    Koh <- K[-missing.var.list,missing.var.list]
    Km  <- Ko - Koh %*%solve(Kh)%*% t(Koh)
    cat(paste0("\ndistance Km-Ko = ",mean((Km-Ko)^2)))
    #  K  <- params
    
  }
  if(method=="gen"){
    ##### 1 Genevieve ######
    diag(adjmat) <- (1)
    params <- params*adjmat
    mineig <- 0.1+0.8*runif(1,min=0,max=1)
    K  <- params+(mineig- min(eigen(params, symmetric=TRUE)$values))*diag(p)
  }
  
  if(method=="chris"){
    ##### 2 Christophe ######
    params <- params*adjmat
    degree <- rowSums(abs(params))
    coeff  <- rep(alpha.observed,p)
    
    if (!is.null(missing.var.list))  coeff[missing.var.list]<-  alpha.hidden
    diag(params)<- coeff *  pmax(1,degree)
    params<-eigen(params)$vectors %*% diag( pmax(eigen(params)$values, max(eigen(params)$values)/runif(p,min=5,max=6))) %*% t(eigen(params)$vectors)
    K <-params
    D <-  diag(1/sqrt(diag(K)))
    K <- D%*%K%*%D
    
    if (!is.null(missing.var.list)) {
      Kh  <- K[missing.var.list,missing.var.list]
      Ko  <- K[-missing.var.list,-missing.var.list]
      Koh <- K[-missing.var.list,missing.var.list]
      Km  <- Ko - Koh %*%solve(Kh)%*% t(Koh)
      #print(mean((Km-Ko)^2))
    }
  }
  
  return(K)
}


chooseMissingVar <- function(adjmat,missing.var.number){
  degrees <-colSums(adjmat)
  return(order(degrees,decreasing=TRUE)[1:missing.var.number])
}


getKm<-function(K,nb.missing.var=1){
  p  <-ncol(K)
  Ko <-K[1:(p-nb.missing.var),1:(p-nb.missing.var)]
  Kh <-K[(p-nb.missing.var+1):p,(p-nb.missing.var+1):p]
  Koh<-K[1:(p-nb.missing.var),(p-nb.missing.var+1):p]
  Km <- Ko - Koh %*% solve(Kh) %*% t(Koh)
}

getAm<-function(A,nb.missing.var=1){
  p  <-ncol(A)
  Ao <-A[1:(p-nb.missing.var),1:(p-nb.missing.var)]
  Aoh<-A[1:(p-nb.missing.var),(p-nb.missing.var+1):p]
  Am <- Ao - Aoh %*% t(Aoh)
}



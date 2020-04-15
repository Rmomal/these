library(rcdd)
pivot.fractional <- function(A, trace=FALSE, inverse=TRUE) {
  A<-d2q(as.matrix(A))
  # trace: afficher ou pas la progression
  if (trace){cat("\nPIVOT DE GAUSS POUR ÉCHELONNER UNE MATRICE")}
  n <- nrow(A) # nombre de ligne de A
  if (inverse) {
    A<-cbind(A,diag(rep(1,n)))
    m<-2*n    # nombre de colonnes
  } else m<-n 
  det.factors<-rep("1",n+1) # La derniere case est pour le signe (suivant
  # le nombre de permutation des lignes
  
  for (i in 1:(n)) { # on considère le ième pivot
    if(trace){cat("\n - ITÉRATION",i)}
    
    # VÉRIFICATION DE L'EXISTENCE D'UN PIVOT NON NUL DANS LA COLONE COURANTE
    if (A[i,i] == "0") { # est-ce que le pivot est nul ?
      # Le nouveau pivot candidat est le plus grand élément dans le reste de la colonne
      A.coli<-qabs(A[(i+1):n,i])
      j <- which(qmax(A.coli)==A.coli) + i 
      if (A[j,i] != "0") { # si cet élément est non nul, on peut s'en servir comme pivot
        if (trace) {cat("\n\t+ Échange des lignes",i,j)}
        A[c(i,j),] <- A[c(j,i),] # échange des ligne i et j
        det.factors[n+1]<-qxq(det.factors[n+1],"-1")
      } else { # sinon on n'a pas trouvé de pivot non nul: on s'arrête là
        return(A)
      }
    }
    
    # ÉCHELONNEMENT DE LA COLONNE COURANTE
    # Normalisation de la ligne du pivot 
    det.factors[i]<-A[i,i]
    A[i,]<-qdq(A[i,],rep(A[i,i],m))   # C'est la seule division du programme... 
    
    if (inverse) {set <- setdiff(1:n,i)  # Alors on réduit la matrice
    }  else {set <- setdiff(1:n,1:i) }
    if (trace) {cat("\n\t+ Élimination de la variable",i, "dans les lignes", set)}
    for (j in set) {
      A[j, ] <- qmq(A[j, ] , qxq(rep(A[j,i],m), A[i, ])) # A[i,i]=1
    }
    
    # AFFICHAGE DE L'ÉTAT COURANT DU SYSTÈME
    if (trace) {
      cat("\n\t+ État du système:\n")
      print(A)
    }
  }
  if (inverse) {A.inv<-A[,(n+1):(2*n)]
  } else {A.inv=NULL}
  
  return(list(A=A,
              det.factors=det.factors,
              A.inv=A.inv))
}

inverse.fractional<-function(A){
  A.inv<-pivot.fractional(A,trace=FALSE,inverse=TRUE)$A.inv
  print(A.inv) # A commenter
  q2d(A.inv)
}

det.fractional<-function(A){
  factors<-pivot.fractional(A,trace=FALSE,inverse=FALSE)$det.factors
  det.A<-qprod(factors)
  determinant <-q2d(det.A)
  if (determinant > .Machine$double.xmax)
    determinant <- .Machine$double.xmax 
  if (determinant < .Machine$double.xmin)
    determinant <- .Machine$double.xmin
  return(determinant)}

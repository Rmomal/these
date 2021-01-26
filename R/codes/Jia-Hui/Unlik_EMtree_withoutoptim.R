# 07/09
# Tests to set illegal links
library(EMtree)
library(PLNmodels)
library(tibble)
library(ggplot2)
#----------------
# Function:

# new function :
direct_EMtree<-function(PLN.Cor,unlinked=NULL,n, maxIter=30, cond.tol=1e-10,
                     eps1 = 1e-3,eps2=1e-3, plot=FALSE){
  CorY=cov2cor(PLN.Cor$model_par$Sigma)
  alpha.psi = Psi_alpha(CorY, n, cond.tol=1e-10)
  psi = alpha.psi$psi
  # new beta.init
  beta.init = matrix(1, p, p); diag(beta.init) = 0; 
  if(!is.null(unlinked)) beta.init[unlinked, unlinked]=0
  beta.init = beta.init / sum(beta.init)
  
  # new FitBeta function
  beta.tol = 1e-4
  beta.min = 1e-16
  beta.old = beta.init / sum(beta.init)
  log.psi = log(psi+(psi==0))
  iter = 0
  logpY = rep(0, maxIter)
  beta.diff = diff.loglik = 2 * eps2
  T1<-Sys.time()
  stop=FALSE
  while (((beta.diff > eps1) || (diff.loglik>eps2) ) && iter < maxIter && !stop){
 # while ((beta.diff > eps1) && iter < maxIter && !stop){
    iter = iter+1
     P = EdgeProba(beta.old*psi)
    
    P[beta.old==0]=0
    M = Meila(beta.old) #Meila(beta.old)
    # for(i in 1:10){
  #  lambda = SetLambda(P, M)# ; print(lambda)
    # got rid of optim dependency for the update of beta
    gamma=log(F_Sym2Vec(P)) - log(F_Sym2Vec(M))#+lambda)
    # browser()
    # gamma=gamma-mean(gamma)
   # beta=P/(M+lambda)
    beta=F_Vec2Sym(exp(gamma))
   #  M = Meila(beta) #Meila(beta.old)
   # # beta[which(beta< beta.min)] = beta.min
   #  }
    # P = EdgeProba(beta*psi)
    # P[beta==0]=0
    logpY[iter] = -NegLikelihood(F_Sym2Vec(beta+(beta==0)),log.psi,P)
    beta.diff = max(abs(beta.old-beta))
    beta.old = beta
    diff.loglik=logpY[iter]-logpY[iter-1]
    print(c(beta.diff,diff.loglik))
    if(iter > 1){diff.loglik =  abs(diff.loglik)}else{diff.loglik=1}
  }
  time<-difftime(Sys.time(),T1)
  logpY = logpY[1:iter]
  if(plot){
    g<-tibble(p=logpY) %>% rowid_to_column() %>%
      ggplot(aes(rowid,p))+geom_point()+geom_line()+theme_minimal()+labs(x="Iter",y="Likelihood")
    print(g)
  }
  P = EdgeProba(beta*psi)
  cat("\nConvergence took",round(time,2), attr(time, "units")," and ", iter," iterations.")
  return(list(edges_prob=P, edges_weight=beta, logpY=logpY,maxIter=iter, norm.cst = SumTree(beta), timeEM=time))
}
# also need the SetLambda and NegLikelihood functions, which are not yet exported by EMtree
SetLambda <- function(P, M, eps = 1e-6, start=1){
  # F.x has to be increasing. The target value is 0
  F.x <- function(x){
    if(x!=0){
      1 - sum(P / (x+M))
    }else{
      1 - (2*sum(P[upper.tri(P)] / M[upper.tri(M)]))
    }
  }
  x.min = ifelse(F.x(0) >0,-start,1e-4);
  t1<-Sys.time()
  while(F.x(x.min)>0 ){ 
    x.min = x.min -1 
    if(difftime(Sys.time(),t1)>1) stop("Could not set lambda.")}
  x.max = start
  while(F.x(x.max)<0){x.max = x.max * 2}
  x = (x.max+x.min)/2
  f.min = F.x(x.min)
  f.max = F.x(x.max)
  f = F.x(x)
  while(abs(x.max-x.min) > eps){
    if(f > 0) {
      x.max = x
      f.max = f
    } else{
      x.min = x
      f.min = f
    }
    x = (x.max+x.min)/2;
    f = F.x(x)
  }
  return(x)
}
NegLikelihood <- function(beta.vec, log.psi, P){
  M = Meila(F_Vec2Sym(beta.vec))
 # lambda = SetLambda(P, M)
  return(- sum(F_Sym2Vec(P)*(log(beta.vec)+F_Sym2Vec(log.psi))) +
           log(SumTree(F_Vec2Sym(beta.vec))))#+
           #lambda*(sum(beta.vec)-0.5))
}
F_Vec2Sym <- function(A.vec){
  # Makes a symmetric matrix from the vector made of its lower tirangular part
  n = (1+sqrt(1+8*length(A.vec)))/2
  A.mat = matrix(0, n, n)
  A.mat[lower.tri(A.mat)] = A.vec
  A.mat = A.mat + t(A.mat)
  return(A.mat)
}
F_Sym2Vec <- function(A.mat){
  # Makes a vector from the lower triangular par of a symmetric matrix
  return(A.mat[lower.tri(A.mat)])
}
EdgeProba<-function (W, verbatim = FALSE){
  it = -1
#  logWcum = logSumTree(W)$det
  Wcum=SumTree(W)
  if (!isSymmetric(W)) {
    cat("Pb: W non symmetric!")
  }
  while (!is.finite(Wcum)) {
    it = it + 1
    borne = 30 - it
    if (verbatim)  message(cat("W corrected, bound=", borne))
    W.log = log(F_Sym2Vec(W))
    W.center = W.log - mean(W.log)
    W.center[which(W.center < (-borne))] = -borne
    W = F_Vec2Sym(exp(W.center))
    Wcum = SumTree(W)
  }
  p = nrow(W)
  P = matrix(0, p, p)
 
  sapply(1:(p - 1), function(j) {
    sapply((j + 1):p, function(k) {
      W_jk = W
      W_jk[j, k] = W_jk[k, j] = 0
     # P[k, j] <<- 1 - SumTree(W_jk)/Wcum
      P[k, j] <<- 1 - det(Laplacian(W_jk)[-1,-1])/Wcum
     # P[k, j] <<- 1 - exp(logSumTree(W_jk)$det - logWcum)
      P[j, k] <<- P[k, j]
    })
  })
  P[which(P < 1e-5)] = 1e-5
  diag(P) = 0
  return(P)
}
hist(EdgeProba(newWeight,verbatim=TRUE))

#----------------
# Simulated example:
set.seed(1)
pal_nodes= c("#adc9e0","#e7bd42") ; pal_edges = "#31374f"
n=200; p=15

simu=data_from_scratch("erdos",p=p)
Y=simu$data # count data
G=1*(simu$omega!=0) ; diag(G) = 0 # original dependency graph
draw_network(G, layout="kk", btw_rank=4, nodes_size=c(3, 5), 
             pal_nodes= pal_nodes, pal_edges = pal_edges, nodes_label = 1:p)$G
# Nodes 6, 8 and 15 show the biggests betweenness centrality scores.
unlink=c(6,8,15)

## First the PLN fit:
PLNfit=PLN(Y~1, control = list(trace=0))

## a classic fit of EMtree:
newFit=new_EMtree(PLN.Cor=PLNfit,n=n, unlinked = NULL, plot=TRUE, maxIter = 40,
                    eps1 = 1e-3)
optimFit=optim_EMtree(PLN.Cor=PLNfit,n=n, plot=TRUE,eps1 = 1e-10,maxIter = 8)
newprob=newFit$edges_prob
newprob2=EdgeProba(newFit$edges_weight)
optimProb2=optimFit$edges_prob #quand on enlève M de l'optim
par(mfrow=c(1,1))
plot(F_Sym2Vec(optimProb), F_Sym2Vec(optimProb2))
abline(0,1)
plot(newprob,optimProb)
abline(0,1)
newWeight=newFit$edges_weight
optimWeight=optimFit$edges_weight
plot(newWeight,optimWeight)
abline(0,1)
dataweights=data.frame(optim = F_Sym2Vec(optimWeight), new=F_Sym2Vec(newWeight))
dataweights %>% ggplot(aes(log(optim), log(new)))+geom_point()+theme_light()+
  geom_abline()
hist(log(newWeight), breaks=30)
hist(log(optimWeight), breaks=30)
par(mfrow=c(2,1))
hist((optimProb), breaks=40)
hist((optimProb2), breaks=40)
hist((newprob), breaks=30)
optim_network=1*(optimFit$edges_prob>0.1)
new_network=1*(newFit$edges_prob>0.05)
draw_network(optim_network, layout="kk", groupes =(1:p)%in%unlink, nodes_size=c(3, 5), 
             pal_nodes= pal_nodes, pal_edges = pal_edges,nodes_label = 1:p)$G
draw_network(new_network, layout="kk", groupes =(1:p)%in%unlink, nodes_size=c(3, 5), 
             pal_nodes= pal_nodes, pal_edges = pal_edges,nodes_label = 1:p)$G
auc(pred =optim_network,G )
auc(pred =new_network,G )
## now let's unlink 6, 8 and 15:
unlinkedFit=new_EMtree(PLN.Cor=PLNfit,unlinked = unlink, plot=TRUE)
unlinked_network=1*(unlinkedFit$edges_prob>1e-8)
draw_network(unlinked_network, layout="kk", groupes =(1:p)%in%unlink, nodes_size=c(3, 5), 
             pal_nodes= pal_nodes, pal_edges = pal_edges, nodes_label = 1:p)$G



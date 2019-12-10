#' em.litree
#'
#' The function \code{em.litree} implements the algorithm of Gaussian Graphical Model Inference with missing variable described
#' in Robin et. al (2017). The underlying model is based on the aggregation of spanning trees,
#'  and the estimation procedure on the Expectation-Maximization algorithm.
#'  We treat the graph structure and the unobserved nodes as missing variables and compute posterior probabilities of edge appearance.
#' To provide a complete methodology, we also propose three model selection criteria to estimate the number of missing nodes.
#'
#' @param X X is a data matrix
#' @param k is a vector of missing variable number. It should be an integer greater than 0 or a vector of such integers
#' @param criterion is the name of the criterion used for selecting the best model ({"ICL_T", "ICL_ZT", "BIC" }). ICL_T by default
#' @param max.iter is the maximum number of iterations (20 by default)
#' @param eps is the precision used for stopping the algorithm
#'
#' @return A list of three items, \code{criteria} a dataframe whose columns are the three critera and each lines corresponds
#' to a given number of missing variables, a list of models and the best model according the specified criterion
#'
#' @references  Geneviève Robin, Christophe Ambroise, Stéphane Robin (Submitted on 26 May 2017).
#' Graphical model inference with unobserved variable via latent tree aggregation. Arxiv Paper.
#' \url{https://arxiv.org/abs/1705.09464}
#'
#' @export
#' @import Matrix
#' @importFrom pracma pinv
#' @importFrom igraph plot.igraph
#' @examples
#' data(cyto)
#' res.raf.full <- em.litree(X.raf,0:2)
#'
em.litree<- function(X,k=0:2, criterion = "ICL_T", max.iter = 20,eps = 0.1){
  #########################
  # INPUT PARAMETERS
  #########################
  # X is a data matrix
  # k is a vector of missing variable number. It should be an integer greater than 0 or a vector of such integers
  # criterion, is the name of the criterion used for selecting the best model ({"ICL_T", "ICL_ZT", "BIC" })
  ##########################
  # OUTPUT PARAMETERS
  ##########################
  # A list of three items
  #
  if (any((k-round(k))!=0) | any(k<0)) stop("k is a vector of missing variable number. It should be an integer greater than 0")

  results<- vector(mode = "list",length=length(k))
  initial.param <- vector(mode = "list",length=2)
  # Looping all the components of  vector k
  # ---------------------------------------------------------------------------
  results<-lapply(k,function(nb.missing.var){
   if (nb.missing.var==0) {
     initial.param$Sigma0<-cov(X)
     initial.param$K0<-pinv(initial.param$Sigma0)
   } else {
     initial.param<-initEM(X,cliquelist = findCliques(X,nb.missing.var+1)[1:nb.missing.var])
    }
   em.latent.trees.res<-treeAgr.EM(S = cov(X),k=nb.missing.var, K0 = initial.param$K0,
                                   Sigma0 = initial.param$Sigma0,
                                   pii=0.5, n=nrow(X), max.iter = max.iter,eps = eps)
   criteria<-modelChoiceCriteria(X,em.latent.trees.res)
   # more or less working with nb.missing.var >= 1 (criteria computed but often big
   list(criteria=criteria,em.res=em.latent.trees.res)})

   # Post processing of the models and criteria
   # --------------------------------------------------------------------------
   # Collect all criteria  and models in separated structure
   criteria<- data.frame(do.call("rbind",lapply(results, function(x) (unlist(x$criteria)))))
   models <- lapply(results, function(x) ((x$em.res)))
   # find the best number of missing variables according the chosen criteria
   best.model.index <- switch(criterion,
          ICL_T = which.max(criteria$ICL_T),
          ICL_ZT = which.max(criteria$ICL_ZT),
          BIC = which.max(criteria$BIC),
          stop("Non existing criterion name"))
  return(list( criteria = criteria , models = models, best.model = models[[best.model.index]] ))
}


modelChoiceCriteria <- function(X, EM.res,pii=0.5){
  K=EM.res$K
  Sigma=EM.res$Sigma
  log_gamma=EM.res$gamma
  alpha=EM.res$alpha
  P=EM.res$P
  diag(log_gamma)=0
  n=nrow(X)
  q=nrow(log_gamma)
  p=ncol(X)
  r=nrow(Sigma)-p
  #pii=log(matrix(2/(p+r),p+r,p+r))/n
  pii_mat=log(matrix(pii,p+r,p+r))
  diag(pii_mat)=0
  diag(P)=0
  delta=-exp(log_gamma)
  diag(delta)=rowSums(exp(log_gamma))
  print(det(delta[2:q,2:q]))
  print(q)
  print(p)
  H.T=log(det(delta[2:q,2:q]))-sum(sum(log_gamma*alpha)) # Entropy of trees
  #H.T=sum(log(d))-sum(sum(log_gamma*alpha)) # Entropy of trees
   if (r==0){ # no missing variable
    E.H.Xh=0
     } else
   if (r==1){ # one missing variable
   	E.H.Xh=r*log(2*pi*exp(1))/2-(1/2)*sum(log(K[(p+1),(p+1)])) # Entropy of hidden variables
  } else {
    E.H.Xh=r*log(2*pi*exp(1))/2-(1/2)*sum(log(diag(K[(p+1):(p+r),(p+1):(p+r)]))) # Entropy of hidden variables
  }
  pen.r=(log(n))*(p*(p+1)/2+r*p+r)
  delta0 = -exp(pii_mat)+diag(rowSums(exp(pii_mat)))
  Z0 = det(delta0[2:q,2:q])
  #mat=alpha*(pii-Z0+P- (p+r)*log(2*pi)/2)
  mat=alpha*(pii_mat+P)
  E.c=sum(sum(mat)) - Z0-(p+r)*log(2*pi)/2# Expectation of complete likelihood
  logPY = E.c+E.H.Xh+H.T # Marginal log-likelihood of observed variables
  if (r==0) ICL_ZT <- NA else ICL_ZT <- E.c-pen.r
  return(list(ICL_T=logPY-H.T-pen.r, ICL_ZT= ICL_ZT, BIC=logPY-pen.r))
}




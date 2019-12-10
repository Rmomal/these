#' R6 Class for running experiment of Gaussian Graphical Model
#'
#' The class aims at comparing multiple Gaussian Graphical Model Inference procedure when a ground
#' truth is available
#'
#' @section Usage:
#' \preformatted{experiment = GGMexperiment$new(X.list, adjmat, methods=c("glasso","em.latent.trees"))
#' }
#'
#' @section Arguments:
#' \code{- X.list} A list of data.frame, which will be used as input data for the inference
#'
#' \code{- adjmat} The ground truth adjacency matrix used for the evaluation (roc curve)
#'
#' \code{- nb.missing.var} Number of missing variable (0 by default)
#'
#' \code{- methods} A vector of characters list the methods to be tested (see GGMfit for the list of possible methods)
#'
#' \code{- fit.number} Number of evaluation point (20 by default)
#'
#' \code{- nb.sample}  Number of data frame in X.list (set automatically while initializing the object)
#'
#' \code{- K.score.array}  Array of prediction of edges. The array is 3 dimensional (2 first dimension for storing a results, third dimension for compiling all results)
#'
#' \code{- prediction} A dataframe with 3 columns (prediction, label and method) used for methods evaluation
#'
#'
#' @section Methods:
#' \code{$new(X.list=NULL,adjmat=NULL,nb.missing.var=0,methods="glasso",fit.number=20)} Initialize the experiment
#'
#' \code{$run(bagging=FALSE)} Running the experiment (with or without using bagging)
#'
#' \code{$roc.plot()} Plot the roc curves (all methods on the same plot)
#' @name GGMexperiment
#' @examples
#' \dontrun{
#' star.graph <- graphModel$new(type = "starerdos",size=30, p.or.m = 0.05)
#' star.model <- GGMmodel$new(graph=star.graph)
#' plot(star.model)
#' star.model$randomSample(n=50)
#' testingGlasso<-GGMexperiment$new(X.list = list(star.model$getX()), adjmat = star.model$getAdjmat())
#' testingGlasso$run()
#' print(glasso.auc<-testingGlasso$auc())
#' testingGlasso$roc.plot()
#' }
#'
NULL

#' @importFrom R6 R6Class
#' @import ggplot2
#' @importFrom MASS mvrnorm
#' @importFrom  ROCR   prediction   performance
#' @importFrom bnlearn   chow.liu
#' @import saturnin
#' @import Matrix
#' @importFrom cluster  silhouette
#' @export
GGMexperiment <- R6Class("GGMexperiment",
                          public=list(X.list = NULL,
                                      adjmat =NULL,
                                      nb.missing.var =NULL,
                                      methods = NULL,
                                      fit.number = NULL,  # number of fit (number of different lambdas for glasso or different cutoffs for LITree)
                                      nb.sample =NULL,
                                      K.score.array= NULL,
                                      prediction = NULL,
                                      initialize = function(X.list=NULL,adjmat=NULL,nb.missing.var=0,methods="glasso",fit.number=20){
                                        self$adjmat<-adjmat

                                        if ((is.data.frame(X.list))||(is.matrix(X.list)))  {X.list<-list(X.list)}
                                        p<-ncol(X.list[[1]])
                                        self$X.list<-X.list
                                        self$nb.sample <- length(X.list)
                                        self$methods<-methods
                                        self$K.score.array <-array(0,dim=c(p,p,length(X.list)*length(methods)),dimnames=list(1:p,1:p,rep(methods,self$nb.sample)))
                                        self$nb.missing.var <- nb.missing.var
                                        self$fit.number<-fit.number
                                      },
                                      run = function(bagging=FALSE){
                                        if (is.null(self$X.list)) stop("First initialize the experiment with a list of samples")
                                        run<-0

                                        for (X in self$X.list){
                                          for (method in self$methods)
                                           {
                                            run<-run+1
                                            GGM.fit.with.X<- GGMfit$new(X,method=method,nb.missing.var=self$nb.missing.var,fit.number=self$fit.number)
                                            GGM.fit.with.X$run(bagging=bagging)
                                            self$K.score.array[,,run]<-GGM.fit.with.X$K.score
                                          }}
                                        #p<-dim(self$K.score.array)[1]
                                        nb.run<-dim(self$K.score.array)[3]
                                        methods<-dimnames(self$K.score.array)[[3]]
                                        # ploting the mean ROC curve for each method
                                        K.list.of.vect<-vector("list",nb.run)
                                        names(K.list.of.vect)<-self$methods
                                        for (run in 1:nb.run){
                                          K.list.of.vect[[run]]<-as.vector(removeDiag(self$K.score.array[,,run]))
                                        }

                                        adjmat.vect<- as.vector(removeDiag(self$adjmat))
                                        self$prediction<-Reduce("rbind",Map( function(x,y){
                                          pred <- data.frame(predictions = x, labels = adjmat.vect, method=y)},
                                          K.list.of.vect, methods))

                                      },
                                      save = function(){},
                                      load = function(){},
                                      roc.plot=function(){
                                        # separate the results of each methods in different data.frame
                                        pred.per.method<-split(self$prediction,self$prediction$method)
                                        # Compute the performance of each method
                                        perf<-do.call("rbind", lapply(pred.per.method,function(pred){
                                          method<-pred[1,]$method
                                          perf<- performance(prediction(pred$predictions,pred$labels),'tpr', 'fpr')
                                          perf<- data.frame(TPR= perf@y.values[[1]], FPR = perf@x.values[[1]],method=method)
                                        }))
                                            p <- ggplot(perf, aes(x=FPR, y=TPR))+
                                            geom_point(aes(colour=method,shape =method))+
                                            coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
                                            theme_bw()+
                                            geom_abline(slope = 1, color='grey')+xlab("FPR")+ylab("TPR")
                                          return(p)},
                                      auc = function(){
                                        pred.per.method<-split(self$prediction,self$prediction$method)
                                        auc<-do.call("rbind", lapply(pred.per.method,function(pred){
                                          method<-pred[1,]$method
                                          auc<- performance(prediction(pred$predictions,pred$labels),'auc')@"y.values"[[1]]
                                          auc.line<-data.frame(auc=auc, method=method)
                                        }))
                                        }
                                        ))





#' R6 Class for estimation of Gaussian Graphical Model
#'
#' The Class \code{GGMfit} estimate a Gaussian Graphical Model. It can use 5 different models and associated estimation procedure
#'  with or without bagging.
#'  Three of the algorithms (\code{em.latent.trees}, \code{em.glasso}, \code{em.chow.liu}) take into account missing variables while the two other do not (\code{glasso},
#'  \code{chow.liu}).
#"
#'
#' @section Usage:
#' \preformatted{experiment = GGMfit$new(X,method="glasso",nb.missing.var=0,fit.number=20,...)
#' }
#'
#' @section Arguments:
#' \code{- X} A data.frame, which will be used as input data for the inference
#'
#' \code{- nb.missing.var} Number of missing variable (0 by default)
#'
#' \code{- method} The  name of the  estimation procedure (\code{em.latent.trees}, \code{em.glasso}, \code{em.chow.liu}, \code{glasso},
#'  \code{chow.liu})
#'
#' \code{- fit.number} Number of evaluation point (20 by default)
#'
#'
#' \code{- K.score}  Array of prediction of edges. The array is 2 dimensional homogeneous in dimension to the adjacency matrix of the graph to be inferred.
#'
#'
#'
#' @section Methods:
#' \code{$new(X,method="glasso",nb.missing.var=0,fit.number=20,...)} Initialize the experiment
#'
#' \code{$run(bagging=FALSE,nb.bootstrap=29)} Running the experiment (with or without using bagging)
#'
#' @name GGMfit
#' @examples
#' \dontrun{
#' star.graph <- graphModel$new(type = "starerdos",size=30, p.or.m = 0.05)
#' star.model <- GGMmodel$new(graph=star.graph)
#' star.model.missing <- GGMmodel$new(graph=star.graph,nb.missing.var= 1)
#' dim(star.model.missing$getAdjmat())
#' dim(star.model.missing$getAdjmatCond())
#' star.model.missing$randomSample(n=60)
#' dim(star.model.missing$getX())
#' dim(star.model.missing$getXobs())
#' star.model.missing.fit <-GGMfit$new(star.model.missing$getXobs(),fit.number = 20,method="glasso")
#' star.model.missing.fit$run()
#' star.model.missing.fit2 <-GGMfit$new(star.model.missing$getXobs(),nb.missing.var= 1,fit.number = 20,method="em.latent.trees")
#' star.model.missing.fit2$run()
#' }
NULL

#' @import glasso
#' @import bnlearn
#' @export
GGMfit <- R6Class("GGMfit",
                   public=list(
                     X = NULL,
                     nb.missing.var =NULL,
                     method = NULL,
                     fit.number = NULL,
                     K.score = NULL,
                     initialize= function(X,method="glasso",nb.missing.var=0,fit.number=20,...){
                       self$X<-X
                       self$nb.missing.var <- nb.missing.var
                       self$method <- method
                       self$fit.number <- fit.number
                       p<-ncol(X)
                     },
                     run=function(bagging=FALSE,nb.bootstrap=29){
                       if (bagging==TRUE)
                         self$runBagging(nb.bootstrap=nb.bootstrap)
                       else
                         self$K.score <- private$runSimple(self$X)
                     },
                     runBagging=function(nb.bootstrap){
                        p<-ncol(self$X)
                        X.list<- lapply(vector(mode = "list", length = nb.bootstrap), function(elt) elt<- self$X[sample(1:nrow(self$X),nrow(self$X),replace=TRUE),] )
                        K.score.array<-array(0,dim=c(p,p,length(X.list)),dimnames=list(1:p,1:p,1:length(X.list)))
                        run<-0
                        for (X in X.list){
                            run<-run+1
                            K.score.array[,,run]<-private$runSimple(X)
                        }
                        self$K.score <- apply(K.score.array,c(1,2),mean)
                     }
                  ),
                  private=list(
                    runSimple =function(X){
                      self$nb.missing.var -> nb.missing.var
                      self$fit.number -> fit.number
                      p<-ncol(X)
                      switch(self$method,
                             em.latent.trees={
                               if (nb.missing.var ==0) stop("em.latent.trees works only with a number of missing variables (nb.missing.var >0) ")
                               initial.param<-initEM(X,cliquelist = findCliques(X,nb.missing.var+1)[1:nb.missing.var])
                               em.latent.trees.res<-treeAgr.EM(S = cov(X),k=nb.missing.var, K0 = initial.param$K0, Sigma0 = initial.param$Sigma0, pii=0.5, n=nrow(X), max.iter = 20,eps = 0.1)
                               K.score <- abs(em.latent.trees.res$alpha)
                             },
                             em.chow.liu={
                               if (nb.missing.var ==0) stop("em.chow.liu works only with a number of missing variables (nb.missing.var >0) ")
                               initial.param<-initEM(X,cliquelist = findCliques(X,nb.missing.var+1)[1:nb.missing.var])
                               em.chow.liu.res<-tree.EM(S = cov(X), k=nb.missing.var, initial.param$K0, initial.param$Sigma0, max.iter = 10)
                               K.score <- abs(em.chow.liu.res$K)
                             },
                             glasso={ # Meinhausen and Buhlman method (approx =TRUE) otherwise classical GLASSO
                               if (nb.missing.var !=0) stop("Glasso works only with nb.missing.var=0")
                               S = cov(X)
                               log.lambda.min <- -5
                               log.lambda.max <- log(get.lambda.l1(S))
                               log.lambda <- seq(log.lambda.min, log.lambda.max, length.out = fit.number )
                               MB.res <- lapply(exp(log.lambda), function(lambda) glasso(S, lambda, trace = FALSE, approx = FALSE,
                                                                                         penalize.diagonal = FALSE))
                               adjmat.array <- simplify2array(Map("*",exp(log.lambda),lapply(MB.res, function(x){ (abs(x$wi)>0)*1})))
                               # Let us replace each edge by the  largest Glasso lambda where it disappears (or a sum related to this)
                               K.score <- apply(adjmat.array,c(1,2),sum)
                               # equivalent to the two preceding lines
                               # self$K.score <-  Reduce("+",Map("*",lapply(MB.res, function(x){ (abs(x$wi)>0)*1}),exp(log.lambda)))
                             },
                             em.glasso = {
                               if (nb.missing.var ==0) stop("em.glasso works only with a number of missing variables (nb.missing.var >0) ")
                               initial.param<-initEM(X,cliquelist = findCliques(X,nb.missing.var+1)[1:nb.missing.var])
                               log.lambda.min <- -3
                               log.lambda.max <- log(get.lambda.l1(initial.param$Sigma0))
                               log.lambda <- seq(log.lambda.min, log.lambda.max, length.out =  fit.number)
                               EM.Glasso.res <- lapply(exp(log.lambda), function(lambda) convex.EM(cov(X), k=nb.missing.var, lambda, initial.param$K0, initial.param$Sigma0, max.iter = 10, eps=1e-3))
                               adjmat.list <- lapply(EM.Glasso.res, function(res) res$adjmat)
                               K.score <- Reduce("+",Map("*",adjmat.list,exp(log.lambda)))

                             },
                             chow.liu={
                               if (nb.missing.var !=0) stop("Chow-Liu works only with nb.missing.var=0")
                               K.score <- ChowLiu2(S=cov(X))
                             },
                             # recursive.grouping={
                             #   if (nb.missing.var !=1) stop("Recursive grouping works only with nb.missing.var=1")
                             #   S<-cov(X)
                             #   D <- diag(1/sqrt(diag(S)))
                             #   Cor <- D%*%S%*%D
                             #   p <- nrow(S)
                             #   tiny  <-  1e-10
                             #   Dist <- -log(abs(Cor)+tiny)
                             #   Dist <- (Dist>0)*Dist
                             #   diag(Dist) <- 0
                             #   self$K.score <- RG2(Dist, n=nrow(X))$adjmat
                             # },
                             stop("Method not available !")
                      )
                      return(K.score / max(K.score))

                    }
                  ))


#' R6 Class for simulating Gaussian Graphical Model
#'
#' The Class \code{GGMmodel} simulate Gaussian Graphical Model. It can use many different models
#"
#'
#' @section Usage:
#' \preformatted{experiment = GGMmodel$new(graph=NULL, prop.positive.cor=1, type="erdos",size=30, p.or.m =0.1,eta=0.2,extraeta=eta/5,nb.missing.var=0,alpha.hidden= 2,alpha.observed = 1.2)
#' }
#'
#' @section Arguments:
#' \code{- prop.positive.cor} A real number indicating the proportion of positive correlation (1 by default)
#'
#' \code{- graph} A graph object
#'
#' \code{- K}  K  Precision matrix derived from the Graph
#'
#' \code{- Sigma} Covariance matrix derived from the Graph
#'
#' \code{- missing.var.list}   indices of the missing variables if any
#'
#' \code{- X} Simulated data using a Gaussian model with zero mean vector and covariance matrix of the object
#'
#'
#' @section Methods:
#' \code{$new(graph=NULL, prop.positive.cor=1, type="erdos",size=30, p.or.m =0.1,eta=0.2,extraeta=eta/5,nb.missing.var=0,alpha.hidden= 2,alpha.observed = 1.2)}
#' Initialize the model
#'
#' \code{$getAdjmat()} returns the adjacency matrix of the graph,
#'
#' \code{$getAdjmatCond()} returns the conditional adjacency matrix if there is missing data
#'
#' \code{$getAdjmatMarg=()} returns the marginal adjacency matrix if there is missing data
#'
#'  \code{$randomSample(n=100)}  generates a random sample of size n accissible through the argument X
#'
#'  \code{$getX()} access the generated Data
#'
#'  \code{$getXobs()} access the observed part of the generated Data
#'
#'  \code{getXmis()} access the missing part of the generated Data
#'
#'  \code{plot()} plot the generated Data matrix
#'
#' @examples
#' \dontrun{
#' star.graph <- graphModel$new(type = "starerdos",size=30, p.or.m = 0.05)
#' star.model <- GGMmodel$new(graph=star.graph)
#' plot(star.model)
#' star.model$randomSample(n=50)
##' }
#' @name GGMmodel
NULL

#' @export
GGMmodel  <-  R6Class("GGMmodel",
                       public = list(
                         prop.positive.cor = NULL,  # proporition of positive correlation in the correlation matrix
                         graph=NULL,  #   graph : a graph.model object with an adjancency matrix and simulation parameters
                         K=NULL,   #   K  : precision matrix
                         Sigma=NULL,   #   Sigma  : covariance matrix
                         missing.var.list  = NULL,  # indices of the missing variables
                         X=NULL,       # Data
                         initialize=function(graph=NULL, prop.positive.cor=1, type="erdos",size=30, p.or.m =0.1,eta=0.2,extraeta=eta/5,nb.missing.var=0,alpha.hidden= 2,alpha.observed = 1.2)
                         {
                           self$prop.positive.cor = prop.positive.cor
                           self$missing.var.list=NULL

                           if (!is.null(graph)){
                             if (class(graph)[1]!="graphModel") stop("graph object should be of class graphModel")
                             self$graph=graph
                           } else {
                             self$graph = graph.model$new(type=type,size=size, p.or.m =p.or.m,eta=eta,extraeta=extraeta)
                           }
                           # Dealing with missing variable: choosing the nodes of highest degree
                           if (nb.missing.var > 0) {
                             self$missing.var.list = chooseMissingVar(self$graph$adjmat,nb.missing.var)
                           # Reorder the variables  to have the missing variables as the last variables (for evaluation when 1 var is missing)
                           new.order<-c(c(1:size)[-self$missing.var.list] , self$missing.var.list)
                           self$missing.var.list<-(size - nb.missing.var+1):size
                           self$graph$adjmat<-self$graph$adjmat[new.order,new.order]
                           }
                           # Finding a graph compatible covariance matrix taking missing variable into account
                           self$K <- covarianceFromGraph(adjmat =self$graph$adjmat, prop.positive.cor=prop.positive.cor, missing.var.list=self$missing.var.list,alpha.hidden= alpha.hidden,alpha.observed =alpha.observed)
                           self$Sigma=solve(self$K)

                           # Transformation into a correlation matrix
                           D    <-  diag(1/sqrt(diag(self$Sigma)))
                           self$Sigma  <-  D%*%self$Sigma%*%D

                         },
                         getAdjmat=function(){return(self$graph$adjmat)},
                         getAdjmatCond=function(){if (is.null(self$missing.var.list)) stop("No missing data")
                                                  return(self$graph$adjmat[-self$missing.var.list,-self$missing.var.list])
                                                      },
                         getAdjmatMarg=function(){return(getAm(self$graph$adjmat,length(self$missing.var.list)))},
                         randomSample=function(n=100){self$X <- mvrnorm(n, mu=rep(0,nrow(self$Sigma)), self$Sigma) },
                         getX=function(){if (is.null(self$X)) stop("No simulated data available") else return(self$X) },
                         getXobs=function(){
                           if (is.null(self$missing.var.list)) stop("No missing data")
                           if (is.null(self$X)) stop("No simulated data available") else return(self$X[,-self$missing.var.list]) },
                         getXmis=function(){
                           if (is.null(self$missing.var.list)) stop("No missing data")
                           if (is.null(self$X)) stop("No simulated data available") else return(self$X[, self$missing.var.list]) },
                         plot=function(...){image(self$K, ...)}
                       )
)

#' R6 Class for simulation of Graphs
#'
#'
#' @section Usage:
#' \preformatted{experiment = graphModel$new(type="erdos",size=30, p.or.m =0.1,eta=0.2,extraeta=eta/5)}
#'
#'
#' @section Arguments:
#'
#' \code{- type="erdos"} type of simulated graph ('erdos', 'starerdos', 'star', 'tree')
#'
#' \code{- size=NULL}  number of nodes (variables)
#'
#' \code{- p.or.m=NULL}  probability of an edge (if smaller than 1) or number of edges (if greater than 1)
#'
#' \code{eta=NULL}       factor affected to ...
#'
#' \code{extraeta=NULL} factor affected  to ...
#'
#' \code{adjmat=NULL}  simulated adjacency matrix
#'
#'
#' @section Methods:
#'
#' \code{$new(type="erdos",size=30, p.or.m =0.1,eta=0.2,extraeta=eta/5)}  creates a random graph
#'
#' \code{$plot()} plot the graph using the igraph plot function
#'
#' @examples
#' \dontrun{
#' star.graph <- graphModel$new(type = "starerdos",size=30, p.or.m = 0.05)
#' plot(star.graph)
##' }
#' @name graphModel
NULL

#' @import igraph
#' @export
graphModel <- R6Class(classname = "graphModel",
                       public= list(type="erdos",
                                    size=NULL,
                                    p.or.m=NULL,
                                    eta=NULL,
                                    extraeta=NULL,
                                    adjmat=NULL,
                                    initialize = function(type="erdos",size=30, p.or.m =0.1,eta=0.2,extraeta=eta/5){
                                      self$type=type
                                      self$size=size
                                      switch(type,
                                             sparse={},
                                             dense={self$p.or.m=1},
                                             GGMselect={self$eta=eta; self$extraeta=extraeta},
                                             #simone={self$eta=eta; self$extraeta=extraeta},
                                             star={},
                                             doublestar={},
                                             hmm ={},
                                             regular={},
                                             erdos={self$p.or.m=p.or.m},
                                             starerdos={self$p.or.m=p.or.m},
                                             complete={self$p.or.m=1},
                                             tree={},
                                             stop("Graph type not recognized !")
                                      )
                                      self$adjmat = buildGraph(graph.type=self$type, p=self$size, eta=self$eta,extraeta=self$extraeta, p.or.m=self$p.or.m, erdos.type="gnp")},
                                    plot=function(...){ plot.igraph(graph_from_adjacency_matrix(self$adjmat),...)}
                       )
)


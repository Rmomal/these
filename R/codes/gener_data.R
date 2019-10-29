

##############
# DATA
##############

generator_graph<-function(d = 20, graph = "tree", g = NULL, prob = NULL, dens=0.3, vis = FALSE,
                          verbose = TRUE,r=5){
  gcinfo(FALSE)
  if (verbose)
    #cat("Generating data from the multivariate normal distribution with the",
    #     graph, "graph structure....")
    if (is.null(g)) {
      g = 1
      if (graph == "hub" || graph == "cluster") {
        if (d > 40)
          g = ceiling(d/20)
        if (d <= 40)
          g = 2
      }
    }
  theta = matrix(0, d, d)
  if (graph == "cluster") {

    theta<-SimCluster(d,3,dens,r) #prob=5/d ?

  }
  if (graph == "scale-free") {

    out = .C("SFGen", dd0 = as.integer(2), dd = as.integer(d),
             G = as.integer(theta), PACKAGE = "huge")
    theta = matrix(as.numeric(out$G), d, d)
  }
  if(graph=="tree"){
    theta<-SpannTree(d)
  }
  if(graph=="erdos"){
    theta<- erdos(d=d,p=prob)
    if(sum(theta)<4){
      while(sum(theta)<4){
        theta<- erdos(d=d,p=prob)
      }
    }
  }

  #browser()
  if (verbose)
    cat("done.\n")
  rm(vis, verbose)
  gc()
  return(theta = Matrix(theta, sparse = TRUE))
}
generator_param<-function(G,v=1){
  cste = 1
  omega = diag(rep(cste, ncol(G))) + G*v
  while (min(eigen(omega)$values) < 1e-6){
    cste = 1.1*cste
    omega = diag(rep(cste, ncol(G))) + G*v
  }
  #browser()
  sigma = cov2cor(solve(omega))
  #omega = solve(sigma)
  sim=list(sigma=sigma,omega=omega,cste=cste)
  return(sim)
}
generator_PLN<-function(Sigma,covariates, n=NULL){
  # vraies abondances, log normales
  p<-ncol(Sigma) # nb esp??ces
  if(!is.null(covariates)){
    n<-nrow(covariates)
    c<-ncol(covariates)# nb covariables

    m<- model.matrix(~X1+X2+X3,covariates)[,-1]
    mc<-ncol(m)
    beta<-matrix(runif(p*mc),mc,p)
    prod=m %*% beta
  }else{
    prod=0
  }

  Z<- rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
  Y = matrix(rpois(n*p, exp(Z+prod )), n, p)
  return(list(Y,cor(Z)))
}
#data_from_scratch: generate data under the PLN model with a certain type of dependency structure,
#                   and draw the latter structure if asked.
# type: type of graph, either "tree", "erdos", "cluster" or "scale-free"
# p: wanted number of columns (species)
# covar: covariates
# prob: edge probability for erdos graphs
# dens: density of edges for cluster graphs
# r: within/between connectiviy ratio for cluster graphs
data_from_scratch<-function(type, p=20, r=5, covar,prob=log(p)/p,dens=log(p)/p, draw=FALSE){
  # make graph
  graph<- generator_graph(graph=type,d=p,prob=prob,dens=dens,r=r)
  param<-generator_param(as.matrix(graph))
  data<-generator_PLN(param$sigma,covar)[[1]]
  if(draw){ as_tbl_graph(as.matrix(graph)) %>%
      ggraph(layout="kk")+
      geom_edge_link()+
      geom_node_point(size=3, color="blue")
  }
  return(list(data=data,omega= param$omega))
}

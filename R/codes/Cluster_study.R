##############
# preliminaries
#############
library(tidyverse)
library(gridExtra)
library(grid)
library(sna)
# some useful functions
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
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right"),title) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight),top=title),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth),top=title))
  
  grid.newpage()
  grid.draw(combined)
  # return gtable invisibly
  invisible(combined)
}

################
# studied functions
##############

# << !!! code de Stephane !!!
SimCluster <- function(p, k, d, r){
  # k = nb clusters, d = graph density, r = within/between ratio connection probability
  beta = d / (r / k + (k - 1) / k)
  alpha = r * beta
  while (alpha > 1) {
    r = .9 * r
    beta = d / (r / k + (k - 1) / k)
    alpha = r * beta
  }
  Z = t(rmultinom(p, 1, rep(1 / k, k)))
  ZZ = Z %*% t(Z)
  diag(ZZ) = 0
  ZZ = F_Sym2Vec(ZZ)
  G = F_Vec2Sym(rbinom(p * (p - 1) / 2, 1, alpha * ZZ + beta * (1 - ZZ)))
  gplot(G, gmode='graph', label=1:p, vertex.col=Z%*%(1:k),
        main=paste0(round(r,2),"// alpha ",round(alpha,2),"// beta ",round(beta,2)))
  
  return(G)
}
# !!! code de Stephane !!! >>
#see values of r(output) alpha and beta for input values of d and r
varyd<-function(d,r,k=3){
  beta = d / (r / k + (k - 1) / k)
  alpha = r * beta
  while (alpha > 1) {
    r = .9 * r
    beta = d / (r / k + (k - 1) / k)
    alpha = r * beta
  }
  return(list(r=r,alpha=alpha,beta=beta))
}

varyr<-function(d,r,k=3){
  beta = d / (r / k + (k - 1) / k)
  alpha = r * beta
  while (alpha > 1) {
    d = .9 * d
    beta = d / (r / k + (k - 1) / k)
    alpha = r * beta
  }
  return(list(d=d,alpha=alpha,beta=beta))
}

# vary  d and r, and create a nice graph
graph.r<-function(r){
  #seqd<-seq(1 / p, 8 / p, 0.05)
  seqd<-seq(0,1,0.1)
  # param<-data.frame(sapply(seqd, function(y)
  #   sapply(r, function(x)
  #     varyd(d = y, x))))
  param<-data.frame(sapply(r, function(y)
    sapply(seqd, function(x)
      varyr(r = y, d=x))))
  #colnames(param)<-seqd
  colnames(param)<-r
#listparam<-c("r","alpha","beta")
  listparam<-c("d","alpha","beta")
  param$par<-rep(listparam,length(seqd))  
  for (par in listparam){
    # param_r<-(param[param$par==par,]) %>% 
    #   gather(d,val,-par)
    param_r<-(param[param$par==par,]) %>% 
      gather(r,val,-par)
    param_r$val<-as.numeric(as.character(unlist(param_r$val)))
    param_r$r<-as.numeric(as.character(param_r$r))
    param_r$d<-seqd
    #title<-ifelse(par=="r","effective r",par)
    title<-ifelse(par=="d","effective d",par)
    assign(paste0("graph_",par),ggplot(param_r,aes(r,val,color=d))+
             geom_point()+
             labs(y=title))+
      theme(legend.position="none")
  }
  grid_arrange_shared_legend(graph_alpha,
                           graph_beta,graph_d,nrow=3,ncol=1,title="",position="right")
}


#######
# RUN
#######
# I. Cluster viz
p=10
par(mfrow=c(2,2))
SimCluster(p,3,5/p,2)
SimCluster(p,3,5/p,5)
SimCluster(p,3,5/p,9)
SimCluster(p,3,5/p,15)

# II. parameters viz
r<-1:30
p<-10
graph.r(r)


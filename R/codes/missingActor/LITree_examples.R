devtools::install_github("cambroise/LITree")
# bnlearn, gtable, lazyeval, rlang, tibble, cli, assertthat, ggplot2, igraph, R6, Rcpp,
# RcppEigen, saturnin, simone, blockmodels, RcppArmadillo
# no package pracma

library(LITree)
# simulation of graph (ground truth)
star.graph <- graphModel$new(type = "starerdos",size=30, p.or.m = 0.05)
plot(star.graph, vertex.color="orange",vertex.frame.color="orange", edge.color="steelblue4",
     edge.arrow.size=0.3,edge.arrow.width=1)

# simulation of a GGM
star.model <- GGMmodel$new(graph=star.graph)
plot(star.model)
star.model$randomSample(n=50) # generate random sample
star.model$getX() #access the sample

# Estimate of a GGM, with or without missing variables
star.model.missing <- GGMmodel$new(graph=star.graph,nb.missing.var= 1)
dim(star.model.missing$getAdjmat())
dim(star.model.missing$getAdjmatCond())
star.model.missing$randomSample(n=60)
dim(star.model.missing$getX())
dim(star.model.missing$getXobs())
star.model.missing.fit <-GGMfit$new(star.model.missing$getXobs(),fit.number = 20,method="glasso")
# possible methods: "em.latent.trees", "em.glasso", "em.chow.liu", "glasso", "chow.liu"
star.model.missing.fit$run()
star.model.missing.fit2 <-GGMfit$new(star.model.missing$getXobs(),nb.missing.var= 1,fit.number = 20,
                                     method="em.latent.trees")
star.model.missing.fit2$run()
star.model.missing.fit2$
# comment comparer les deux modèles

# compare several GGM inferences when a ground thruth is available
testingGlasso<-GGMexperiment$new(X.list = list(star.model$getX()), adjmat = star.model$getAdjmat())
testingGlasso$run()
print(glasso.auc<-testingGlasso$auc())
testingGlasso$roc.plot()

X.list<- lapply(vector(mode = "list", length = 10), function(x) {star.model$randomSample(30)
  star.model$getX()  } )
testingGlasso<-GGMexperiment$new(X.list = X.list, adjmat = star.model$getAdjmat(),methods=c("glasso","chow.liu"))
testingGlasso$run()
testingGlasso$roc.plot()
print(auc<-testingGlasso$auc())

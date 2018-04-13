G<-generator_graph(d = 7, graph = "erdos", g = NULL, prob = 0.2)
net<-net_from_matrix(G,1e-16,FALSE)
plot(net)
lambda = 1
omega = diag(rep(lambda, ncol(G))) + G
while (min(eigen(omega)$values) < 0){
  lambda = 1.1*lambda
  omega = diag(rep(lambda, ncol(G))) + G
}

sigma = cov2cor(solve(omega))

sigma<-param_18$sigma
X<-generator_data(n=200,sigma)
#saveRDS(X,"failexample.rds")
inf_treeggm<-TreeGGM(X,"FALSE")$P
infS<-TreeSteph(X)
heatmap(as.matrix(omega),Rowv=NA,Colv=NA)
heatmap(inf_treeggm,Rowv=NA,Colv=NA)

hist(inf_treeggm)
hist(infS$P)

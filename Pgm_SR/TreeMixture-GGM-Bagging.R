# Edge probability in a mixture of trees : GGM
# beta = coef for 'prior' edge probabilities : p(T) = prod_{ij \in T} beta_ij / B
# phi = coef for conditional distribution of : P(Y | T) = prod_{ij \in T} phi_ij

rm(list=ls()); par(pch=20, mfrow=c(2, 2), mex=3/4); 
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('/home/robin/PgmR/Network/FunctionsTree.R')
source('Functions/FunctionsInference.R')
source('Functions/FunctionsSimul.R')
source('Functions/FunctionsDivers.R')
library(sna); library(igraph); library(mvtnorm); library(ROCR); library(glasso); library(vegan); library(mclust)

# Parms
p = 20; n = 100; 
Type = 'Cluster'; Connex = T
r = 5; # ratio within/between for clusters models
B = 50; BB = 0
P = p*(p-1)/2
SimDir = paste0('../Simu/', Type, '-connex', Connex, '-p', p, '-n', n, '/')

# # Data
# Y = as.matrix(read.table("../Data/RAF/1. cd3cd28.csv", h=T, sep=','))

d.list = c(1, 2, 3, 4, 6, 8, 10)/p; d.nb = length(d.list)
# d.list = c(8, 10)/p; d.nb = length(d.list)
# d.list = c(1, 3)/p; d.nb = length(d.list)
# d.list = c(2)/p; d.nb = length(d.list)

for(dd in 1:d.nb){
   d = d.list[dd]
   
   Perf.1step = Perf.EM = Perf.GL = Perf.EM.bagging = rep(0, B)
   Pred.EM.bagging = matrix(0, P, B)
   Cum.Pred.1step = Cum.Pred.EM = Cum.Pred.GL = Cum.Pred.EM.bagging = matrix(0, B*P, 2)
   P.1step = P.EM = matrix(0, B, P)
   par(mfrow=c(6, 4), mex=.5)
   for (b in 1:B){
      # b = 1
      cat(b, '')
      # Graph 
      if(Type == 'Tree'){G = SimTree(p); G.igraph = graph.adjacency(G, mode='undirected')}
      if(Type == 'Erdos'){
         G = SimErdos(p, d); G.igraph = graph.adjacency(G, mode='undirected')
         if(Connex){while (clusters(G.igraph)$no > 1){G = SimErdos(p, d); G.igraph = graph.adjacency(G, mode='undirected')}}
         if(!Connex){while (clusters(G.igraph)$no > sqrt(p)){G = SimErdos(p, d); G.igraph = graph.adjacency(G, mode='undirected')}}
      }
      if(Type == 'Cluster'){
         G = SimCluster(p, 3, d, r); G.igraph = graph.adjacency(G, mode='undirected')
         if(Connex){while (clusters(G.igraph)$no > 1){G = SimCluster(p, 3, d, r); G.igraph = graph.adjacency(G, mode='undirected')}}
         if(!Connex){while (clusters(G.igraph)$no > sqrt(p)){G = SimCluster(p, 3, d, r); G.igraph = graph.adjacency(G, mode='undirected')}}
      }
      nb.triangles = sum(adjacent.triangles(G.igraph))/3
      # Position = gplot(G, gmode='graph')
      
      # Parms
      lambda = 1.1; Omega = diag(rep(lambda, p)) + G; 
      while (min(eigen(Omega)$values) < 0){lambda = 1.1*lambda; Omega = diag(rep(lambda, p)) + G}
      Sigma = solve(Omega)
      beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)
      P.unif = Kirshner(beta.unif)$P
      
      # Sim data
      Y = rmvnorm(n, sigma=Sigma); Cor = cor(Y); Cov = cov(Y)
      phi = 1/sqrt(1 - Cor^2); diag(phi) = 0
      # phi = (1 - Cor^2)^(-n/2); diag(phi) = 0
      
      # 1 step
      PedgeY.1step = Kirshner(beta.unif*phi)$P
      Cum.Pred.1step[((b-1)*P+1):(b*P), ] = cbind(F_Sym2Vec(PedgeY.1step), F_Sym2Vec(G))
      # boxplot(F_Sym2Vec(PedgeY.1step) ~ F_Sym2Vec(G), ylim=c(0, 1))
      P.1step[b, ] = F_Sym2Vec(PedgeY.1step)
      
      # EM
      FitEM = FitBetaStatic(beta.init=beta.unif, phi=phi)
      # plot(FitEM$logpY)
      PedgeY.EM = Kirshner(FitEM$beta*phi)$P
      Cum.Pred.EM[((b-1)*P+1):(b*P), ] = cbind(F_Sym2Vec(PedgeY.EM), F_Sym2Vec(G))
      P.EM[b, ] = F_Sym2Vec((PedgeY.EM))
      # plot(F_Sym2Vec(PedgeY1), F_Sym2Vec(PedgeY2), col=1+F_Sym2Vec(G), xlab='', ylab='')
      Pred.1step = prediction(F_Sym2Vec(PedgeY.1step), F_Sym2Vec(G)); 
      Perf.1step[b] = performance(Pred.1step, measure='auc')@y.values[[1]][1]
      Pred.EM = prediction(F_Sym2Vec(PedgeY.EM), F_Sym2Vec(G)); 
      Perf.EM[b] = performance(Pred.EM, measure='auc')@y.values[[1]][1]
      
      # Glasso
      FitGL = FitGlasso(cov(Y))
      Pred.GL = prediction(F_Sym2Vec(abs(FitGL)), F_Sym2Vec(G)); 
      Cum.Pred.GL[((b-1)*P+1):(b*P), ] = cbind(F_Sym2Vec(FitGL), F_Sym2Vec(G))
      Perf.GL[b] = performance(Pred.GL, measure='auc')@y.values[[1]][1]
      
      # Thresholding lambda
      ScoreGL = F_Sym2Vec(abs(FitGL)); if(sum(is.na(ScoreGL))>0){ScoreGL = ScoreGL[-which(is.na(ScoreGL))]}
      MC = Mclust(log10(ScoreGL), G=2, modelsNames='E')
      thres.GL = ScoreGL[which.min(abs(MC$z[, 1]-.5))]
      
      # Edge prediction
      boxplot(log10(F_Sym2Vec(PedgeY.1step)) ~ F_Sym2Vec(G)); abline(h=log10(2/p))
      hist(log10(F_Sym2Vec(PedgeY.1step)), breaks=p, main='', xlab='', ylab=''); abline(v=log10(2/p))
      boxplot(log10(F_Sym2Vec(PedgeY.EM)) ~ F_Sym2Vec(G), col=2); abline(h=log10(2/p))
      hist(log10(F_Sym2Vec(PedgeY.EM)), breaks=p, main='', xlab='', ylab='', col=2); abline(v=log10(2/p))
      boxplot(log10(F_Sym2Vec(abs(FitGL))) ~ F_Sym2Vec(G), col=4); abline(h=log10(thres.GL))
      hist(log10(F_Sym2Vec(abs(FitGL))), breaks=p, main='', xlab='', ylab='', col=4); abline(v=log10(thres.GL))
      
      # # Bagging EM
      # PedgeY.EM.bagging = matrix(0, p, p); # cat('(')
      # for (bb in 1:BB){
      #    # cat(bb, '')
      #    # sample = sample(1:n, round(n/2))
      #    sample = sample(1:n, n, replace=T)
      #    Y.sample = Y[sample, ]; Cor.sample = cor(Y.sample);
      #    phi.sample = 1/sqrt(1 - Cor.sample^2); diag(phi.sample) = 0
      #    FitEM.sample = FitBetaStatic(beta.init=beta.unif, phi=phi.sample)
      #    PedgeY.EM.sample = Kirshner(FitEM.sample$beta*phi.sample)$P
      #    PedgeY.EM.bagging = PedgeY.EM.bagging + PedgeY.EM.sample
      # }
      # # cat(') ')
      # PedgeY.EM.bagging = PedgeY.EM.bagging / B
      # boxplot(F_Sym2Vec(PedgeY.EM.bagging) ~ F_Sym2Vec(G), ylim=c(0, 1), col=3)
      # Cum.Pred.EM.bagging[((b-1)*P+1):(b*P), ] = cbind(F_Sym2Vec(PedgeY.EM.bagging), F_Sym2Vec(G))
      # Pred.EM.bagging = prediction(F_Sym2Vec(PedgeY.EM.bagging), F_Sym2Vec(G));
      # Perf.EM.bagging[b] = performance(Pred.EM.bagging, measure='auc')@y.values[[1]][1]
      
      # ROC
      Perf.1step.ROC = performance(Pred.1step, 'tpr', 'fpr')
      plot(Perf.1step.ROC, col=1, xlab='', ylab=''); cutoff.1step = which.min(abs(Perf.1step.ROC@alpha.values[[1]]-2/p)); 
      points(Perf.1step.ROC@x.values[[1]][cutoff.1step], Perf.1step.ROC@y.values[[1]][cutoff.1step], col=1)
      Perf.EM.ROC = performance(Pred.EM, 'tpr', 'fpr')
      plot(Perf.EM.ROC, col=2, add=T); cutoff.EM = which.min(abs(Perf.EM.ROC@alpha.values[[1]]-2/p)); 
      points(Perf.EM.ROC@x.values[[1]][cutoff.EM], Perf.EM.ROC@y.values[[1]][cutoff.EM], col=2)
      Perf.GL.ROC = performance(Pred.GL, 'tpr', 'fpr')
      plot(Perf.GL.ROC, col=4, add=T); cutoff.GL = which.min(abs(Perf.GL.ROC@alpha.values[[1]]-thres.GL)); 
      points(Perf.GL.ROC@x.values[[1]][cutoff.GL], Perf.GL.ROC@y.values[[1]][cutoff.GL], col=4)
      # plot(performance(Pred.EM.bagging, 'tpr', 'fpr'), col=3, add=T)
      
      # PR
      Perf.1step.PR = performance(Pred.1step, 'prec', 'rec'); plot(Perf.1step.PR, col=1, xlab='', ylab='')
      points(Perf.1step.PR@x.values[[1]][cutoff.1step], Perf.1step.PR@y.values[[1]][cutoff.1step], col=1)
      Perf.EM.PR = performance(Pred.EM, 'prec', 'rec'); plot(Perf.EM.PR, , col=2, add=T)
      points(Perf.EM.PR@x.values[[1]][cutoff.EM], Perf.EM.PR@y.values[[1]][cutoff.EM], col=2)
      Perf.GL.PR = performance(Pred.GL, 'prec', 'rec'); plot(Perf.GL.PR, , col=4, add=T)
      points(Perf.GL.PR@x.values[[1]][cutoff.GL], Perf.GL.PR@y.values[[1]][cutoff.GL], col=4)
      # plot(performance(Pred.EM.bagging, 'prec', 'rec'), col=3, add=T)
      
      # Predicted graph with EM (thres = 2/p)
      CompGraph(G, 1*(PedgeY.EM>2/p))
   }
   SimResName = paste0(SimDir, 'TreeMixture-GGM-', Type, '-p', p, '-d_', round(100*d), '-n', n, '.Rdata')
   save(P.1step, P.EM, Perf.1step, Perf.EM, Perf.GL, file=SimResName)
   
   
   pdf(paste0(SimDir, 'TreeMixture-GGM-', Type, '-p', p, '-d_', round(100*d), '-n', n, '.pdf'))
   par(mfrow=c(2, 3))
   gplot(G, gmode='graph', main=paste0('d=', round(sum(G)/2/P, 3), ' t=', nb.triangles))
   hist(colSums(G), breaks = (-1:max(colSums(G)))+.5)
   # boxplot(cbind(Perf.1step, Perf.EM, Perf.GL, Perf.EM.bagging), log='y', col=c(0, 2, 4, 3))
   boxplot(cbind(Perf.1step, Perf.EM, Perf.GL), log='y', col=c(0, 2, 4, 3))
   plot(as.vector(P.1step), as.vector(P.EM), log='xy', col=1+as.vector(rep(1, B)%o%F_Sym2Vec(G)))
   plot(performance(prediction(Cum.Pred.1step[, 1], Cum.Pred.1step[, 2]), 'tpr', 'fpr'))
   plot(performance(prediction(Cum.Pred.EM[, 1], Cum.Pred.EM[, 2]), 'tpr', 'fpr'), add=T, col=2)
   plot(performance(prediction(Cum.Pred.GL[, 1], Cum.Pred.GL[, 2]), 'tpr', 'fpr'), add=T, col=4)
   # plot(performance(prediction(Cum.Pred.EM.bagging[, 1], Cum.Pred.EM[, 2]), 'tpr', 'fpr'), add=T, col=3)
   plot(performance(prediction(Cum.Pred.1step[, 1], Cum.Pred.1step[, 2]), 'prec', 'rec'))
   plot(performance(prediction(Cum.Pred.EM[, 1], Cum.Pred.EM[, 2]), 'prec', 'rec'), add=T, col=2)
   plot(performance(prediction(Cum.Pred.GL[, 1], Cum.Pred.GL[, 2]), 'prec', 'rec'), add=T, col=4)
   # plot(performance(prediction(Cum.Pred.EM.bagging[, 1], Cum.Pred.EM.bagging[, 2]), 'prec', 'rec'), add=T, col=3)
   dev.off()
}

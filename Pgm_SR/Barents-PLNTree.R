# PLNtree for Barents fish data

rm(list=ls()); par(pch=20)
library(PLNmodels); library(sna)
source('../pack1/R/codes/FunctionsMatVec.R')
source('../pack1/R/codes/FunctionsTree.R')
source('../pack1/R/codes/FunctionsInference.R')

# Data
data.dir = '../Data_SR/'
data.name = 'BarentsFish'
data.name = 'BarentsFish_Group'
load(paste0(data.dir, data.name, '.Rdata'))
n = nrow(Data$count); p = ncol(Data$count); species = colnames(Data$count)
freq.sel.thres = 0.80

# Functions
TreeGGM <- function(CorY){
   p = ncol(CorY);
   psi = 1/sqrt(1 - CorY^2); diag(psi) = 0
   beta.unif = matrix(1, p, p); 
   FitEM = FitBetaStatic(Y, beta.init=beta.unif, phi=psi)
   Pedge.marg = Kirshner(FitEM$beta)$P
   Pedge.cond = Kirshner(FitEM$beta*psi)$P
   return(list(beta=FitEM$beta, log.pY=FitEM$logpY, psi=psi, Pedge.marg=Pedge.marg, Pedge.cond=Pedge.cond))
}
F_ResampleTreePLN <- function(Y, X, O, v=0.8, B=1e2){
   # v = 0.8; B = 1e2
   # Resampling edge probability
   # Y, X, O: same as for PLN
   # v = (subsample size) / (total sample size)
   # B = nb resamples
   # Out = Pmat = B x p(p-1)/2 matrix with edge probability for each resample
   # Y = Data$count; X = matrix(1, n, 1); O = matrix(0, n, p)
   n = nrow(Y); p = ncol(Y); P = p*(p-1)/2; V = round(v*n); Pmat = matrix(0, B, P); 
   for (b in 1:B){
      cat('\n', b, '')
      sample = sample(1:n, V, replace=F)
      Y.sample = Y[sample, ]; X.sample = X[sample, ]; O.sample = O[sample, ]
      # if(ncol(Y)>1){Y.sample = Y[sample, ]}else{Y.sample = matrix(Y[sample], V, 1)}; 
      # if(ncol(X)>1){X.sample = X[sample, ]}else{X.sample = matrix(X[sample], V, 1)}; 
      # if(ncol(O)>1){O.sample = O[sample, ]}else{O.sample = matrix(O[sample], V, 1)}; 
      PLN.sample = PLN(Y.sample ~ -1 + X.sample + offset(O.sample))
      Sigma.sample = PLN.sample$model_par$Sigma
      Pmat[b, ] = F_Sym2Vec(TreeGGM(cov2cor(Sigma.sample))$Pedge.cond)
   }
   return(Pmat)
}

# # Fit VEM-PLN with 1 = no covariates, 2 = depth+temp, 3 = all covariates
# VEM = list()
# VEM[[1]] = PLN(Data$count ~ 1)
# VEM[[2]] = PLN(Data$count ~ Data$covariates[, 1:2])
# VEM[[3]] = PLN(Data$count ~ Data$covariates)
# save(VEM, file=paste0(data.dir, data.name, '-VEM.Rdata'))
load(paste0(data.dir, data.name, '-VEM.Rdata'))
M = length(VEM)

# # Fit TreeGGM
# EMtree = list()
# invisible(sapply(1:M, function(m){
#    EMtree[[m]] <<- TreeGGM(cov2cor(VEM[[m]]$model_par$Sigma))
#    plot(F_Sym2Vec(EMtree[[m]]$Pedge.marg), F_Sym2Vec(EMtree[[m]]$Pedge.cond),
#         xlim=c(0, 1), ylim=c(0, 1))
# }))
# save(EMtree, file=paste0(data.dir, data.name, '-EMtree.Rdata'))
load(paste0(data.dir, data.name, '-EMtree.Rdata'))
invisible(sapply(1:M, function(m){cat(VEM[[m]]$loglik, EMtree[[m]]$log.pY[length(EMtree[[m]]$log.pY)], '\n')}))
# Pedge = cbind(F_Sym2Vec(EMtree[[1]]$Pedge.marg), F_Sym2Vec(EMtree[[2]]$Pedge.marg), 
              # F_Sym2Vec(EMtree[[3]]$Pedge.marg))
Pedge = cbind(F_Sym2Vec(EMtree[[1]]$Pedge.cond), F_Sym2Vec(EMtree[[2]]$Pedge.cond), 
              F_Sym2Vec(EMtree[[3]]$Pedge.cond))
colnames(Pedge) = c('M.null', 'M.env', 'M.all')
# plot(as.data.frame(Pedge))
# plot(as.data.frame(qlogis(Pedge)))
par(mfrow=c(2, 2))
for (m1 in (1:(M-1))){for (m2 in ((m1+1):M)){plot(qlogis(Pedge[, m1]),qlogis(Pedge[, m2]))}}
invisible(sapply(1:M, function(m){
   if(m==1){plot(sort(Pedge[, m]), type='b', ylim=c(0, 1)); abline(h=2/p)}else{points(sort(Pedge[, m]), col=m, type='b')}
   }))

# Resampling
Stab.sel = list()
X = list(); X[[1]] = matrix(1, n, 1); X[[2]] =  cbind(rep(1, n), Data$covariates[, 1:2]);
X[[3]] = cbind(rep(1, n), Data$covariates)
invisible(sapply(1:M, function(m){
   Stab.sel[[m]] <<- F_ResampleTreePLN(Data$count, X[[m]], matrix(0, n, p), B=5e2)
   }))
save(Stab.sel, file=paste0(data.dir, data.name, '-StabSel.Rdata'))
load(paste0(data.dir, data.name, '-StabSel.Rdata'))

# Edge selection and comparisons
par(mfrow=c(2, 2))
edge.sel = freq.sel = list()
invisible(sapply(1:M, function(m){
   freq.sel[[m]] <<- colMeans(1*(Stab.sel[[m]] > 2/p))
   hist(freq.sel[[m]], breaks=p)
   edge.sel[[m]] <<- 1*(freq.sel[[m]] > freq.sel.thres)
   }))
invisible(sapply(1:M, function(m){
   if(m==1){plot(sort(freq.sel[[m]]), type='b'); abline(h=freq.sel.thres)}else{points(sort(freq.sel[[m]]), col=m, type='b')}
   }))

par(mfrow=c(2, 2))
invisible(sapply(1:(M-1), function(m1){sapply((m1+1):M, function(m2){
   print(table(edge.sel[[m1]], edge.sel[[m2]]))
   edge.col = 1 + (2*edge.sel[[m1]] + edge.sel[[m2]])
   plot(qlogis(Pedge[, m1]), qlogis(Pedge[, m2]), col=edge.col)
   abline(h=qlogis(2/p), v=qlogis(2/p))
})}))

# Networks
node.coord = gplot(F_Vec2Sym(edge.sel[[M]]), gmode='graph') 
par(mfrow=c(2, 2))
invisible(sapply(1:M, function(m){
  gplot(F_Vec2Sym(edge.sel[[m]]), gmode='graph', coord=node.coord, label=species, main=sum(edge.sel[[m]])) 
}))

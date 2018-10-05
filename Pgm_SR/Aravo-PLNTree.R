# PLNtree for Barents fish data

rm(list=ls()); par(pch=20)
library(PLNmodels); library(sna); library(ade4)
source('../pack1/R/codes/FunctionsMatVec.R')
source('../pack1/R/codes/FunctionsTree.R')
source('../pack1/R/codes/FunctionsInference.R')

# Data
data.dir = '../Data_SR/'
data(aravo); data.name = 'Aravo'
Data = list(count=as.matrix(aravo$spe), covariates=aravo$env)
Data$covariates$Aspect = as.factor(Data$covariates$Aspect)
Data$covariates$Form = as.factor(Data$covariates$Form)
Data$covariates$ZoogD = as.factor(Data$covariates$ZoogD)
colnames(Data$count) = aravo$spe.names
n = nrow(Data$count); p = ncol(Data$count); species = colnames(Data$count)

# Algo parms
freq.sel.thres = 0.80; iter.max = 1e2; B.resample = 5e2
VEM.fit = F; Tree.fit = T; Tree.res = T

# Fit VEM-PLN with 1 = no covariates, 2 = all covariates
if(VEM.fit){
   VEM = list()
   VEM[[1]] = PLN(count ~ 1, data=Data)
   VEM[[2]] = PLN(count ~ covariates$Aspect + covariates$Slope + covariates$Form + covariates$PhysD + covariates$ZoogD + covariates$Snow, data=Data)
   save(VEM, file=paste0(data.dir, data.name, '-VEM.Rdata'))
}
load(paste0(data.dir, data.name, '-VEM.Rdata'))
M = length(VEM)

# Fit TreeGGM
if(Tree.fit){
   EMtree = list()
   invisible(sapply(1:M, function(m){
      EMtree[[m]] <<- TreeGGM(cov2cor(VEM[[m]]$model_par$Sigma), n, step='FALSE', maxIter=300)
      plot(F_Sym2Vec(EMtree[[m]]$P), F_Sym2Vec(EMtree[[m]]$probaCond),
           xlim=c(0, 1), ylim=c(0, 1))
   }))
   save(EMtree, file=paste0(data.dir, data.name, '-EMtree.Rdata'))
}
load(paste0(data.dir, data.name, '-EMtree.Rdata'))
par(mfrow=c(M, 1)); sapply(1:M, function(m){plot(EMtree[[m]]$L)})

# Compare edge probabilities
invisible(sapply(1:M, function(m){cat(VEM[[m]]$loglik, EMtree[[m]]$log.pY[length(EMtree[[m]]$log.pY)], '\n')}))
Pedge = cbind(F_Sym2Vec(EMtree[[1]]$P), F_Sym2Vec(EMtree[[2]]$P))
colnames(Pedge) = c('M.null', 'M.all')
par(mfrow=c(2, 2))
for (m1 in (1:(M-1))){for (m2 in ((m1+1):M)){plot(qlogis(Pedge[, m1]),qlogis(Pedge[, m2]))}}
invisible(sapply(1:M, function(m){
   if(m==1){plot(sort(Pedge[, m]), type='b', ylim=c(0, 1)); abline(h=2/p)}else{points(sort(Pedge[, m]), col=m, type='b')}
}))

# Resampling
if(Tree.res){
   Stab.sel = list()
   X = list(); X[[1]] = matrix(1, n, 1); X[[2]] =  cbind(rep(1, n), Data$covariates[, 1:2]);
   X[[3]] = cbind(rep(1, n), Data$covariates)
   invisible(sapply(1:M, function(m){
      Stab.sel[[m]] <<- F_ResampleTreePLN(Data$count, X[[m]], matrix(0, n, p), B=B.resample, maxIter=300, cond.tol=1e-8)
   }))
   save(Stab.sel, file=paste0(data.dir, data.name, '-StabSel.Rdata'))
}
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


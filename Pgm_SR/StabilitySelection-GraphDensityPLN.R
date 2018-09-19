# Network density estimation via stability selection

rm(list=ls()); par(pch=20, mfrow=c(1, 1))
source('Functions/FunctionsSimul.R')
source('../pack1/R/codes/FunctionsMatVec.R')
source('../pack1/R/codes/FunctionsTree.R')
source('../pack1/R/codes/FunctionsInference.R')
# source('../pack1/R/codes/TreeMixture-RML.R')
library(mvtnorm); library(PLNmodels); library(sna)

# Functions
TreeGGM <- function(CorY){
  p = ncol(CorY);
  phi = 1/sqrt(1 - CorY^2); diag(phi) = 0
  beta.unif = matrix(1, p, p); 
  FitEM = FitBetaStatic(Y, beta.init=beta.unif, phi=phi)
  return(list(P=Kirshner(FitEM$beta)$P, L=FitEM$logpY))
}
F_ResampleTreePLN <- function(Y, X, O, v=0.8, B=1e2){
  # v = 0.8; B = 1e2
  # Resampling edge probability
  # Y, X, O: same as for PLN
  # v = (subsample size) / (total sample size)
  # B = nb resamples
  # Out = Pmat = B x p(p-1)/2 matrix with edge probability for each resample
  V = round(v*n); Pmat = matrix(0, B, P); 
  for (b in 1:B){
    cat('\n', b, '')
    sample = sample(1:n, V, replace=F)
    Y.sample = Y[sample, ]; X.sample = X[sample, ]; O.sample = O[sample, ];
    # Use PLN 'inception' option to accelerate: problem with the dimension test (cannot update PLN$n)
    # PLN.init = PLN;
    # PLN.init$var_par = list(M=PLN$var_par$M[sample, ], S=PLN$var_par$S[sample, ])
    # control$inception = PLN.init
    # PLN.sample = PLN(Y.sample ~ -1 + X.sample + offset(O.sample), control=control)
    PLN.sample = PLN(Y.sample ~ -1 + X.sample + offset(O.sample))
    Sigma.sample = PLN.sample$model_par$Sigma
    Pmat[b, ] = F_Sym2Vec(TreeGGM(cov2cor(Sigma.sample), "FALSE", FALSE)$P)
    # if(b%%10==0){boxplot(log10(Pmat[1:b, ]) ~ Gmat[1:b, ])}
  }
  return(Pmat)
}

# Sim parms
p = 15; n = 50; c = 3; d = 2.5/p; P = p*(p-1)/2
seed = 1; set.seed(seed)

# Model parms
G = SimErdos(p, d); Gvec = F_Sym2Vec(G)
gplot(G, gmode='graph')
D = rep(1, p); lambda = 1; Omega = lambda*D + G; 
while (min(eigen(Omega)$values) < 1e-2){lambda = 1.1*lambda; Omega = lambda*diag(D) + G}
Sigma = solve(Omega); O = matrix(rnorm(n*p), n, p)

# Data simulation
Z = rmvnorm(n, sigma=Sigma); 
X = matrix(rnorm(n*c), n, c)
beta = matrix(rnorm(c*p), c, p)
Lambda = exp(O+X%*%beta+Z)
Y = matrix(rpois(n*p, Lambda), n, p); 
# plot(Lambda, Y); abline(0, 1)

# Algo parms
tol.def = list(trace=0, ftol_rel=1e-6, ftol_abs=1e-6, xtol_rel=1e-4, xtol_abs=1e-4)
control = tol.def

# PLN inference
PLN = PLN(Y ~ -1 + X + offset(O)); 
Sigma.hat = PLN$model_par$Sigma
Pvec = F_Sym2Vec(TreeGGM(cov2cor(Sigma.hat), "FALSE", FALSE)$P)
Porder = order(Pvec); plot(Pvec[Porder], col=1+Gvec[Porder]); abline(h=2/p)

# Resampled edge probabilities
Pmat1 = F_ResampleTreePLN(Y, X, O, B=5e2)
Pmat2 = F_ResampleTreePLN(Y, X, O, v=0.5, B=5e2)
Pmat = Pmat1

# Stability selection
Pfreq = 1*(Pmat > 2/p) # Threshold at tree density
Psel = (colMeans(Pfreq) > .5) # Keep edges selected more than half of the time
table(Gvec, Psel)
cat(sum(Gvec), sum(Psel))
Porder = order(Pvec); plot(Pvec[Porder], col=1+Gvec[Porder]); 
abline(h=2/p); abline(h=quantile(Pvec, prob=1-sum(Psel)/P), col=2)



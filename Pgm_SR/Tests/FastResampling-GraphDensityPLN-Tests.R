# Fast resampling for PL-network density estimaition
# Based on Schur's trick used in https://arxiv.org/pdf/1505.07281.pdf, p23

rm(list=ls())
source('Functions/FunctionsSimul.R')
source('Functions/FunctionGGMpvalResampling.R')
library(mvtnorm); library(PLNmodels); library(sna)

# Sim parms
p = 30; n = 100; d = .1; P = p*(p-1)
par(pch=20); B = 1e3
seed = 1; set.seed(seed)

# Model parms
G = SimErdos(p, d); 
D = rep(1, p); lambda = 1; Omega = lambda*D + G; 
while (min(eigen(Omega)$values) < 1e-2){lambda = 2*lambda; Omega = lambda*diag(D) + G}
Sigma = solve(Omega)
O = matrix(rnorm(n*p), n, p)

# Data
Z = rmvnorm(n, sigma=Sigma); 
Y = matrix(rpois(n*p, exp(O+Z)), n, p)

# Test
par(mfrow=c(2, 3))
k = 3; LM = lm(Z[, k] ~ -1 + Z[, -k], x=T);
ZpZ = t(Z)%*%Z; invZpZ = solve(ZpZ)
invXpX = summary(LM)$cov.unscaled;
invxpx = invZpZ[-k, -k] - outer(invZpZ[-k, k], invZpZ[k, -k]) / invZpZ[k, k]; plot(invxpx, invXpX); abline(0, 1);
xpy = ZpZ[-k, k]; invypy = invZpZ[k, k]
Beta.k = summary(LM)$coef[, 1]; beta.k = invxpx%*%xpy; plot(beta.k, Beta.k); abline(0, 1);
Sd.k = summary(LM)$coef[, 2]; sd.k = sqrt(diag(invxpx)) / sqrt(invypy) / sqrt(n-p+1); plot(sd.k, Sd.k); abline(0, 1);
Stat.k = summary(LM)$coef[, 3]; stat.k = beta.k/sd.k; plot(stat.k, Stat.k); abline(0, 1)

# Test statistics using LM on Z
Stat.lm = matrix(NA, p, p)
invisible(sapply(1:p, function(k){
   LM = lm(Z[, k] ~ -1 + Z[, -k])
   Stat.lm[-k, k] <<- summary(LM)$coef[, 3]
}))

# Inference based on PLN
PLN = PLN(Y ~ -1 + offset(O))
M = PLN$var_par$M; S = PLN$var_par$S

# Theoretical p-values (based on Z)
Stat.th = F_StatAll(t(Z)%*%Z, n, p)
Pval.th = 2*pt(abs(Stat.th), df=(n-p-1), lower.tail=F)

# Test statistics based on PLN inference
ZpZ = t(M)%*%M + diag(colSums(S)); 
Stat = F_StatAll(ZpZ, n, p)

# Inferred covariance matrix
par(mfrow=c(2, 2))
gplot(G, gmode='graph');
plot(Sigma, cov(Z), col=(1+G)); abline(0, 1)
plot(Sigma, ZpZ/n, col=(1+G)); abline(0, 1)
plot(cov(Z), ZpZ/n, col=(1+G)); abline(0, 1)
gplot(G, gmode='graph');
plot(Omega, solve(cov(Z)), col=(1+G)); abline(0, 1)
plot(Omega, solve(ZpZ/n), col=(1+G)); abline(0, 1)
plot(solve(cov(Z), solve(ZpZ/n)), col=(1+G)); abline(0, 1)

# Bootstrap p-values estimation
Pval.bt = F_ResamplePval(M, S, Stat, B)

# Fast bootstrap p-values estimation
Pval.fast = F_FastResamplePval(M, S, Stat, B)

# Comparison of p-values
par(mfrow=c(3, 3))
hist(Pval.th, breaks=p); hist(Pval.bt, breaks=p); hist(Pval.fast, breaks=p)
plot(0, col=0);
plot(Pval.th, Pval.bt, col=(1+G)); abline(0, 1)
plot(Pval.th, Pval.fast, col=(1+G)); abline(0, 1)
plot(0, col=0);
plot(qnorm(Pval.th), qnorm(Pval.bt), col=(1+G)); abline(0, 1)
plot(qnorm(Pval.th), qnorm(Pval.fast), col=(1+G)); abline(0, 1)

P0.th = 2*sum(Pval.th >= .5, na.rm=T); P1.th = P - P0.th
P0.bt = 2*sum(Pval.bt >= .5, na.rm=T); P1.bt = P - P0.bt
P0.fast = 2*sum(Pval.fast >= .5, na.rm=T); P1.fast = P - P0.fast
cat(P0.th, P1.th, '/', P0.bt, P1.bt, '/', P0.fast, P1.fast, '/', sum(G), '\n')

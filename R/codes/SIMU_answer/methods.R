library(mvtnorm) 
library(clusterGeneration) 
library(mgcv) 
library(statmod) 
library(corpcor)
library(ecoCopula)

########################
# Counts GCGM
########################
simulate_gcgm_p<-function(N,P,sig,env,mmean,effect){
  #simulate latent Gaussian random variables with covariance matrix sig 
  Z=rmvnorm(n=N,mean = rep(0,P),sigma = sig)
  #transform to uniform (copula) scale
  U=pnorm(Z)
  #add environmental heterogeneity and transform to Poisson scale 
  meanmat=cbind(1,env)%*%rbind(mmean,effect) 
  Y=qpois(U,lambda=exp(meanmat))
  Y
}
# We use the following values for the parameters
n_spp = 9
n_sites =50
env = matrix(rep(c(0,1),n_sites/2), ncol = 1)
env_unif= 0 #for uniform env 
env_het=5 # for heterogeneous env
effect=rnorm(n_spp,0,env_unif)
mmean=rnorm(n_spp, 0) #random Gaussian
# the covariance matrix is generated using library(clusterGeneration)
# The "onion" method was used as it gave precision matrices
# which were most distinct from the correlation matrix
sig=genPositiveDefMat(n_spp,covMethod="onion")$Sigma
# The "true" interactions are the negative values of the precision matrix, i.e.
my_gr=-solve(sig)

my_samp=simulate_gcgm_p(N = n_sites,P=n_spp,sig=sig,env=env,mmean=mmean,effect=effect)

########################
# Count Markov adapt. harris 2016
########################
make_coefficients = function(n_spp, p_neg, mean_alpha){
  # Exponential distribution has lots of mass near 0 but has
  # a long tail.
  true_beta_magnitudes = rexp(choose(n_spp, 2), rate = 1)
  
  # Multiply some proportion of the interactions
  # by -1
  b = true_beta_magnitudes * sample(
    c(-1, 1), 
    size = length(true_beta_magnitudes), 
    prob = c(p_neg, 1 - p_neg),
    replace = TRUE
  )
  
  # Species' intercepts are normally distributed
  a = rnorm(n_spp, mean_alpha)
  
  
  # Return the simulated values.
  # The rosalia function stores pairwise parameters in the upper
  # triangle of an n-by-n matrix and stores the species intercepts
  # along the diagonal, so these values are named accordingly.
  c(alpha  = a, beta = b)
}
simulate_markov_p = function(n_spp, n_sites, n_gibbs, f, rdist, p_neg){
  par = make_coefficients(n_spp, p_neg, 0) # mean_alpha=0
  truth = par[-(1:n_spp)]
  alpha = par[1:n_spp]
  # Turn the interaction values into an n-by-n matrix
  # Start with empty matrix; fill in upper triangle;
  # then fill in lower triangle with its transpose
  beta = matrix(0, n_spp, n_spp)
  beta[upper.tri(beta)] = truth
  beta = beta + t(beta)
  
  #environment is binary or normally distributed
  #env = matrix(rep(c(0,1),n_sites/2), ncol = 1)
  env = matrix(rnorm(n_sites * n_env), ncol = n_env)
  alpha_env = matrix(rnorm(n_spp * n_env, sd = sd), nrow = n_env)
  
  # Simulate the landscape from known process with Gibbs sampling
  # Landscape starts as if betas were all zero. Each species' occurrence 
  # probability or abundance depends on its alpha value and on the
  # environment (assuming alpha_env is not set to zero).
  
  x = matrix(
    f(rep(1, n_sites) %*% t(alpha) + env %*% alpha_env), nrow = n_sites,
    ncol = n_spp
  )
  # Gibbs sampling
  for(i in 1:n_gibbs){
    # Each round of Gibbs sampling updates one species (column) across all sites 
    #according to its conditional probability (i.e. conditional on environment 
    # and the other species that are present).
    for(j in 1:n_spp){
      x[,j] = rdist( nrow(x),
                     f(x %*% beta[ , j] + alpha[j] + env %*% alpha_env[,j]))
    } 
  }
  out=list(x=x,beta=beta)
  out 
}
#We use the following values for the parameters
n_spp = 9
n_sites =50
n_gibbs = 5000
f = exp# (log link)
p_neg=1 # for abundance data in harris, 0.75 for binary data 
simu=simulate_markov_p(n_spp, n_sites, n_gibbs, f, rdist=rpois, p_neg)
########################
# Estimation models
########################
data_fam = "Poisson"
#1. gllvm
library(gllvm)
# the data is in a matrix my_samp
# when the enviroment is heterogeneous we specify the following model 
env = matrix(rep(c(0,1),N/2), ncol = 1) #environmental covariate
het_gllvm=gllvm(my_samp, NULL, formula = ~ env, num.lv = 5, family = "negative.binomial") 
# when the enviroment is uniform we specify the following 
unif_gllvm=gllvm(my_samp, NULL, formula = ~ 1 ,num.lv = 5, family = data_fam) # data_fam = "binomial" for binary
# data_fam = "Poisson" for counts
# the covariance matrix obtained with the getResidualCor function is compared to direct species associations
resgllvm=getResidualCor(het_gllvm)


#2 ggmlog
library(glasso)
# the data is in a matrix my_samp 
unif_ggmlog=glasso(cor(log(my_samp+1)),rho = 0))

#3. MRF cov
library(MRFcov)
# the data is in a matrix my_samp
# when the enviroment is heterogeneous we specify the following model
env = matrix(rep(c(0,1),n_sites/2), ncol = 1) #environmental covariate 
clark_dat=data.frame(my_samp) 
het_gllvm=try(MRFcov(cbind(clark_dat,env),n_nodes = n_spp,family = "poisson")) 
het_gllvm$graph
# when the enviroment is uniform we specify the following model 
unif_gllvm=MRFcov(clark_dat,n_nodes = P,family = data_fam)
# data_fam = "binomial" for binary
# data_fam = "Poisson" for counts

# 4. ecocopula

#ecoCopula works by first estimating a model, and then using the fitted model to estimate direct associations.
#For binary and count data we use the manyglm function from the mvabund package to estimate the model.
library(ecoCopula)
# the data is in a matrix my_samp
# when the enviroment is heterogeneous we specify the following model 
env = matrix(rep(c(0,1),N/2), ncol = 1) 
#environmental covariate 
my_mod=manyglm(my_samp~env, family="data_fam")
# when the enviroment is uniform we specify the following model 
my_mod=manyglm(my_samp~env, family="data_fam")#
# data_fam = "binomial" for binary
# data_fam = "negativ.binomial" for counts (default in manyglm)


#The direct species associations are then estimated with the cgr function from the ecoCopula package.
cgr(my_mod,lambda=0)



##############
library(mvabund)
library(ecoCopula)
data(spider)
abund <- mvabund(spider$abund)
X <- spider$x
spider_mod=manyglm(abund~X)
spid_graph=cgr(spider_mod)
plot(spid_graph,pad=1)
spid_graph$best_graph

# score ecoCopula
str(spid_graph$all_graphs)

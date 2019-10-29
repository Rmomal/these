################################################################################
# File: gcoda.R
# Aim : Conditional dependence network inference for compositional data
#-------------------------------------------------------------------------------
# Author: Fang Huaying (Peking University)
# Email : hyfang@pku.edu.cn
# Date  : 13MAY2016
#-------------------------------------------------------------------------------
require(huge);
#-------------------------------------------------------------------------------
# Function parameter description:
# function: gcoda
#   Input:
#                 x ------ n x p data matrix (row/column is sample/variable)
#                          n samples & p compositional variables
#            counts ------ Is the compositional data matrix a count matrix?
#                          Default: FALSE
#            pseudo ------ pseudo count if counts = TRUE
#                          Default: 0.5
#  lambda.min.ratio ------ lambda(max) / lambda(min)
#                          Default: 1e-4
#           nlambda ------ number of tuning parameters
#                          Default: 15
#        ebic.gamma ------ gamma value of EBIC
#                          Default: 0.5
#   Output:
#      A list structure contains:
#            lambda ------ lambda sequence for compuation of EBIC score
#           nloglik ------ negative log likelihood for lambda sequence
#                df ------ number of edges for lambda sequence
#              path ------ sparse pattern for lambda sequence
#              icov ------ inverse covariance matrix for lambda sequence
#        ebic.score ------ EBIC score for lambda sequence
#             refit ------ sparse pattern with best EBIC score
#          opt.icov ------ inverse covariance matrix with best EBIC score
#        opt.lambda ------ lambda with best EBIC score
#-------------------------------------------------------------------------------
gcoda <- function(x, counts = F, pseudo = 0.5, lambda.min.ratio = 1e-1,
                  nlambda = 30, ebic.gamma = 0.5, covar=NULL,shuffleRes=FALSE, shufflevar) {

  # Counts or fractions?
  if(counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
  }
  n <- nrow(x);
  p <- ncol(x);
  # Log transformation for compositional data
  # covariates taken into account ?
  if(is.null(covar)){
    S <- var(log(x) - rowMeans(log(x)))
  }else{

    x<-log(x) - rowMeans(log(x))
    # m<- model.matrix(~X1+X2+X3,covar)
    #m<- model.matrix(~X1,covar)
    string<-paste("x", paste(covar, collapse=" + "), sep=" ~ ")
    formula<-as.formula(string)
    X = as.matrix(lm(formula, x=T)$x)
    model<-lm(x~X)
    U<-model$residuals
   
    if(shuffleRes){ #browser()
    U<- shuffle(U,shufflevar)}
    S <- var(U)
  }

  # Generate lambda via lambda.min.ratio and nlambda
  lambda.max <- max(max(S - diag(p)), -min(S - diag(p)));
  lambda.min <- lambda.min.ratio * lambda.max;
  lambda <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda));
  # Store fit result for gcoda via a series of lambda
  fit <- list();
  fit$lambda <- lambda;
  fit$nloglik <- rep(0, nlambda);
  fit$df <- rep(0, nlambda);
  fit$path <- list();
  fit$icov <- list();
  # Compute solution paths for gcoda
  icov <- diag(p);
  for(i in 1:nlambda) {
    out.gcoda <- gcoda_sub(A = S, iSig = icov, lambda = lambda[i]);
    icov <- out.gcoda$iSig;
    fit$nloglik[i] <- out.gcoda$nloglik;
    fit$icov[[i]] <- icov;
    fit$path[[i]] <- 0 + (abs(icov) > 1e-6);
    diag(fit$path[[i]]) <- 0;
    fit$df[i] <- sum(fit$path[[i]]) / 2;
  }
  # Continue if edge density is too small
  imax <- nlambda + 15;
  emin <- p * (p - 1) / 2 * 0.618;
  while(fit$df[i] <= emin && i <= imax && lambda[i] > 1e-6) {
    cat(i, "Going on!\n");
    lambda <- c(lambda, lambda[i]/2);
    i <- i + 1;
    fit$lambda[i] <- lambda[i];
    out.gcoda <- gcoda_sub(A = S, iSig = icov, lambda = lambda[i]);
    icov <- out.gcoda$iSig;
    fit$nloglik[i] <- out.gcoda$nloglik;
    fit$icov[[i]] <- icov;
    fit$path[[i]] <- 0 + (abs(icov) > 1e-6);
    diag(fit$path[[i]]) <- 0;
    fit$df[i] <- sum(fit$path[[i]]) / 2;
  }
  # Compute EBIC score for lambda selection
  fit$ebic.score <- n * fit$nloglik + log(n) * fit$df +
    4 * ebic.gamma * log(p) * fit$df;
  fit$opt.index <- which.min(fit$ebic.score);
  fit$refit <- fit$path[[fit$opt.index]];
  fit$opt.icov <- fit$icov[[fit$opt.index]];
  fit$opt.lambda <- fit$lambda[fit$opt.index];
  return(fit);
}
#-------------------------------------------------------------------------------
# Optimization for gcoda with given lambda
gcoda_sub <- function(A, iSig = NULL, lambda = 0.1, tol_err = 1e-4,
                      k_max = 500) {
  p <- ncol(A);
  if(is.null(iSig)) {
    iSig <- diag(p);
  }

  err <- 1;
  k <- 0;
  fval_cur <- Inf;

  while(err > tol_err && k < k_max) {
    iSig_O <- rowSums(iSig);
    iS_iSig <- 1 / sum(iSig_O);
    iSig_O2 <- iSig_O * iS_iSig;
    A_iSig_O2 <- rowSums(A * rep(iSig_O2, each = p));
    A2 <- A - A_iSig_O2 - rep(A_iSig_O2, each = p) +
      sum(iSig_O2 * A_iSig_O2) + iS_iSig;
    iSig2 <- huge_glasso_mod(S = A2, lambda = lambda);

    fval_new <- obj_gcoda(iSig = iSig2, A = A, lambda = lambda);
    xerr <- max(abs(iSig2 - iSig) / (abs(iSig2) + 1));
    err <- min(xerr, abs(fval_cur - fval_new)/(abs(fval_new) + 1));

    if(is.na(err)){
      cat("BREAK")
      break
    }

    k <- k + 1;
    iSig <- iSig2;
    fval_cur <- fval_new;
  }
  nloglik <- fval_cur - lambda * sum(abs(iSig));

  if(k >= k_max) {
    cat("WARNING of gcoda_sub:\n", "\tMaximum Iteration:", k,
        "&& Relative error:", err, "!\n");
  }

  return(list(iSig = iSig, nloglik = nloglik));
}
#----------------------------------------
# Objective function value of gcoda (negative log likelihood + penalty)
obj_gcoda <- function(iSig, A, lambda) {
  p <- ncol(A);
  iSig_O <- rowSums(iSig);
  S_iSig <- sum(iSig_O);
  nloglik <- - log(det(iSig)) + sum(iSig * A) + log(S_iSig) -
    sum(iSig_O * rowSums(A * rep(iSig_O, each = p))) / S_iSig;
  pen <- lambda * sum(abs(iSig));
  return(nloglik + pen);
}
#----------------------------------------
# Modified huge::huge.glasso for quick preparation
# Input S must be covariance matrix
require(huge);
sourceCpp("/Users/raphaellemomal/these/R/codes/hugeglasso.cpp")
huge_glasso_mod <- function(S, lambda) {
  icov <- diag(1/(diag(S) + lambda));
  z <- which(rowSums(abs(S) > lambda) > 1);
  q <- length(z);
  if (q > 0) {
    out.glasso= hugeglasso_sub(S = as.matrix((S[z, z])), W = as.matrix((S[z, z])), T = as.matrix(diag(as.double(q))),
                               d= as.integer(q), ilambda = as.double(lambda), df = as.integer(0), 
                               scr=TRUE)
    icov[z, z] = matrix(out.glasso, ncol = q);
  }

  return(icov);
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 1 Basic example (no edges in the conditional dependence network)
# 1.1 Generate logistic normal variables
# n <- 100;
# p <- 20;
# x <- matrix(rnorm(n * p), nrow = n);
# x.frac <- exp(x) / rowSums(exp((x)));
# totCount <- round(runif(n = n,  min = 1000, max = 2000));
# x.count <- x.frac * totCount;
# # 1.2 Run gCoda
# # using fraction
# res_gcoda_frac <- gcoda(x = x.frac, counts = F);
# # using counts
# res_gcoda_count <- gcoda(x = x.count, counts = T);
# # 1.3 Get the estimation of the inverse covariance matrix
# {
#   cat("gCoda using fraction data:\n");
#   print(round(res_gcoda_frac$opt.icov, 2));
#   cat("gCoda using count data:\n");
#   print(round(res_gcoda_count$opt.icov, 2));
# }

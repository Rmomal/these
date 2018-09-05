# Fast resampling for PL-network density estimaition
# Based on Schur's trick used in https://arxiv.org/pdf/1505.07281.pdf, p23

# Functions
F_Statk <- function(xpy, invxpx, invypy, n, p){
   # Computes the test statistics of the regression coefficients for variable k as a response
   # Denoting Z = [X, Y]
   # xpy = X'Y, invxpx = (X'X)^(-1), invypy
   stat = rep(NA, p-1)
   beta = as.vector(invxpx%*%xpy)
   sd = as.vector(sqrt(diag(invxpx)) / sqrt(invypy) / sqrt(n-p+1))
   stat = beta/sd
   return(stat)
}
F_StatAll <- function(zpz, n, p){
   # Computes the test statistics of the regression coefficients for all variables 
   # ZpZ = (approximation of) Z' * Z; n, p = dimensions
   invzpz = solve(zpz)
   stat = matrix(NA, p, p)
   sapply(1:p, function(k){
      xpy = zpz[-k, k]; invypy = invzpz[k, k]
      invxpx = invzpz[-k, -k] - outer(invzpz[-k, k], invzpz[k, -k]) / invzpz[k, k]
      stat[-k, k] <<- F_Statk(xpy, invxpx, invypy, n, p)
   })
   return(stat)
}
F_ResamplePval <- function(m, s, stat, B=1e3){
   # Computes the p-values for all test statistics by resampling
   # m = matrix of approximate mean of Z; s = matrix of approximate variance
   # Stat = matrix of test statistics; B = number of resampling
   n = nrow(m); p = ncol(m)
   pval = matrix(0, p, p); diag(pval) = NA
   cat('k =')
   for (k in 1:p){
      cat(k, '')
      for (b in 1:B){
         # Shuffling column k
         mbk = m; mbk[, k] = mbk[order(runif(n)), k]
         # Pseudo covariance matrix
         zpz = t(mbk)%*%mbk + diag(colSums(s)); 
         xpy = zpz[-k, k]
         invzpz = solve(zpz)
         invxpx = invzpz[-k, -k] - outer(invzpz[-k, k], invzpz[k, -k]) / invzpz[k, k]
         invypy = invzpz[k, k]
         statbk = F_Statk(xpy, invxpx, invypy, n, p)
         pval[-k, k] = pval[-k, k] + 1*(abs(statbk) > abs(stat[-k, k]))
      }
   }
   cat('\n')
   return(pval / B)
}
F_FastResamplePval <- function(m, s, stat, B=1e3){
   # Computes the p-values for all test statistics by resampling in a faster manner
   # m = matrix of approximate mean of Z; s = matrix of approximate variance
   # Stat = matrix of test statistics; B = number of resampling
   n = nrow(m); p = ncol(m)
   pval = matrix(0, p, p); diag(pval) = NA
   zpz = t(m)%*%m + diag(colSums(s)); invzpz = solve(zpz)
   cat('k =')
   for (k in 1:p){
      cat(k, '')
      for (b in 1:B){
         # Shuffling column k
         m.bk = m; m.bk[, k] = m.bk[order(runif(n)), k]
         # Pseudo covariance matrix
         zpz.bk = t(m.bk)%*%m.bk + diag(colSums(s)); 
         xpy = zpz.bk[-k, k]
         invxpx = invzpz[-k, -k] - outer(invzpz[-k, k], invzpz[k, -k]) / invzpz[k, k]
         invypy = 1/(zpz.bk[k, k] - zpz.bk[k, -k]%*%invxpx%*%zpz.bk[-k, k])
         statbk = F_Statk(xpy, invxpx, invypy, n, p)
         pval[-k, k] = pval[-k, k] + 1*(abs(statbk) > abs(stat[-k, k]))
      }
   }
   cat('\n')
   return(pval / B)
}


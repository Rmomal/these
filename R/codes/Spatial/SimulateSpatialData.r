##--------------------------------------------------------------------------------------------------------
## SCRIPT : Simulation de donn?es pour ?valuer l'impact de pseudo-absence
## 
## Authors : Matthieu Authier & Clara P?ron
## Last update : 2018-04-17
## R version 3.4.3 (2017-11-30) -- "Kite-Eating Tree"
##--------------------------------------------------------------------------------------------------------

### charger les biblioth?ques de fonctions n?cessaires
lapply(c("sp", "mvtnorm", "tidyverse", "fields", "raster"), library, character.only = TRUE)

rm(list = ls())

set.seed(20190426)
sqrt_n <- 10
n <- sqrt_n^2
grid <- data.frame(lon = rep(seq(0, 1, length.out = sqrt_n), times = sqrt_n), 
                   lat = rep(seq(0, 1, length.out = sqrt_n), each = sqrt_n)
                   )

### covariance functions
cov_exp <- function(distance, range = 1, sill = 1) {
  sill * sill * exp(- distance / range)
}
cov_matern <- function(distance, range = 1, sill = 1) {
  sill * sill * (1 + distance * sqrt(3) / range) * exp(-distance * sqrt(3) / range)
}

### effective range
n_pts <- 1e3
upper <- 6
data.frame(distance = rep(seq(0, upper, length.out = n_pts), 2),
           autocorrelation = c(cov_exp(distance = seq(0, upper, length.out = n_pts)),
                               cov_matern(distance = seq(0, upper, length.out = n_pts))
                               ),
           covariance = rep(c("Exponential", "Matern"), each = n_pts)
           ) %>% 
  ggplot(aes(x = distance, y = autocorrelation)) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  geom_vline(xintercept = c(exp(1), 3), linetype = "dotted", color = "red") +
  scale_x_continuous(breaks = c(0:upper)) +
  facet_wrap(~covariance, ncol = 1) +
  theme_bw()

## effective range is the distance for which autocorrelation is approx 0.05
## for exponential covariance it is approx 3 * range
## for matern (3/2) covariance it is approx exp(1) * range
# biggest difference between the two covariance fcts is behaviour at the origin

### covariance matrix
Omega <- as.matrix(cov_matern(distance = dist(grid, method = "euclidean")))
diag(Omega) <- 1
## cholesky decomposition
chol_omega <- t(chol(Omega))

# non-spatial covariate data
# X <- cbind(rep(1, n), replicate(4, rnorm(n)))

# spatial covariate data --> simulate values for 4 inputs
X <- cbind(rep(1, n), sapply(1:4, function(i) { as.numeric(chol_omega %*% rnorm(n)) })) # can also use drop(mvtnorm::rmvnorm(1, rep(0.0, n), Omega))
grid <- cbind(grid, X[, -1])
names(grid)[3:6] <- paste("x", 1:4, sep = "")

# here take the 4 inputs and use their values and their squared values as predictors
X <- cbind(X, X[, -1]^2)

sp_data <- function(df, varname) {
  spdf <- df[, c("lon", "lat", varname)]
  coordinates(spdf) = ~ lon + lat
  gridded(spdf) <- TRUE
  return(spdf)
}

plot(sp_data(df = grid, varname = "x1"))
plot(sp_data(df = grid, varname = "x2"))
plot(sp_data(df = grid, varname = "x3"))
plot(sp_data(df = grid, varname = "x4"))

# covariate effects
beta <- c(0, 0.5 * rnorm(4), 0.2 * rnorm(4))

grid$linpred <- as.numeric(X %*% beta)
plot(sp_data(df = grid, varname = "linpred")) # linear predictor

# response variable
grid$y_obs <- rpois(nrow(grid), exp(grid$linpred))
plot(sp_data(df = grid, varname = "y_obs")) # observed data

### add spatial noise
sill <- 1
grid$noise <- sill * as.numeric(chol_omega %*% rnorm(n))
plot(sp_data(df = grid, varname = "noise")) # spatial autocorrelation

grid$y_obs <- rpois(nrow(grid), exp(grid$linpred + grid$noise))
plot(sp_data(df = grid, varname = "y_obs")) # observed data

grid$latent <- exp(grid$linpred + grid$noise)
plot(sp_data(df = grid, varname = "latent")) # latent

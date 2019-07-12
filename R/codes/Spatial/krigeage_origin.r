
lapply(c("tidyverse", "lubridate", "rgdal", "geoR", "fields", "mvtnorm", "MASS", "reshape", "WhatIf"), library, character.only = TRUE)
theme_set(theme_minimal(base_size = 14))

##--------------------------------------------------------------------------------------------------------
## SCRIPT : Fonctions pour Krigeage sur données d'observations par transect
##
## Authors : Matthieu Authier & Fabien Barriau
## Last update : 2019-06-27
## R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
## Copyright (C) 2018 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

# variance inter-segments: possibilite de calculer la variance en tenant compte du design des transects (variance intra-transects)
emp_vario_intra <- function(y,
                            esw,
                            l,
                            center,
                            projected = TRUE,
                            id_transect = NULL,
                            breaks = NULL,
                            plot = TRUE,
                            pairs.min = 100,
                            n_sim = 100,
                            distmat = NULL,
                            region=NULL,
                            region.label=NULL,# region label to filter data accordingly
                            esp="species",  # only for title
                            perTransect = TRUE # should distances be calculated transect-wise ?
) {
  # y = données de comptage
  # esw = effective half-strip width en km
  # l = longueur en km des legs ou segments en km
  # center = matrice de dim(length(y), 2) avec les coordonnées XY des centres de chaque segment
  # id_transect = numero de code du transect (pour en tenir sinon laisser null)
  # breaks = break points for empirical variogram
  # distmat = distance matrix
  if(!is.null(region)){
    yorigin=y
    data=cbind(y,center, id_transect,l,region.label) %>% filter(region.label==region)
    y=data[,1]
    center=data[,2:3]
    id_transect=data[,4]
    l=data[,5]
    rm(data)
    
  }
  if(is.null(breaks)) {
    writeLines("Automatic, but potentially crappy, choice for distance bins of 2.5 km, capped at 100 km")
    breaks <- c(-0.01, 0.001, 1, seq(2.5, 100, 2.5))
  }
  
  ### fonction pour calculer le taux de rencontre moyen
  TauxObs <- function(y, esw, l){
    my_control <- list(epsilon = 1e-6, maxit = 100000, trace = FALSE)
    # utiliser quasi poisson pour eventuelle surdispersion
    
    model_m <- glm(y ~ 1 + offset(2 * esw * l), family = "quasipoisson", control = my_control)
    m <- mean(exp(rnorm(1000, as.numeric(coef(model_m)), as.numeric(sqrt(vcov(model_m))))))
    return(round(m, 8))
  }
  
  ### calculer la matrice de distance
  if(is.null(distmat)) {
    # tenir compte du design de l'échantillonnage
    if(is.null(id_transect)) { id_transect <- rep(1, length(y)) }
    else {
      offset_inter_transect <- 1000 # effort en km a rajouter à chaque changement de transect pour les rendre spatiallement indépendants
      if(!is.numeric(id_transect)) { id_transect <- as.numeric(as.factor(id_transect)) }
      ### trier dans l'ordre croissant
      croissant <- order(id_transect) ; id_transect <- id_transect[croissant]
      y <- y[croissant] ; l <- l[croissant] ; center <- center[croissant, ]
    }
    
    D <- ifelse(id_transect[-1] - id_transect[-length(y)] != 0, 1, 0)
    xD <- diffinv(D) # cumul de D
    if(perTransect){
      # classes de distances pour couple de points
      dl <- sqrt((center[-1, 1] - center[-length(y), 1])^2 + (center[-1, 2] - center[-length(y), 2])^2)
      xl <- diffinv(dl) + xD * offset_inter_transect
      convert2km <- ifelse(projected, 1/1e3, (60 * 1.852)) # convertir en km (si WGS84, 1°lat = 60'= 60 * 1.852 km soit environ 110km)
      distmat <- fields::rdist(xl) * convert2km
      distclas <- cut(distmat, breaks = breaks)
    }else{
      convert2km <- ifelse(projected, 1/1e3, (60 * 1.852)) # convertir en km (si WGS84, 1°lat = 60'= 60 * 1.852 km soit environ 110km)
      distmat <- fields::rdist(center,center)* convert2km
      breaks=c(-0.01, 0.01,1:4,seq(5,100,10), seq(105, min(300,max(distmat)),20))
      distclas <- cut(distmat, breaks = breaks)
    }
  }
  
 
  
  # calcul de la moyenne du taux d'observation m pour correction Poisson du Variogramme
 # m <- TauxObs(y, esw, l)
  m<-mean(y/(2*esw*l))
  # calcul du variogramme
  calc_vario <- function(y, l, esw, m) {
    return(0.5 *((outer(l, l, "*") * 2 * esw / outer(l, l, "+")) *(outer(y /(l * 2 * esw), y /(l * 2 * esw), "-"))^2 - m))
  }
  z2 <- calc_vario(y = y, l = l, esw = esw, m = m)
  
  # constante de pondération 
  cstpond <- tapply(outer(l, l, "*") * 2 * esw / outer(l, l, "+"), 
                    distclas, 
                    sum, na.rm = TRUE
  )
  
  #on divise par la constante de pondération
  vario1 <- tapply(z2, distclas, sum, na.rm = TRUE) / cstpond
  distvario1 <- tapply(distmat, distclas, mean, na.rm = TRUE)
  nbcouples1 <- tapply(z2, distclas, function(x) { length(which(!is.na(x))) })
  
  # enlever les estimations trop bruitées, ie pas assez de couples de points
  noise <- c(which(is.na(nbcouples1)), which(nbcouples1 < pairs.min))
  noise <- as.numeric(noise[order(noise)])
  vario1 <- ifelse(vario1 < 0, 0, vario1)
  
  if(length(noise) != 0) {
    vario1 <- vario1[-noise]
    distvario1 <- distvario1[-noise]
    nbcouples1 <- nbcouples1[-noise]
  }
  
  # mettre les résultats dans un objet vario de classe "variogram", package geoR
  vario <- list(u = as.numeric(distvario1), # distance
                v = as.numeric(vario1),     # estimated variogram values
                n = as.numeric(nbcouples1), # nombre de couples de points
                sd = rep(0, length(vario1)),
                bins.lim = breaks[-c(1, noise)],
                ind.bin = ifelse(as.numeric(nbcouples1) > pairs.min, TRUE, FALSE),
                var.mark = var(y /(2 * esw * l)),
                beta.ols = m,
                output.type = "bin",
                max.dist = max(breaks),
                estimator.type = "classical", # it should rather be Monestiezal
                n.data = length(y),
                lambda = 1,
                trend = "cte",
                pairs.min = pairs.min,
                nugget.tolerance = 1e-12,
                direction = "omnidirectional",
                tolerance = "none",
                uvec = as.numeric(distvario1),
                call = NULL
  )
  class(vario) <- "variogram"
  
  # test graphique
  if(plot) {
    # calculer une enveloppe nulle
    # Quelques fonctions simples pour faire du bootstrap et déterminer s'il faut utiliser une loi négative binomiale
    overdispersion <- function(x,m) { 
      # esp2=TauxObs(y^2,esw,l)
      # varx = esp2-m^2
      return(var(x)/mean(x))}
    permute_index <- function(x) { return(x[sample(c(1:length(x)), length(x), replace = TRUE)])}
    bootstrap_pval <- function(x, simulator, statistic, n_sim, alpha, truth) {
      tboots <- replicate(n_sim, statistic(simulator (x),m)) 
      ci_lower <- 2*statistic(x,m) - quantile(tboots, 1 - alpha/2, na.rm=TRUE)
      ci_upper <- 2*statistic(x,m) - quantile(tboots, alpha/2, na.rm=TRUE)
      return(as.numeric(ifelse(ci_lower > truth, 0, ifelse(ci_upper < truth, 0, 1))))
    }
    
    if(overdispersion(y,m) < 1) { negbin <- FALSE }
    else { negbin <- ifelse(bootstrap_pval(yorigin, permute_index, overdispersion, 10000, 0.10, 1) == 0, TRUE, FALSE) }
   # cat("mean",mean(y),"var",var(y),"estim mean",m,"overdisp",overdispersion(y,m))
    
    # null_env <- function(distri) {
    #   # simuler des données
    #   if(distri) {
    #     x <- rnbinom(length(y), size = mean(y)/(overdispersion(y) - 1), mu = mean(y)) 
    #     model_x <- glm(x ~ 1 + offset(2 * esw * l), family = "quasipoisson", control = list(epsilon = 1e-6, maxit = 10000, trace = FALSE))
    #   }
    #   else{ 
    #     x <- rpois(length(y), mean(y))
    #     model_x <- glm(x ~ 1 + offset(2 * esw * l), family = "poisson", control = list(epsilon = 1e-6, maxit = 10000, trace = FALSE))
    #   }
    #   z0 <- calc_vario(y = x, l = l, esw = esw, m=m
    #                   # m = mean(exp(rnorm(10000, as.numeric(coef(model_x)), as.numeric(sqrt(vcov(model_x))))))
    #   )
    #   vario0 <- tapply(z0, distclas, sum, na.rm = TRUE) / cstpond
    #   if(length(noise) != 0) { vario0 <- vario0[-noise] }
    #   return(vario0)
    # }
    
    null_env2<-function(b){
      set.seed(b)
      shufflebyTransect=lapply(unique(id_transect), function(x){
        indices = which( id_transect==x)
        return(sample(y[indices],length(indices), replace=FALSE))
      })
      x=unlist(shufflebyTransect)
   #   mnew=TauxObs(x,esw,l)
      mnew=mean(x/(2*esw*l))
      z0 <- calc_vario(y = x, l = l, esw = esw, m=mnew )
      vario0 <- tapply(z0, distclas, sum, na.rm = TRUE) / cstpond
      if(length(noise) != 0) { vario0 <- vario0[-noise] }
      return(vario0)
    }
    null<-do.call(rbind, lapply(1:n_sim, function(x){null_env2(x)}))
    # null <- t(replicate(n_sim, null_env(distri = negbin), simplify = "array"))
    
    # resumer les simulations et rajouter les données sur la dernière ligne
    sumnull <- rbind(apply(null, 2, quantile, probs = c(0.1, 0.5, 0.9), na.rm = TRUE), 
                     as.numeric(vario1)
    )
    sumnull= sumnull %>% t() %>%data.frame() %>%  mutate(distcl=as.character(as.numeric(distvario1)))
    colnames(sumnull)<- c("lower", "median", "upper", "emp","distcl")
    
    # reprendre les simulations
    null <- as.data.frame(null) ; names(null) <- as.character(as.numeric(distvario1))
    null$sim <- as.character(1:n_sim)
    null <- melt(null, id = "sim") ; null$variable <- as.numeric(distvario1)[as.numeric(null$variable)]
    
    # graphique avec ggplot2
    theme_set(theme_bw(base_size = 14))
    g <- ggplot() + 
      # geom_line(data = null,
      #           aes(x = variable, y = value, group = sim), 
      #           alpha = 0.05
      # ) +
      geom_ribbon(data = sumnull,
                  aes(x = as.numeric(distcl), ymin = lower, ymax = upper), 
                  fill = "midnightblue", 
                  alpha = 0.3, linetype = "dashed", size = 1
      ) +
      geom_line(data = sumnull,
                aes(x = as.numeric(distcl), y = emp), 
                color = "firebrick1", size = 1
      ) +
      labs(x="Distance",y="Semi-variance", title=paste0(esp," in ", region)) +
      theme(plot.title = element_text(lineheight = 0.8, face = "bold"), 
            axis.text = element_text(size = 12)
      ) +coord_cartesian(ylim=c(0, max(sumnull$emp)))
    
    return(list(vario = vario, g = g))
  }
  else { return(list(vario = vario)) }
}

### fonction pour prédire
ordinary_kriging <- function(y, esw, l, model, distmat_x, distmat_xy, fast_inversion = TRUE, saturate = TRUE) {
  # distmat_x = matrice des distances entre points du jeu de données
  # distmat_xy = matrice des distances entre points du jeu de données et points à prédire
  
  ### fonction pour calculer le taux de rencontre moyen
  TauxObs <- function(y, esw, l){
    my_control <- list(epsilon = 1e-6, maxit = 100000, trace = FALSE)
    # utiliser quasi poisson pour eventuelle surdispersion
    model_m <- glm(y ~ 1 + offset(2 * esw * l), family = "quasipoisson", control = my_control)
    m <- mean(exp(rnorm(1000, as.numeric(coef(model_m)), as.numeric(sqrt(vcov(model_m))))))
    return(round(m, 6))
  }
  
  ### covariance functions  
  cov_matern <- function(h, sill, range) {
    return(sill *(1 + h * sqrt(3) / range) * exp(- h * sqrt(3) / range))
  }
  
  cov_exp <- function(h, sill, range) {
    return(sill * exp(- h / range))
  }
  
  cov_cauchy <- function(h, sill, range) {
    return(sill /(1 + h / range)^2)
  }
  
  cov_spherical <- function(h, sill, range) {
    return(sill * ifelse(h < range, 1 - 3 / 2 * h / range + 1 / 2 *(h / range)^3, 0))
  }
  
  # calcul du taux d'obs moyen m
  m <- TauxObs(y, esw, l)
  
  # on prend le modele de covariance associe au model
  form <- model$cov.model
  
  # on prend l'ajustement 2 avec parametres: 
  sill <- as.numeric(model$cov.pars[1])
  range <- as.numeric(model$cov.pars[2])
  
  # nombre de points(segments)
  nbpts <- length(y)
  
  # krigeage
  
  if(form == "matern") {    
    # matrice des covariances entre toutes les données
    a <- cov_matern(distmat_x, sill, range)
    
    # matrice des covariances entre données et prédictions
    cmat <- cov_matern(distmat_xy, sill = sill, range = range)
    
  }
  
  if(form == "exponential") {
    # matrice des covariances entre toutes les données
    a <- cov_exp(distmat_x, sill, range)
    
    # matrice des covariances entre données et prédictions
    cmat <- cov_exp(distmat_xy, sill = sill, range = range)
    
  }
  
  if(form == "cauchy") {    
    # matrice des covariances entre toutes les données
    a <- cov_cauchy(distmat_x, sill, range)
    
    # matrice des covariances entre données et prédictions
    cmat <- cov_cauchy(distmat_xy, sill = sill, range = range)
    
  }
  
  if(form == "spherical") {
    # matrice des covariances entre toutes les données
    a <- cov_spherical(distmat_x, sill, range)
    
    # matrice des covariances entre données et prédictions
    cmat <- cov_spherical(distmat_xy, sill = sill, range = range)
    fast_inversion <- FALSE
  }
  
  if(fast_inversion) {
    diag(a) <- diag(a) + m /(l * 2 * esw)
    
    A1 <- chol2inv(chol(a))
    B1 <- matrix(rep(1, nrow(a)), ncol = 1)
    B2 <- A1 %*% B1
    A4 <- -1 / t(B2) %*% B1
    A2 <- B2 %*% A4
    A1 <- A1 + A2 %*% t(B2)
    
    inv_a <- rbind(cbind(A1, -A2), cbind(-t(A2), A4))
  }
  
  else {
    # construire la matrice de covariance + une colonne et une ligne de 1
    a <- rbind(cbind(a, rep(1, nbpts)), c(rep(1, nbpts), 0))
    
    #ajouter la constante sur la diagonale
    diag(a) <- diag(a) + c(m /(l * 2 * esw), 0) #on n'oublie pas le 0
    
    # ecrire inverse de a
    inv_a <- solve(a)
  }
  
  # equation de krigeage
  lambda <- cbind(cmat, rep(1, nrow(distmat_xy))) %*% inv_a
  
  # moyenne
  mean_pred <- as.numeric(lambda %*% matrix(c(y /(l * 2 * esw), 0), ncol = 1))
  mean_pred <- ifelse(mean_pred < 0, 0, mean_pred)
  if(saturate) {
    mean_pred <- ifelse(mean_pred > quantile(mean_pred, probs = 0.999), 
                        quantile(mean_pred, probs = 0.999), 
                        mean_pred
    )
  }
  # variance
  var_pred <- as.numeric(sill - diag(cmat %*% t(lambda[, 1:nbpts])) + lambda[, nbpts + 1])
  var_pred <- ifelse(var_pred < 0, 0, var_pred)
  return(data.frame(mean_pred = round(mean_pred, 3),
                    var_pred = round(var_pred + model$nugget, 3)
  )
  )
}

### faire un plot de la semi-variance
semivarPlot <- function(variogram, model, distance) {
  
  semivario_function <- function(h, sill, range, form) {
    if(form == "matern") { f <- sill *(1 -(1 + h * sqrt(3) / range) * exp(- h * sqrt(3) / range)) }
    if(form == "exponential") { f <- sill *(1 - exp(- h / range)) }
    if(form == "cauchy") { f <- sill *(1 - 1 /(1 + h / range)^2) }
    if(form == "spherical") { f <- sill *(1 - ifelse(h < range, 1 - 3 / 2 * h / range + 1 / 2 *(h / range)^3, 0)) }
    return(f)
  }
  
  obs_df <- with(variogram, data.frame(x = u, y = v, n = n))
  pred_df <- data.frame(x = distance,
                        y = model$nugget +
                          semivario_function(distance,
                                             model$cov.pars[1],
                                             model$cov.pars[2],
                                             form = model$cov.model
                          )
  )
  theme_set(theme_bw(base_size = 14))
  g <- ggplot() + 
    geom_point(data = obs_df,
               aes(x = x, y = y, size = n),
               color = "red"
    ) +
    geom_line(data = pred_df,
              aes(x = x, y = y), 
              color = "midnightblue", size = 1
    ) +
    xlab("Distance") + ylab("Semi-variance") +
    scale_size(trans = "log") +
    guides(size = "none") +
    theme(plot.title = element_text(lineheight = 0.8, face = "bold"), 
          axis.text = element_text(size = 12)
    )
  return(g)
}

fit_variomodel <- function(variogram, form) {
  return(geoR::variofit(variogram, 
                        cov.model = form, 
                        fix.nugget = FALSE, 
                        nugget = mean(variogram$v) / 4, 
                        fix.kappa = TRUE, kappa = ifelse(form == "matern", 1.5, 1),
                        weights = "npairs", 
                        minimisation.function = "optim",
                        ini.cov.pars = expand.grid(seq(0, mean(variogram$v)*2, l = 10),
                                                   seq(as.numeric(quantile(variogram$u, prob = 0.1)), 
                                                       as.numeric(quantile(variogram$u, prob = 0.8)), 
                                                       l = 10
                                                   )
                        )
  )
  )
}

### Gower's distance
# function to automate somewhat computations over 
# all possible combinations of variables
make_cfact <- function(calibration_data, test_data, var_name = NULL) {
  if(is.null(var_name)) { var_name = names(calibration_data) }
  ## standardize new data to predict from
  ### useful functions
  rescale <- function(x) { return((x - mean(x)) / sd(x)) }
  rescale2 <- function(ynew, y) { return((ynew - mean(y, na.rm = TRUE)) /(sd(y, na.rm = TRUE))) }
  # this simplifies computation A LOT!
  make_X <- function(calibration_data, test_data, var_name){
    X <- sapply(var_name,
                function(k) { rescale2(ynew = test_data[, k], y = calibration_data[, k])}
    )
    X <- as.data.frame(X); names(X) <- var_name
    return(X)
  }
  # compute counterfactuals
  cfact <- WhatIf::whatif(formula = NULL,
                          data = make_X(calibration_data = calibration_data, test_data = calibration_data, var_name),
                          cfact = make_X(calibration_data = calibration_data, test_data = test_data, var_name),
                          return.distance = TRUE
  )
  return(cfact$dist)
}

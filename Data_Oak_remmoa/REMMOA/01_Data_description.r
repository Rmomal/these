# nettoyer l'environnement de travail
rm(list = ls())

setwd(WorkDir <- "/home/momal/Git/these/Data/REMMOA")
DataDir <- paste(WorkDir, "output", sep = "/")

load(paste(DataDir, "20180620_PLN_REMMOANC.RData", sep = "/"))

### covariables
str(X)

### matrice des detections
str(N)
apply(N, 2, sum)

### matrice des comptages
str(Y)
apply(Y, 2, sum)

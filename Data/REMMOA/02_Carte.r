##--------------------------------------------------------------------------------------------------------
## SCRIPT : Carte Nouvelle Cal?donie
##
## Authors : Matthieu Authier
## Last update : 2018-06-20
## R version 3.4.3 (2017-11-30) -- "Kite-Eating Tree"
##--------------------------------------------------------------------------------------------------------

lapply(c("rgdal","maptools", "fields", "dplyr", "reshape", "ggplot2", "ggthemes", "broom", "rgdal", "readxl", "lubridate"), library, character.only=TRUE)

### clean up
rm(list = ls())
setwd(WorkDir <- "/home/momal/Git/these/Data/REMMOA")
ShapeDir <- paste(WorkDir, "shape", sep = "/")
OutDir <- paste(WorkDir, "output", sep = "/")

load(paste(OutDir, "20180620_PLN_REMMOANC.RData", sep = "/"))

### polygone de la zone d'?tude
poly_NC <- readOGR(dsn = paste(ShapeDir, "SecteurNC_UTM58.shp", sep = "/"), layer = "SecteurNC_UTM58") 

### repr?sentation de l'effort d'?chantillonnage
theme_set(theme_bw(base_size = 12))
effort_plot <- ggplot(data = X, 
                      aes(x = X, y = Y, group = Transect.Label)
                      ) + 
  geom_polygon(data = tidy(poly_NC), 
               aes(x = long, y = lat, group = group), 
               fill = grey(1.0), color = grey(0.0), size = 0.5
               ) +
  coord_equal() + 
  geom_line() + 
  xlab("Eastings") + ylab("Northings") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right"
        )
effort_plot

### visualiser les d?tections de petits delphinin?s
effort_plot +
  geom_point(data = subset(cbind(X, N), SMADEL != 0),
             aes(x = X, y = Y, size = SMADEL)
             ) +
  scale_size(name = "Nb D?tections", breaks = 1:10) +
  ggtitle("Petits Delphinin?s")

### visualiser les d?tections de P?trels/Puffins bruns
effort_plot +
  geom_point(data = subset(cbind(X, N), BROPET != 0),
             aes(x = X, y = Y, size = BROPET)
             ) +
  scale_size(name = "Nb D?tections", breaks = c(1, seq(5, 35, 5)), trans = "sqrt") +
  ggtitle("P?trels & Puffins Bruns")
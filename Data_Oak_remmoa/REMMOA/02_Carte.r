##--------------------------------------------------------------------------------------------------------
## SCRIPT : Carte Nouvelle Cal?donie
##
## Authors : Matthieu Authier
## Last update : 2018-06-20
## R version 3.4.3 (2017-11-30) -- "Kite-Eating Tree"
##--------------------------------------------------------------------------------------------------------

lapply(c("rgdal","maptools", "dplyr", "reshape", "ggplot2", "ggthemes",
         "broom", "rgdal", "readxl", "lubridate"), library, character.only=TRUE)

### clean up
rm(list = ls())
setwd(WorkDir <- "/Users/raphaellemomal/these/Data_Oak_remmoa/REMMOA")
ShapeDir <- paste(WorkDir, "shape", sep = "/")
OutDir <- paste(WorkDir, "output", sep = "/")

load(paste(OutDir, "20180620_PLN_REMMOANC.RData", sep = "/"))

### polygone de la zone d'etude
poly_NC <- readOGR(dsn = paste(ShapeDir, "SecteurNC_UTM58.shp", sep = "/"), layer = "SecteurNC_UTM58") 

### representation de l'effort d'echantillonnage
theme_set(theme_bw(base_size = 12))
effort_plot <- ggplot(data = X, 
                      aes(x = X, y = Y, group = Transect.Label)
                      ) + 
  geom_polygon(data = tidy(poly_NC), 
               aes(x = long, y = lat, group = group), 
               fill = grey(1.0), color = "darkblue", size = 0.5
               ) +
  coord_equal() + 
  geom_line(color="gray70") + 
  xlab("Eastings") + ylab("Northings") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right"
        )
effort_plot

### visualiser les detections (N) de petits delphinin?s
map_species<-function(species,mat, count=FALSE){ #SMADEL
  legend=ifelse(!count, "Nb Detections", "Obs count")
  species=enquo(species)
  data=cbind(X, mat) %>% filter(!!species!=0)
  effort_plot +
  geom_point(data = data,
             aes(x = X, y = Y, size = !!species, color = !!species)
             ) +
  scale_size(name = legend) + 
    scale_colour_gradient(name = legend,low="blue",high="firebrick1")+#,low="firebrick1",high="darkorchid4")+
  guides(size=FALSE)+
  ggtitle(species)
}
g1<-map_species(BROTER,Y, count=TRUE)
g2<-map_species(BROTER,N)
grid.arrange(g1,g2,nrow=1, ncol=2)
map_species(ZIPHIUS,Y)


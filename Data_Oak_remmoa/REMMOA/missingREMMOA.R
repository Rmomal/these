###############@###########@#######@
# data REMMOA
load(file="/Users/raphaellemomal/these/Data_Oak_remmoa/REMMOA/output/20180620_PLN_REMMOANC.RData")

colnames(X)[7:8]<-c("Xcart","Ycart")
mat_effort=outer(X$Effort,rep(1,32))
X[,c(6:10,14:20)]<-apply(X[,c(6:10,14:20)], 2, function(x) as.numeric(x))
esw=c(250,285,340,300,310,470,270,310,260,170,430,215,210,170,340,230,240,210,145,165,160,200,200,200,200,200,200,200,200,200,200,200)
species=colnames(Y)
detec=data.frame(esw=as.numeric(esw)/1e3,species)
mat_offset=outer(X$Effort, 2*as.numeric(detec$esw))
###########@###########@###########@
# missing actor

get_model<-function(data, vec, O){
  t1<-Sys.time()
  string<-paste(deparse(substitute(data)), paste(vec, collapse=" + "), sep=" ~ ")
  formula<-as.formula(string)
  mat = as.matrix(lm(formula, x=T)$x)

  model<-PLN(data ~ -1+mat + offset(log(O)))#, control=list("ftol_rel" = 1e-6, "ftol_abs"=0,"xtol_rel" = 1e-4        "xtol_abs"=1e-4,"lower_bound" = 1e-4))
  t2<-Sys.time()
  print(difftime(t2,t1))
  return(model)
}

mat_offset=outer(X$Effort, 2*as.numeric(detec$esw))
m0<-get_model(as.matrix(Y),"1", mat_offset) # plus offset

sigmaREMMOA=m0$model_par$Sigma

sigmaREMMOA[sigmaREMMOA<0]<-0 # pb d'estimation

missingREMMOA=litree(1,sigmaREMMOA,nrow(Y))
image(sigmaREMMOA)
condEsp=getCondEsp(missingREMMOA, 33)
Xobs = m0$var_par$M
Xpred = as.vector(condEsp%*%t(Xobs)) # pb couche latente
###########@###########@###########@
# conditonal expectation on map
lapply(c("rgdal","maptools", "dplyr", "reshape", "ggplot2", "ggthemes",
         "broom", "rgdal", "readxl", "lubridate"), library, character.only=TRUE)

WorkDir <- "/Users/raphaellemomal/these/Data_Oak_remmoa/REMMOA"
ShapeDir <- paste(WorkDir, "shape", sep = "/")
poly_NC <- readOGR(dsn = paste(ShapeDir, "SecteurNC_UTM58.shp", sep = "/"), layer = "SecteurNC_UTM58")

### representation de l'effort d'echantillonnage
theme_set(theme_bw(base_size = 12))
effort_plot <- ggplot(data = X,
                      aes(x = Xcart, y = Ycart, group = Transect.Label)
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

map_species<-function(species,mat,titre = "Obs count", low="black", high="firebrick1", filter=0.2, vecBounds=NULL){ #SMADEL

  species=enquo(species)
  if(is.null(vecBounds)){
    dataplot=cbind(X, mat) %>% filter(!!species>filter)
  }else{
    Ylow=vecBounds[1] ;  Yhigh=vecBounds[2] ;     Xlow=vecBounds[3];    Xhigh=vecBounds[4];
    data=cbind(X, mat) %>% filter( Ycart>Ylow, Ycart < Yhigh,
                                  Xcart>Xlow, Xcart<Xhigh)%>% filter(!!species>filter)
  }
  # browser()
  legend=titre
  effort_plot +
    geom_point(data = dataplot,
               aes(x = Xcart, y = Ycart, size = !!species, color = !!species)
    ) +
    scale_size(name = legend, range=c(0,2)) +
    scale_colour_gradient2(name =legend,low=low,high=high,mid="gray")+#,low="firebrick1",high="darkorchid4")+
    #   scale_color_brewer(palette = "YlOrRd")
    guides(size=FALSE)+
    ggtitle(species)
}
missingNode=tibble(miss=Xpred)
g1=map_species(miss,missingNode, filter=-0.5, titre="condEsp<-0.5")

g2= map_species(miss,missingNode, filter=0.2, titre="condEsp>0.2")

grid.arrange(g1,g2,nrow=1,ncol=2)

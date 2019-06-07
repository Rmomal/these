library(sp)
data(meuse)
# no trend:
coordinates(meuse) = ~x+y
data(meuse.grid)

gridded(meuse) <- TRUE
class(meuse)
bbox(meuse)
variogram(zinc~1, meuse)
# residual variogram w.r.t. a linear trend:
variogram(log(zinc)~x+y, meuse)
# directional variogram:
vv=variogram(log(zinc)~x+y, meuse, alpha=c(0,45,90,135))
variogram(log(zinc)~1, meuse, width=90, cutoff=1300)

# GLS residual variogram:
v = variogram(log(zinc)~x+y, meuse)
v.fit = fit.variogram(v, vgm(1, "Sph", 700, 1))
v.fit
plot(v.fit, cutoff=2000)
set = list(gls=1)
v
g = gstat(NULL, "log-zinc", log(zinc)~x+y, meuse, model=v.fit, set = set)
vv=variogram(g)

ggplot(vv, aes(x=dist, y=gamma))+geom_point(col="#00A1D5FF")+geom_line(col="#00A1D5FF")+theme_minimal()

##########
# test
#########
data_from_scratch("cluster", p=40, r=10, covar=NULL)

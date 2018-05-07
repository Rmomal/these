# tests poilog
library(poilog)

data.dir = '/home/momal/Git/these/Data/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))
theme_set(theme_bw())
# Parms
Y = as.matrix(Data$count);
O = Data$offset; X = Data$covariates

# Selection des especes
YO = Y/O
Rank = rank(colSums(YO))
plot(cumsum(sort(colSums(YO))))
Seuil = 20
Ycum = colSums(Y); Order = order(Ycum)
#plot(cumsum(Ycum[Order]), col = 1+(Rank[Order]<Seuil))
Y = Y[, Rank > Seuil]; O = O[, Rank >Seuil];  n = nrow(Y); p = ncol(Y)

poilogMLE(Y[1,])
bipoilogMLE(Y[,1:2])


#####################
### simulate observations

n <- rpoilog(S=80,mu=1,sig=2)

### obtain estimates of parameters
est <- poilogMLE(n)

### similar, but now with bootstrapping ###
## Not run: est <- poilogMLE(n,nboot=10)

### change start values and request tracing information 
### from optimization procedure
est <- poilogMLE(n,startVals=c(2,3),
                 control=list(maxit=1000,trace=1, REPORT=1))

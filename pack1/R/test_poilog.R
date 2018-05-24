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

Yj<-poilogMLE(Y[,1])
Yk<-poilogMLE(Y[,2])
Yjk<-bipoilogMLE(Y[,1:2])
n<-0:nrow(Y)
log(dpoilog(n,mu=Yjk$par)

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


#####################
### plot density for given parameters 
barplot(dpoilog(n=0:20,mu=2,sig=1),names.arg=0:20)

### draw random deviates from a community of 50 species 
rpoilog(S=50,mu=2,sig=1)

### draw random deviates including zeros 
rpoilog(S=50,mu=2,sig=1,keep0=TRUE)

### draw random deviates with sampling intensity = 0.5 
rpoilog(S=50,mu=2,sig=1,nu=0.5)

### how many species are likely to be observed 
### (given S,mu,sig2 and nu)? 
hist(replicate(1000,length(rpoilog(S=30,mu=0,sig=3,nu=0.7))))

### how many individuals are likely to be observed
### (given S,mu,sig2 and nu)? 
hist(replicate(1000,sum(rpoilog(S=30,mu=0,sig=3,nu=0.7))))

###########################
### change in density of n2 for two different values of rho (given n1=10)   
barplot(rbind(dbipoilog(n1=rep(10,21),n2=0:20,mu1=0,mu2=0,sig=1,sig2=1,rho=0.0),
              dbipoilog(n1=rep(10,21),n2=0:20,mu1=0,mu2=0,sig=1,sig2=1,rho=0.8)),
        beside=TRUE,space=c(0,0.2),names.arg=0:20,xlab="n2",col=1:2)
legend(35,0.0012,c("rho=0","rho=0.8"),fill=1:2)


### draw random deviates from a community of 50 species 
rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7)

### draw random deviates including zeros
rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7,keep0=TRUE)

### draw random deviates with sampling intensities nu1=0.5 and nu2=0.7 
rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7,nu1=0.5,nu2=0.7)

### draw random deviates conditioned on a certain number of species 
rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7,nu1=0.5,nu2=0.7,condS=TRUE)


### how many species are likely to be observed in at least one of the samples
### (given S,mu1,mu2,sig1,sig2,rho)? 
hist(replicate(1000,nrow(rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7))),
     main="", xlab = "Number of species observed in at least one of the samples")

### how many individuals are likely to be observed 
### (given S,mu1,mu2,sig1,sig2,rho)? 
hist(replicate(1000,sum(rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7))),
     main="", xlab="sum nr of individuals in both samples")


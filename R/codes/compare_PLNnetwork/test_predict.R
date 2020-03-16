
library(ade4)
library(PLNmodels)
data("baran95")

# data
sample.fit=sample(1:nrow(Y),round(0.8*nrow(Y)))
sample.test=setdiff(1:nrow(Y), sample.fit)
baran=prepare_data(baran95$fau[sample.fit,], baran95$plan[sample.fit,])

# fit model
model<-PLN(Abundance~date, data=baran)
 
# abundances prediction
model$predict(newdata=data.frame(baran95$plan[sample.test,"date"]),type="response") 

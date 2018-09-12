# rm(list=ls()); 
library(mclust)
par(pch=20); 

# Thresholding edge proba from PLN-tree

Pedge = readRDS('../Data_SR/scores_PLN.rds')
p = ncol(Pedge); P = p*(p-1)/2
Pvec = Pedge[upper.tri(Pedge)]
Psort = sort(Pvec); plot(Psort)

# Find elbow
par(mfrow=c(2, 1), mex=.7)
# SS = rep(Inf, P)
# for(i in (2:P)){
#   SS[i] = anova(lm(Psort[1:i] ~ 1 + as.vector(1:i)))[2, 2]
#   SS[i] = SS[i] + anova(lm(Psort[(i+1):P] ~ 1 + as.vector(1:(P-i))))[2, 2]
# }
plot(SS)
i.min = which.min(SS)
LM1 = lm(Psort[1:i.min] ~ 1 + as.vector(1:i.min))
LM2 = lm(Psort[(i.min+1):P] ~ 1 + as.vector((i.min+1):P))
plot(Psort, col=1+((1:P)>i.min))
abline(h=2/p)
abline(LM1$coef[1], LM1$coef[2], col=2)
abline(LM2$coef[1], LM2$coef[2], col=2)

# Mclust
GM = Mclust(log(Psort))
# GM = Mclust(log(Psort), modelNames='E')
print(GM)
plot(log(Psort), col=GM$classification)
colSums(GM$z)


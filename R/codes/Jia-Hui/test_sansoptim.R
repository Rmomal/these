library(EMtree)
library(PLNmodels)
library(nestor)
library(tibble)
library(ggplot2)
library(gridExtra)
# Simulated example:
set.seed(1)
pal_nodes= c("#adc9e0","#e7bd42") ; pal_edges = "#31374f"
n=200; p=15

simu=data_from_scratch("erdos",p=p,n=n)
Y=simu$data # count data
G=1*(simu$omega!=0) ; diag(G) = 0 # original dependency graph
draw_network(G, layout="kk", btw_rank=4, nodes_size=c(3, 5), 
             pal_nodes= pal_nodes, pal_edges = pal_edges, nodes_label = 1:p)$G
# PLN fit:
PLNfit=PLN(Y~1, control = list(trace=0))
PLNfit$n
# several methods for fitting EMtree:
 directFit=direct_EMtree(PLN.Cor=PLNfit,n=n, plot=TRUE, maxIter = 40, eps1 = 1e-3) #direct formula
# withMFit=optimM_EMtree(PLN.Cor=PLNfit,n=n, plot=TRUE,maxIter = 40, eps1 = 1e-6) #M as parameter of optim
optimFit=EMtree(PLN.Cor=PLNfit, plot=TRUE,maxIter = 40) #optim with M computed inside
newGrad=new_EMtree(PLN.Cor=PLNfit, plot=TRUE,n = 200)

plotPerf(newGrad$edges_prob,G,r=0,thresh =2/15)
plotPerf(directFit$edges_prob,G,r=0,thresh =2/15)
plotPerf(withMFit$edges_prob,G,r=0)
plotPerf(optimFit$edges_prob,G,r=0,thresh = 2/15)

grid.arrange(ggimage(directFit$edges_weight),ggimage(optimFit$edges_weight*100),
             ggimage(withMFit$edges_weight), ncol=3)

sd(F_Sym2Vec(directFit$edges_weight))
sd(F_Sym2Vec(withMFit$edges_weight))
sd(F_Sym2Vec(optimFit$edges_weight))
sd(F_Sym2Vec(newGrad$edges_weight))
r1=rank(F_Sym2Vec(directFit$edges_weight))
r2=rank(F_Sym2Vec(withMFit$edges_weight))
r3=rank(F_Sym2Vec(optimFit$edges_weight))
r4=rank(F_Sym2Vec(newGrad$edges_weight))
data=data.frame( r3=r3, r4=r4)
ggplot(data, aes(r4, r3))+geom_point()+theme_light()+geom_abline() +
  coord_cartesian(xlim=c(10,106), ylim=c(10,106))

hist(log(F_Sym2Vec(newGrad$edges_weight)), breaks=30)
mean(log(F_Sym2Vec(newGrad$edges_weight)))
beta=exp(log(F_Sym2Vec(directFit$edges_weight))-mean(log(F_Sym2Vec(directFit$edges_weight))))
hist(EdgeProba(F_Vec2Sym(beta)))
hist((F_Sym2Vec(directFit$edges_prob)), breaks=30)
test=log(F_Sym2Vec(directFit$edges_weight))
test[test<(-5)]=-5
hist(EdgeProba(exp(F_Vec2Sym(test))))
hist(optimFit$edges_prob)
hist(log(F_Sym2Vec(optimFit$edges_weight)), breaks=30)
plotPerf(EdgeProba(exp(F_Vec2Sym(test))),G,r=0,thresh = 2/15)

rank(F_Sym2Vec(directFit$edges_prob))

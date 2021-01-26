# pour comparer les EMtree en proba et fréquence, en direct et optim
library(EMtree)
library(nestor)
n=200 ; p=15 ; S=200
data=data_from_scratch("erdos",p=p, n=n,signed = TRUE,dens = 0.2)
G=1*(data$omega!=0)
diag(G)=0
ggimage(G)
#---------------
# with freqs
resample=ResampleEMtree(data$data, S=S,cores = 3)
new_resample=new_ResampleEMtree(counts = data$data,S =S, cores=3)
freqs = freq_selec(resample$Pmat,0.2)
new_freqs = freq_selec(new_resample$Pmat,0.2)
G_fnew=1*(new_freqs>0.2)
G_f=1*(freqs>0.2)
plot(freqs, new_freqs)
#-----
# with probs
PLN_Y = PLNmodels::PLN(data$data~1)
Prob=EMtree(PLN.Cor=PLN_Y,verbatim=TRUE, plot=TRUE)$edges_prob
newProb=new_EMtree(PLN.Cor=PLN_Y, plot=TRUE,n=n)$edges_prob
Prob_direct=direct_EMtree(PLN.Cor=PLN_Y,plot=TRUE,n=n,eps1 = 1e-3,eps2 = 1e-3,maxIter = 50)$edges_prob
G_p=1*(Prob>2/p)
G_pnew=1*(newProb>2/p)
G_pd=1*(Prob_direct>2/p)
plot(newProb,Prob)
ggimage(newProb)
ggimage(Prob)
EM=EMtree(PLN.Cor=PLN_Y,verbatim=TRUE, plot=TRUE)
newEM=new_EMtree(PLN.Cor=PLN_Y, plot=TRUE,n=n)
directEM=direct_EMtree(PLN.Cor=PLN_Y,plot=TRUE,n=n,eps1 = 1e-3,eps2 = 1e-3) 

data.frame(L=c(EM$logpY, newEM$logpY), x=c(1:3,1:5),
           model=rep(c("ancien", "bon grad"),  c(length(EM$logpY),length(newEM$logpY)))) %>% 
  ggplot(aes(x, L , color=model))+geom_point()+geom_line()+theme_light()
plot(newProb, newEM$edges_weight*Meila(newEM$edges_weight))
plot(Prob, EM$edges_weight*Meila(EM$edges_weight))
plot(Prob_direct, directEM$edges_weight*Meila(directEM$edges_weight))
#_--------------
perf_p=ppvtpr(G_p,G,r=0)[c("PPV","TPR")]
perf_new=ppvtpr(G_pnew,G,r=0)[c("PPV","TPR")]
perf_fnew=ppvtpr(G_fnew,G,r=0)[c("PPV","TPR")]
perf_f=ppvtpr(G_f,G,r=0)[c("PPV","TPR")]
perf_pd=ppvtpr(G_pd,G,r=0)[c("PPV","TPR")]
grid.arrange(ggimage(G)+labs(title="Truth"),
             ggimage(G_p)+labs(title=paste0("PPV: ",perf_p[[1]],"\nTPR: ",perf_p[[2]])),
             ggimage(G_pd)+labs(title=paste0("PPV: ",perf_pd[[1]],"\nTPR: ",perf_pd[[2]])),
             ggimage(G_pnew)+labs(title=paste0("PPV: ",perf_new[[1]],"\nTPR: ",perf_new[[2]])),ncol=2)
grid.arrange(ggimage(G)+labs(title="Truth"),
            ggimage(G_f)+labs(title=paste0("PPV: ",perf_f[[1]],"\nTPR: ",perf_f[[2]])),
            ggimage(G_fnew)+labs(title=paste0("PPV: ",perf_fnew[[1]],"\nTPR: ",perf_fnew[[2]])),ncol=3)

hist(Prob, breaks=30)
hist(newProb, breaks=30)
hist(Prob_direct, breaks=30)
auc(pred=Prob, label=G)
auc(pred=newProb, label=G)
auc(pred=freqs, label=G)
auc(pred=Prob_direct, label=G)
# summary(freqs[G==0])
# summary(Prob[G==0])

grid.arrange(draw_network(G,layout="nicely",btw_rank = 4, nodes_label = 1:15, pal_nodes = c("white","orange"))$G,
             draw_network(Prob,layout="nicely",btw_rank = 4,nodes_label = 1:15, pal_nodes = c("white","orange"))$G,
             draw_network(newProb,layout="nicely",btw_rank = 4,nodes_label = 1:15, pal_nodes = c("white","orange"))$G,
             draw_network(freqs,layout="nicely",btw_rank = 4,nodes_label = 1:15, pal_nodes = c("white","orange"))$G,
             draw_network(Prob_direct,layout="nicely",btw_rank = 4,nodes_label = 1:15, pal_nodes = c("white","orange"))$G,ncol=2)
grid.arrange(draw_network(G,layout="nicely",btw_rank = 4,nodes_label = 1:15, pal_nodes = c("white","orange"))$G,
             draw_network(G_f,layout="nicely",btw_rank = 4,nodes_label = 1:15, pal_nodes = c("white","orange"))$G,
             draw_network(G_fnew,layout="nicely",btw_rank = 4,nodes_label = 1:15, pal_nodes = c("white","orange"))$G ,ncol=2)

#---------------
# remarques: 
#1.
#le calcul direct produit des probabilités moins contrastées 
# donc plus dur de ranger les arêtes.
#Serait mieux pour trouver le graphe, mais la structure dérive vers celle d'un arbre
# assez clairement. La précision est moins bonne qu'ne seuillant els proba
# d'optim ou les fréquences.
# Il faut comprendre ce que fait optim, et pourquoi il semble que les poids des 
# arêtes ne dérivent pas vers les trop petites valeurs.
# Pourquoi il est plus intéressant de choisir une méthode qui donne une vraisemblance
# plus faible (environ 70 pour le calcul direct, environ 25 avec optim !)
#2.
# le calcule de la beetweenness doit vraiment se faire à partir des proba non seuillées
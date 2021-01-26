# source Jia_version.R
library(EMtree)
library(nestor)
library(tidyverse)
auc<-function(pred,label, digits=2){
  prediction<-ROCR::prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  ROC_auc <- round(ROCR::performance(prediction,"auc")@y.values[[1]],digits=digits)
  return(ROC_auc)
}
ppvtpr<-function(probs,G,r, thresh=0.5, digits=2){
  q=ncol(probs) ;
  PPV=round(sum((G!=0)*(probs>thresh))/(sum((G!=0)*(probs>thresh))+ sum((G==0)*(probs>thresh))),digits)#TP/(TP+FP)
  TPR=round(sum((G!=0)*(probs>thresh))/sum(G!=0), digits)
  if(r>0){
    h=(q-r):q
    PPVH=round(sum((G[h,]!=0)*(probs[h,]>thresh))/(sum((G[h,]!=0)*(probs[h,]>thresh))+ sum((G[h,]==0)*(probs[h,]>thresh))),digits)
    PPVO=round(sum((G[-h,-h]!=0)*(probs[-h,-h]>thresh))/(sum((G[-h,-h]!=0)*(probs[-h,-h]>thresh))+  sum((G[-h,-h]==0)*(probs[-h,-h]>thresh))),digits)
    TPRH=round(sum((G[h,]!=0)*(probs[h,]>thresh))/sum(G[h,]!=0), digits)
    TPRO=round(sum((G[-h,-h]!=0)*(probs[-h,-h]>thresh))/sum(G[-h,-h]!=0), digits)
    P=(G[h,]!=0) ; N=(G[h,]==0)
    FP=sum((probs[h,]>thresh)*(G[h,]==0))
    TN=sum((probs[h,]<thresh)*(G[h,]==0))
    FPRH=FP/(FP+TN)
    FNRH=1-TPRH
  }else{
    PPVO=PPV ; TPRO=TPR
    PPVH=TPRH=FPRH=FNRH=NULL
  }
  return(list(PPV=PPV,PPVH=PPVH,PPVO=PPVO,TPR=TPR,TPRH=TPRH,TPRO=TPRO,
              FPRH=FPRH,FNRH=FNRH))
}
n=200 ; p=15 ; S=200
data=data_from_scratch("cluster",p=p, n=n,signed = TRUE,dens = 0.2)
G=1*(data$omega!=0)
diag(G)=0
PLN_Y = PLNmodels::PLN(data$data~1)
EMs=list()
vrais=list()
lapply(1:10, function(i){
  EM=new_EMtree(PLN.Cor=PLN_Y, plot=FALSE,n=n, unif=FALSE)
  L=EM$logpY
  EMs[[i]]<<-EM
 vrais[[i]]<<-(data.frame(L=L, num=i, iteration=1:length(L)))
})

unif_init=new_EMtree(PLN.Cor=PLN_Y, plot=FALSE,n=n, unif=TRUE)
directEM=direct_EMtree(PLN.Cor=PLN_Y,plot=FALSE,n=n,eps1 = 1e-3,eps2 = 1e-3,maxIter = 50)

EMs[[11]]=unif_init
EMs[[12]]=directEM
data=do.call(rbind,vrais)
data=rbind(data, data.frame(L=unif_init$logpY, num=11, iteration=1:unif_init$maxIter))
data=rbind(data, data.frame(L=directEM$logpY, num=12, iteration=1:directEM$maxIter))
data %>% ggplot(aes(iteration, L, color=as.factor(num)))+geom_point()+geom_line()+theme_light()+
  guides(color=FALSE)+scale_color_manual(values=c(rep("black",10),"red", "blue"))

results=data.frame(do.call(rbind,lapply(EMs, function(EM){
  P=EM$edges_prob ; B=EM$edges_weight ; M=Meila(B)
  model=summary(lm(F_Sym2Vec(P)~F_Sym2Vec(B*M)))
  perf=ppvtpr(1*(P>2/p),G,r=0, digits=8)[c("PPV","TPR")]
  c(L=tail(EM$logpY,1),auc=auc(pred=P, label=G,8),PPV=perf$PPV, TPR=perf$TPR, R2=model$r.squared, 
    logB=log(SumTree(B)))
}))) 

results%>% mutate(unif=c(rep("noise",10),"unif","direct")) %>%  ggplot(aes(L,logB, color=unif))+geom_point()+theme_light()
results%>% mutate(unif=c(rep("noise",10),"unif","direct")) %>%  ggplot(aes(L,auc, color=unif))+geom_point()+theme_light()
results%>% mutate(unif=c(rep("noise",10),"unif","direct")) %>%  ggplot(aes(L,R2, color=unif))+geom_point()+theme_light()
results%>% mutate(unif=c(rep("noise",10),"unif","direct")) %>%  ggplot(aes(auc,R2, color=unif))+geom_point()+theme_light()
results%>% mutate(unif=c(rep("noise",10),"unif","direct")) %>%  ggplot(aes(L,PPV, color=unif))+geom_point()+theme_light()
results%>% mutate(unif=c(rep("noise",10),"unif","direct")) %>%  ggplot(aes(L,TPR, color=unif))+geom_point()+theme_light()
grid.arrange(ggimage(1*(EMs[[12]]$edges_prob>0.5)),
ggimage(1*(EMs[[11]]$edges_prob>0.5)), ncol=2)
plotPerf(1*(EMs[[12]]$edges_prob>2/p),G,r=0)
plotPerf(1*(EMs[[11]]$edges_prob>2/p),G,r=0)
best=which.max(results$L[1:11])
grid.arrange(draw_network(1*(EMs[[12]]$edges_prob>2/p), layout="nicely", btw_rank = 3)$G,
draw_network(1*(EMs[[best]]$edges_prob>2/p), layout="nicely",btw_rank = 3)$G,
draw_network(G, layout="nicely",btw_rank = 3)$G, ncol=3)
do.call(rbind,lapply(EMs, function(EM){
  sum(EM$edges_weight)
}))

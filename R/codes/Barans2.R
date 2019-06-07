
library(PLNmodels)
library(ade4)
library(tidyverse)
source('/Users/raphaellemomal/these/pack1/R/codes/FunctionsMatVec.R')
source('/Users/raphaellemomal/these/pack1/R/codes/FunctionsTree.R')
source('/Users/raphaellemomal/these/pack1/R/codes/FunctionsInference.R')

# Data
data.dir = '/Users/raphaellemomal/these/pack1/R/Data_RM/'
data(baran95)
counts = as.matrix(baran95$fau)
covar = as_tibble(baran95$plan)

n = nrow(counts)
p = ncol(counts)

# Algo parms
f= 0.80;iter.max = 1e2;B = 5e2;VEM.fit = T;Tree.fit = T;Tree.res = T;v=0.8

# code from readme
library(rlist)

model<-PLN(counts ~ covar$site)
set.seed(3)
output<-EMtree(model,  maxIter = 10, plot=TRUE)

fillFreq<-function(s){
  models<- list("1",  "covar$site", "covar$date",  c("covar$site","covar$date"))
  res<-lapply(models, function(model){
    ResampleEMtree(counts,model, S=s, maxIter=30,cond.tol=1e-12, cores=1)
  })
  return(res)
}
Faatfreq=list()
# for(i in 1:20){
#   Faatfreq[[i]]<-lapply(c(2,5,10,20,40,80,200), function(x){cat("\n !!! nb S= ",x, " iter ",i) ;fillFreq(x)})
# }
# saveRDS(Faatfreq,"/Users/raphaellemomal/simulations/Simu/PLN.2.0/Faatfreq.rds")
Faatfreq<-readRDS("/Users/raphaellemomal/simulations/Simu/PLN.2.0/Faatfreq.rds")
lapply(10:length(Faatfreq), function(x){
  to.append<-lapply(c(50,100,150), function(y){ cat("\n !!! nb S= ",y, " iter ",x) ;fillFreq(y) })
  Faatfreq[[x]]<<-c(Faatfreq[[x]],to.append)
})

#corriger le append sur les 9 premiers elmts de l'expérience
# lapply(1:9, function(x){
#   toAppend<-Faatfreq[[x]][[length(Faatfreq[[x]])]]
#   Faatfreq[[x]]<<-Faatfreq[[x]][-length(Faatfreq[[x]])]
#   Faatfreq[[x]]<<-c(Faatfreq[[x]],toAppend)
# })
saveRDS(Faatfreq,"/Users/raphaellemomal/simulations/Simu/PLN.2.0/Faatfreq2.rds")

# repérer le niveau avec 200 reechantillons
lapply(seq_along(Faatfreq), function(x){
  nbrow<-lapply(Faatfreq[[x]], function(S){
    nrow(S[[1]]$Pmat)
  })
  names(Faatfreq[[x]])<<-unlist(nbrow)
})

####################################
Faatfreq<-readRDS("/Users/raphaellemomal/simulations/Simu/PLN.2.0/Faatfreq2.rds")

Stabdata=Faatfreq[[1]]$`100`
f=0.9
x=ncol(Stabdata[[1]]$Pmat)
p=(1+sqrt(8*x+1))/2
mat<-data.frame(freq_selec_list(Stabdata,1,p,f))
models=c("Null" , "Site"   ,   "Date"       ,   "Site+Date")
allNets<-tibble(P = list(mat), models =models )  %>%
  mutate(P=map( seq_along(P), function(x) {
    df<-freq_selec_list(Stabdata,x,p,f)
    df[lower.tri(df, diag = TRUE)]<-0
    df<-data.frame(df)
    colnames(df)<-1:ncol(df)
    df
  })) %>%
  mutate(P = map(P,~rownames_to_column(.) %>%
                   gather(key, value , -rowname))) %>%
  unnest()
allNets<-allNets[,c(3,2,1,4)]
allNets %>% group_by(models) %>% summarise(sum=sum(value))
allNetsOrder=rbind(allNets[which(allNets$models=="Null"),],
                   allNets[which(allNets$models=="Date"),],
                   allNets[which(allNets$models=="Site"),],
                   allNets[which(allNets$models=="Site+Date"),])
compar_graphs(allNetsOrder, alpha=FALSE, nb=4)
ggsave("BaransNets.png",width=11.5,height=3.7)
####################################
data<-do.call(rbind,lapply(seq_along(Faatfreq),function(x){
 
  selec200<-freq_selec_pmat(Faatfreq[[x]]$`200`,p,f)
  
  obj<-do.call(rbind, lapply(Faatfreq[[x]], function(freq){
    do.call(rbind, lapply(seq_along(freq), function(model){
      Pmat=freq[[model]]$Pmat
      if(nrow(Pmat)==1){
        selec=1*(Pmat>2/p)>f
      }else{
        selec<- 1*colMeans(1*(Pmat>2/p))>f
      }
      tpfn=c(table(selec, selec200[[model]]))
      res<-data.frame("TN"=tpfn[1],"FP"=tpfn[2],"FN"=tpfn[3],"TP"=tpfn[4],
                      Ssub=nrow(freq[[model]]$Pmat),
                      model=model)
      return(res)
    }))
  }))
  
 # obj <- obj %>% as_tibble() %>% mutate(FDR=FP/(TP+FP))
}))


data=do.call(rbind,data)

saveRDS(data,"/Users/raphaellemomal/simulations/Simu/PLN.2.0/BaransTPFNdata.rds")

plotdata<-data %>% as_tibble() %>% mutate(FDR=FP/(TP+FP), model=rename_factor(as.factor(model), 
                                                                                    "1"="Null","2"="Site",
                                                                                    "3"="Date","4"="Site+Date"))%>%
  filter(!(Ssub %in% c(2,200))) %>%
  group_by(model) %>% mutate(TPnorm=(TP-mean(TP))/sd(TP))

p1<-plotdata %>% ggplot(aes( y=TPnorm, x=as.factor(Ssub), color=model, fill=model))+
  geom_boxplot()+theme_minimal()+
#  scale_color_manual(values=c(brewer.pal(n = 9, name = "YlOrRd")[-c(1:3)],"black"))+labs(x="")+
#  guides(color=FALSE)+
  geom_hline(yintercept = 0,linetype="dashed")+
  stat_summary(aes(y=TPnorm, group=model, fill=model),fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               color="black",pch=22, position=position_dodge(.9))+
  labs(y="Standardized TP",x="S")
p2<-plotdata %>% ggplot(aes( y=TP, x=as.factor(Ssub), color=model, fill=model))+
  geom_boxplot()+theme_minimal()+
  #  scale_color_manual(values=c(brewer.pal(n = 9, name = "YlOrRd")[-c(1:3)],"black"))+labs(x="")+
  #  guides(color=FALSE)+
  #geom_hline(yintercept = 0,linetype="dashed")+
  # stat_summary(aes(y=TPnorm),fun.ymin = function(z) { quantile(z,0.25) },
  # fun.ymax = function(z) { quantile(z,0.75) },
  # fun.y = median,
  # fill="white",color="black",pch=22)+
  labs(y="TP",x="S")

grid_arrange_shared_legend(p2,p1,nrow=1,ncol=2)
ggsave("DisagRate.png",height=3,width=4)
##############
# Shuffle baran data

covar=covar %>% mutate(lvl=as.factor(paste0(date,site)))
levels(covar$lvl)<-1:length(levels(covar$lvl))
counts_shuffle=0*counts
for (bloc in unique(levels(covar$lvl))){
  indices<-which(covar$lvl==bloc)
  counts_shuffle[indices,]<-apply(counts[indices,], 2, function(x){
    sample(x)
  })
}



counts_fullshuffle=apply(counts,2,function(x){sample(x)})

s=200
infShuffleEmtree<-ResampleEMtree(counts_shuffle,c("covar$site","covar$date"), S=s, maxIter=30,cond.tol=1e-15, cores=1)
fillFreq<-function(counts){
  models<- list("1")
  res<-lapply(models, function(model){
    ResampleEMtree(counts,model, S=100, maxIter=300,cond.tol=1e-12, cores=1)
  })
  return(res)
}
shuffleEMtree<-fillFreq(counts_shuffle)
shuffleNull<-fillFreq(counts_fullshuffle)

sum(F_Sym2Vec(freq_selec(infShuffleEmtree$Pmat,p,0.9)))
frequencies_shuffle<-1*colMeans(1*(infShuffleEmtree$Pmat>2/p))
hist(frequencies_shuffle,freq=FALSE, xlim=c(0,1), breaks=30, col="cornflowerblue")
entropie((frequencies_shuffle)/sum(frequencies_shuffle)+1e-10)

covarexp<-covar[,1:2]
colnames(covarexp)<-c("X1","X2")
# infShufflegCoda <- F_Sym2Vec(gcoda(counts_shuffle, counts=T, covar=covarexp)$opt.icov>1e-16)*1
# sum(infShufflegCoda)
# 
# covarexp<-covarexp %>% mutate(X1=as.numeric(X1), X2=as.numeric(X2))
# inShuffleMInt<- F_Sym2Vec(eval_store_mint(counts_shuffle,covarexp,path)>1e-16)*1
# sum(inShuffleMInt)

m<- model.matrix(~X1+X2,covarexp)
U<-t(clr.matrix(counts_shuffle,mar=1))
model<-lm(U~m)
inShuffleSpiec<-spiec.easi(model$residuals, method='mb', lambda.min.ratio=1e-2, nlambda=30,icov.select.params=list(rep.num=B ))
spiecgraph<-F_Sym2Vec(as.matrix(inShuffleSpiec$refit[[1]]))
sum(spiecgraph)
K.score <- Reduce("+",inShuffleSpiec$est$path)
scores<- K.score / max(K.score)
entropie((scores)/sum(scores)+1e-10)

#######################################
# Frenquencies histogram and entropies
#######################################
# get frenquencies

Faatfreq<-readRDS("/Users/raphaellemomal/simulations/Simu/PLN.2.0/Faatfreq2.rds")

Freqs<-do.call(rbind,lapply(seq_along(Faatfreq),function(x){
  obj<-do.call(rbind, lapply(Faatfreq[[x]], function(freq){
    do.call(rbind, lapply(seq_along(freq), function(model){
      Pmat=freq[[model]]$Pmat
      if(nrow(Pmat)==1){
        selec=1*(Pmat>2/p)
      }else{
        selec<- 1*colMeans(1*(Pmat>2/p))
      }
      res<-data.frame(selec,Ssub=nrow(Pmat),
                      model=model, edge=1:ncol(Pmat), sample=x)
      return(res)
    }))
  }))
  return(obj)
}))

entropie<-function(vec){
  -sum(vec*log(vec+1*(vec==0)))
}
pal_edges <- viridisLite::viridis(5, option = "C")[c(2,4,1,3)]
naretes=p*(p-1)/2
plot<-Freqs %>% mutate( model=rename_factor(as.factor(model), 
                                      "1"="Null","2"="Site",
                                      "3"="Date","4"="Site+Date")) %>% 
  group_by(sample, model, Ssub) %>% count(selec) %>% mutate(proportion=n/naretes, 
                                                            bords=ifelse((selec<0.1 | selec >0.9),
                                                                           proportion, 0),
                                                            milieu=ifelse((selec>=0.1 & selec <=0.9),
                                                                         proportion, 0),
                                                            haut=ifelse((selec>0.9),
                                                                          proportion, 0)) %>% 
  summarise(med=median(selec),bords=sum(bords),milieu=sum(milieu),haut=sum(haut),H=entropie(proportion),
            mean=weighted.mean(selec, w=proportion), var=weighted.mean(selec^2, w=proportion)-mean^2) %>%# mutate(Hnorm=H/log(Ssub+1)) %>% filter(Ssub>=10) %>% 
  mutate(H=H/log(Ssub+1))  %>% filter(Ssub>5) %>% 
  ggplot(aes(x=Ssub, color=as.factor(model),pch=as.factor(model)))+scale_color_manual(values=  pal_edges )+
  labs(x="S")+theme_minimal()+guides(color=guide_legend(title=""),pch=guide_legend(title=""))

plot1<-plot+  geom_point(aes(y=bords),size=2)+labs(title="bords")
plot2<-plot+  geom_point(aes(y=milieu),size=2)+labs(title="milieu")
plot3<-plot+  geom_point(aes(y=haut),size=2)+labs(title="haut")

grid_arrange_shared_legend(plot1,plot2,plot3, nrow=3, ncol=1)

saveRDS(Freqs,"tabFreqs.rds")
# ce qu'on observe : les bords s'enrichissent et le milieu s'appauvrit plus le nombre de rééchantillons augmente.
# Cependant, c'est surtout le bord inférieur qui prend de l'ampleur, déjà à partir de 10 rééchantillons le bord supérieur contient
# presque la même proportion de déccouverte qu'avec 200 rééchantillons, d'autant plus visible pour un f'=80%.
# Néammoins les arêtes sélectionnées ne sont pas forcément les mêmes, on voit une sorte de contraintes de nombre de découvertes
# ici.

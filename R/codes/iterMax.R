parameters<-list(c(seq(10,30,2)),c(seq(10,120,10)), c(seq(0.5,5,0.5)/20),
                 c(seq(1,30,5)),c(seq(0.1,0.4,0.05)))
names(parameters)<-c("d","n","prob","r","dens")
Bgraph<-100

variable<-"n"
type<-"erdos"
path<-"/Users/raphaellemomal/simulations/Simu/PLNcov/"
library(parallel)
sapply(parameters[[variable]],function(x){
  print(paste0("type: ", type," // var: ", variable," // valeur :",x))
  mclapply(1:Bgraph,function(nbgraph){
    inf<-from_stored_graphs(type,variable,nbgraph,x) #returns list(inf_treeggm, inf_treeggmCond,K, itermax)
    save_file<-paste0(path,type,"/test100/",variable,"/Scores/")
    saveRDS(inf,paste0(save_file,"Graph",nbgraph,"_EM100_",x,".rds"))
  }, mc.cores=3)
})
# 13h39 - 14h26: 45min
itermax<-c()
sapply(parameters[[variable]],function(x){
  lapply(1:Bgraph,function(nbgraph){
    itermax<<-c(itermax,
                readRDS(paste0(paste0(path,"erdos/test100/",variable,"/Scores/"),"Graph",nbgraph,"_EM100_",x,".rds"))[[4]])
  })
})
par(mfrow=c(1,1))
hist(itermax)
summary(itermax)
p1<-ggplot(data.frame(itermax=itermax),aes(itermax))+
 geom_histogram(fill="#56B4E9", bins=35)



##### coudes
nbgraph<-15
x<-30 # valeur de n
inf<-readRDS(paste0(paste0(path,"erdos/test100/n/Scores/Graph"),nbgraph,"_EM100_",x,".rds"))
infEM<-inf[[1]]
infEMCond<-inf[[2]]
K<-inf[[3]]
inf20<-readRDS(paste0(paste0(path,"erdos/n/Scores/Graph"),nbgraph,"_treeggm_",x,".rds"))
#par(mfrow=c(1,2))
plot(sort(inf20))
points(sort(infEM), col="red")
dataggplot<-data.frame(index=1:length(inf20),inf20=sort(inf20),inf100=sort(infEM))
dataggplot<-gather(dataggplot, key, value, -index)


p2<-ggplot(dataggplot,aes(x=index,value,color=key))+
  geom_point(size=0.6)+
  geom_line(size=0.3)+
  scale_color_manual(values=c("#56B4E9","#fc5f94"))+
  labs(title=paste0("Scores with Erös graph n°",nbgraph,", value of n: ",x))+
  theme(legend.title = element_blank())


##### AUC
length<-length(parameters[[variable]])
df<-data.frame(var=rep(parameters[[variable]],Bgraph),param=rep(1:Bgraph,each=length),
              inf100=rep(0,Bgraph*length), inf20=rep(0,Bgraph*length))

sapply(parameters[[variable]],function(x){
  lapply(1:Bgraph,function(nbgraph){
    inf<-readRDS(paste0(paste0(path,"erdos/test100/n/Scores/Graph"),nbgraph,"_EM100_",x,".rds"))
    infEM<-inf[[1]]
    inf20<-readRDS(paste0(paste0(path,"erdos/n/Scores/Graph"),nbgraph,"_treeggm_",x,".rds"))
    K<-inf[[3]]
    df$inf100[which(df$var==x & df$param==nbgraph)]<<-roc_curve(infEM,K)
    df$inf20[which(df$var==x & df$param==nbgraph)]<<-roc_curve(inf20,K)
    })
})
#diagnostics_itermax<-function(data){
  tab<-df
  lignes<-which(is.na(tab[,1]))
  if (length(lignes)!=0) tab<-tab[-lignes,]
  tab<- gather(tab,key=method,value=value,inf100,inf20)
  tab<-summarise(group_by(tab,var,method),mns=median(value),inf=quantile(value,0.25),sup=quantile(value,0.75))

  tab$var<-as.numeric(as.character(tab$var))
  p3<-ggplot(tab, aes(y=mns,x=as.numeric(as.character(var)),shape=method,color=method))+
    geom_errorbar(aes(ymin=inf, ymax=sup), width=0,position=position_dodge((max(tab$var)-min(tab$var))/100))+
    geom_point(size=2.5)+
    geom_line(size=0.2)+
    labs(y=variable,x="n")+
    scale_shape_manual(values=c(16,8)
                       # ,
                       # breaks=c("treeggm","ggm1step", "spiecEasi","spiecResid" ,"oracle"),
                       # labels=c("EM ","1 step","SpiecEasi","spiecResid", "oracle" ) 
                       )+
    scale_color_manual(values=c("#56B4E9","#fc5f94")
                       # , "#8037c9","#56B4E9" ,"#E69F00")
                       # ,
                       # breaks=c("treeggm","ggm1step", "spiecEasi","spiecResid" ,"oracle"),
                       # labels=c("EM ","1 step","SpiecEasi","spiecResid", "oracle" ) 
                       )+
    scale_y_continuous(limits = c(0.5,1))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())
  #print(p)
#}
  pdf(paste0(path,"erdos/test100/n",x,".pdf"),
      width=10,
      height=5,onefile=TRUE)
  print(grid.arrange(p1,p2,p3,nrow=2,ncol=2))
  dev.off()
  


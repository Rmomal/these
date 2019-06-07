source("/Users/raphaellemomal/simulations/codes/Kit_pour_VEM_EM.R")
library(parallel)
library(tidyverse)
##############
# FUNCTIONS
##############
Make_graph<-function(valeur,variable,type,path,nbgraph){ 
  if( variable=="n"){
    graph<-generator_graph(graph=type,d=20,prob=0.25,dens=0.25)
    param<-generator_param(as.matrix(graph))
    
    file<-paste0(paste0(path,type,"/n/Sets_param/Graph",nbgraph,".rds"))
    saveRDS(param,file)
  }else{
    graph<-switch(variable, "d"=generator_graph(graph=type,d=valeur,prob=5/valeur,dens=5/valeur),
                  "prob"=generator_graph(graph=type,prob=valeur),
                  "r"=generator_graph(graph=type,r=valeur,dens=0.25),
                  "dens"=generator_graph(graph=type,dens=valeur))
    param<-generator_param(as.matrix(graph))
    file<-paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds")
    saveRDS(param,file)
  }
  
}


Make_data<-function(type, variable, covar, path, cores){
  mclapply(parameters[[variable]], function(valeur) {
    lapply(1:100, function(nbgraph) {
      data<-data_from_stored_graphs(type, variable, nbgraph, valeur, covar=covar,path)
      saveRDS(data,paste0(path,type,"/",variable,"/YData/Y_",nbgraph,"_",valeur,".rds"))
    })
  }, mc.cores = cores)
}

################################
scores_EM<-function(type, variable, covar=fixeCovar, path, step="FALSE", Bgraph, cond.alt,cores=1){
  # inferences
  
  expected_length<-switch(step,"FALSE"=6,"TRUE"=5)
  obj<-lapply(parameters[[variable]], function(valeur) {
    if(variable=="n") covar<-readRDS(paste0(path, "fixeCovar_n",valeur,".rds"))
    print(paste0("type: ", type," // var: ", variable," // valeur :",valeur))
    mclapply(1:Bgraph, function(nbgraph) {
      from_stored_graphs(type, variable, nbgraph, valeur,covar=covar,step=step,cond.tol=1e-12,path=path)
    },mc.cores=cores)
  })
  ## Check if success
  lengths_inlist<-sapply(seq_along(parameters[[variable]]), function(i){
    mean(sapply(obj[[i]],function(x){
      length(x[[1]])
    }))
  })
  success<-mean(lengths_inlist)==expected_length
  if(!success){
    # reperer dans quelle valeur, puis dans quel graph
    vals_fail<-which(lengths_inlist!=6)
    graphs_fail<-lapply(vals_fail,function(valeur){
      which(sapply(obj[[valeur]],function(x){
        length(x[[1]])
      })!=expected_length)
    })
    # redo the failing inferences, without parallelization for easier storage too
    lapply(parameters[[variable]][vals_fail], function(valeur) {
      if(variable=="n") covar<-readRDS(paste0(path, "fixeCovar_n",valeur,".rds"))
      index_val<-which(parameters[[variable]]==valeur)
      seq_graphs<-graphs_fail[[which(parameters[[variable]][vals_fail]==valeur)]]
      lapply(seq_graphs, function(nbgraph) {
        obj[[index_val]][[nbgraph]] <<-from_stored_graphs(type, variable, nbgraph, valeur,
                                                          covar=covar,step=step, cond.tol=cond.alt,path=path)
      })
    })
  }
  
  ## save scores
  lapply(seq_along(parameters[[variable]]), function(valeur) {
    lapply(1:Bgraph, function(nbgraph) {
      save_path<-paste0(path,type,"/",variable,"/Scores/")
      val<-parameters[[variable]][valeur]
      if(step=="FALSE"){
        saveRDS(obj[[valeur]][[nbgraph]][[1]]$P,paste0(save_path,"EMmarg_graph",nbgraph,"_",val,".rds"))
        saveRDS(obj[[valeur]][[nbgraph]][[1]]$probaCond,paste0(save_path,"EMCond_graph",nbgraph,"_",val,".rds"))
      }else{
        saveRDS(obj[[valeur]][[nbgraph]][[1]]$P,paste0(save_path,"OneMarg_graph",nbgraph,"_",val,".rds"))
        saveRDS(obj[[valeur]][[nbgraph]][[1]]$probaCond,paste0(save_path,"OneCond_graph",nbgraph,"_",val,".rds"))
      }
    })
  })
  ##save times alpha and maxIter
  names(obj)<-parameters[[variable]]
  tidyObj<-obj  %>% 
    as_tibble() %>%
    gather(valeur_param,listes) %>% 
    mutate(listes = map(listes, ~as_tibble(t(.x)))) %>%
    unnest()
  if(step=="FALSE"){
    tidyObj<-tidyObj %>% 
      mutate(alpha=map(V1,~.x[[5]]),maxIter=map(V1,~.x[[4]])) %>% 
      dplyr::select(-V1) %>% 
      rename(time=V2) %>% 
      unnest()
    saveRDS(tidyObj,paste0(path,type,"/",variable,"/EM_features.rds"))
  }else{
    tidyObj<-tidyObj %>% 
      mutate(alpha=map(V1,~.x[[5]])) %>% 
      dplyr::select(-V1) %>% 
      rename(time=V2) %>% 
      unnest()
    saveRDS(tidyObj,paste0(path,type,"/",variable,"/OneStep_features.rds"))
  }
}


scores_gcoda_spiec_resid<-function(type, variable, method,cores=1,covar=fixeCovar){#method ou spiecResid
  
  obj<-lapply(parameters[[variable]], function(valeur) {
    print(paste0("method: ",method, " // type: ", type," // var: ", variable," // valeur :",valeur))
    if(variable=="n") covar<-readRDS(paste0(path, "fixeCovar_n",valeur,".rds")) 
    mclapply(1:100, function(nbgraph) {
      Y<-readRDS(paste0(path,type,"/",variable,"/YData/Y_",nbgraph,"_",valeur,".rds"))[[1]]
      #imputation par la mediane, parfois NA, bizarre
      # if(sum(is.na(Y))!=0){
      #   index<- which(is.na(Y),arr.ind=TRUE)
      #   if(sum(is.na(Y))==1){
      #     Y[index]<-median(Y[,index[,2]], na.rm=TRUE)
      #   }else{
      #     Y[index]<-apply(Y[,index[,2]],2,function(x) median(x, na.rm=TRUE))
      #   }
      # }
      if(method=="gcodaResid"){
        
        T1<-Sys.time()
        
        out_gcodaResid<-gcoda(Y, counts=T, covar=covar)
        K.score <- Reduce("+",out_gcodaResid$path)
        scores<- K.score / max(K.score)
        
        T2<-Sys.time()
        time<-difftime(T2,T1)
        
      }else{
        T1<-Sys.time()
        # browser()
        U<-t(clr.matrix(Y,mar=1))
        m<- model.matrix(~X1+X2+X3,covar)
        model<-lm(U~m)
        inf<- spiec.easi(model$residuals, icov.select = FALSE, nlambda = 50, verbose = FALSE)
        K.score <- Reduce("+",inf$est$path)
        scores<- K.score / max(K.score)
        
        T2<-Sys.time()
        time<-difftime(T2,T1)
      }
      return(list(scores,time))
    },mc.cores=cores)
  })
  #save scores
  lapply(seq_along(parameters[[variable]]), function(valeur) {
    lapply(1:Bgraph, function(nbgraph) {
      save_path<-paste0(path,type,"/",variable,"/Scores/")
      val<-parameters[[variable]][valeur]
      saveRDS(obj[[valeur]][[nbgraph]][[1]],paste0(save_path,"/",method,"_graph",nbgraph,"_",val,".rds"))
    })
  })
  #save times
  names(obj)<-parameters[[variable]]
  tidyObj<-obj  %>% 
    as_tibble() %>%
    gather(valeur_param,listes) %>% 
    mutate(listes = map(listes, ~as_tibble(t(.x)))) %>%
    unnest()%>% 
    dplyr::select(-V1) %>% 
    rename(time=V2) %>% 
    unnest()
  
  saveRDS(tidyObj,paste0(path,type,"/",variable,"/",method,"_features2.rds"))
}


#######
# RUN
#######
parameters<-list(d=c(seq(10,30,2)),n=c(seq(20,120,10)), prob=c(seq(1,5,0.5)/20), 
                 r=c(seq(1,50,5)), dens=c(seq(1,5,0.5)/20))
Bgraph<-100
path<-"/Users/raphaellemomal/simulations/Simu/PLN_nonfav/"
fixeCovar<-readRDS(paste0(path, "fixeCovar.rds")) #cbind(rep(c(0,1),each=n/2),rnorm(n,8,0.5),round(runif(n)*10))
T1<-Sys.time()
type<-c("erdos","cluster")
for(type in "cluster"){
  cparam<-switch(type,"cluster"=c("n","d","r"),"erdos"=c("n","d"))
  for(variable in cparam){
    sapply(parameters[[variable]],function(valeur){
      cat("variable: ", variable," / valeur: ",valeur,"\n")
      mclapply(1:Bgraph,function(nbgraph){
        Make_graph(valeur,variable,type,path,nbgraph)
      },mc.cores=3)
    })
  }
}

T2<-Sys.time()
difftime(T2,T1)

T1<-Sys.time() # attention LONG, idem
for(type in c("erdos","cluster")){
  cparam<-switch(type,"cluster"=c("n","d","r"),"erdos"=c("n","d"))
  #  cparam<-"n"
  for( variable in cparam){
   # Make_data(type, variable, covar=fixeCovar, path, cores=1)
    scores_EM(type, variable, covar=fixeCovar, path, step="FALSE",Bgraph=Bgraph,cond.alt=1e-6, cores=3)
    scores_EM(type, variable, covar=fixeCovar, path, step="TRUE",Bgraph=Bgraph,cond.alt=1e-6, cores=3)

  }
}
T2<-Sys.time()
difftime(T2,T1)



T1<-Sys.time() #n beaucoup trop long, tester avec des valeurs faciles
for(type in c("erdos","cluster")){
  cparam<-switch(type,"cluster"=c("n","d","r"),"erdos"=c("n","d"))
  #cparam<-"n"
  for( variable in cparam){
    scores_gcoda_spiec_resid(type, variable, method="gcodaResid",cores=3,covar=fixeCovar)
  }
}
T2<-Sys.time()
time_gcoda<-difftime(T2,T1)

T1<-Sys.time() #idem
for(type in c("erdos","cluster")){
  cparam<-switch(type,"cluster"=c("n","d","r"),"erdos"=c("n","d"))
  #cparam<-"n"  
  for( variable in cparam){
    scores_gcoda_spiec_resid(type, variable, method="spiecResid",cores=3,covar=fixeCovar)
  }
}
T2<-Sys.time()
time_spiec<-difftime(T2,T1)





clr.default <- function(x.f, base=exp(1), tol=.Machine$double.eps) {
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

clr.matrix <- function(x.f, mar=2, ...) {
  apply(x.f, mar, clr, ...)
}









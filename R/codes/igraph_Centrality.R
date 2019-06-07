# Detection d'acteurs-clef

rm(list=ls())
source('/home/robin/PgmR/General/FunctionsMatVec.R')
library(sna); library(igraph); library(combinat); library(influenceR)
par(mfrow=c(1, 1), mex=.5, pch=20)

# Parms
StabselBarans <- readRDS("/Users/raphaellemomal/simulations/Illustrations/Barans/StabselBarans.rds")
StabselOak <- readRDS("/Users/raphaellemomal/simulations/Illustrations/Oak/StabselOak.rds")
dataList = list(StabselBarans,StabselOak)
modelGlobalList = list(c("Null","Site","Date","Site+Date"), c("Null","Tree","Tree+D1+D2+D3"))
keyNbList = c(3, 5)
datanames=c("Barans","Oak")
# functions
f = 0.9
freq_selec <- function(Pmat, p, f){return(F_Vec2Sym(1*(colMeans( 1*(Pmat>2/p))>f)))}

# Data 
dev.off()
for (dataNum in 1:length(dataList)){
 
  Stabsel = dataList[[dataNum]]
  dataName=datanames[dataNum]
  nbSpecies = (1+sqrt(8*ncol(Stabsel[[1]]$Pmat)+1))/2
  angle = 2*pi*(1:nbSpecies)/nbSpecies; nodePositions = cbind(cos(angle), sin(angle))
  
  # Models
  modelList = modelGlobalList[[dataNum]]; modelNb = length(modelList); 
  
  # Models
  keyNb = keyNbList[dataNum]; 
  par(mfrow=c(2, 2), mex=.5); 
  for (m in 1:modelNb){

    title = paste(dataName, modelList[m])
    cat('\n', title, '\n')
    net = freq_selec(Stabsel[[m]]$Pmat, nbSpecies, f)
    degree = colSums(net)
    igraphNet = graph.adjacency(net, mode='undirected', diag=F)
    
    # Authority
    authority = authority_score(igraphNet)
    authorityBest = order(authority$vector, decreasing=T)[1:keyNb]
    cat('Authority:', authorityBest, '\n')
    
    # Betweenness
    betweenness = betweenness(net)
    betweennessBest = order(betweenness, decreasing=T)[1:keyNb]
    cat('Betweenness:', betweennessBest, '\n')
    
    # Closeness centrality
    closeness = closeness(net)
    closenessBest = order(closeness, decreasing=T)[1:keyNb]
    cat('Closeness centrality:', order(closenessBest, decreasing=T)[1:keyNb], '\n')
    
    # Key players
    solKeyPlayers = c(); solNb = 0
    cat('Key players:')
    for(i in 1:100){#cat(i, '')
      tmpKeyPlayers = sort(keyplayer(igraphNet, k=keyNb))
      OK = F; s=0;
      while(s<solNb){s = s+1; 
      if(prod(solKeyPlayers[s, ]==tmpKeyPlayers)==1){OK=T; s=solNb}
      }
      if(!OK){
        solKeyPlayers = rbind(solKeyPlayers, tmpKeyPlayers); solNb = solNb+1
        cat(' (', solNb, ') ', tmpKeyPlayers)
        nodeColor = rep(0, nbSpecies);
        nodeColor[solKeyPlayers[solNb, ]] = 4
        nodeColor[authorityBest] = 2
        nodeColor[betweennessBest] = 3
        nodeColor[closenessBest] = 5
        gplot(net, gmode='graph', coord=nodePositions, vertex.col=nodeColor, edge.col=8,
              main=paste0(title, ' (', solNb, ')'))
      }
    }
    cat('\n')
    
    # Store
    Stabsel[[m]]$degree = degree
    Stabsel[[m]]$betweenness = betweenness
    Stabsel[[m]]$authority = authority$vector
    Stabsel[[m]]$closeness = closeness
    
    # Plot
    nodeColor = rep(0, nbSpecies); 
    nodeColor[unique(as.vector(solKeyPlayers))] = 4
    nodeColor[authorityBest] = 2
    nodeColor[betweennessBest] = 3
    nodeColor[closenessBest] = 5
    gplot(net, gmode='graph', coord=nodePositions, vertex.col=nodeColor, edge.col=8, 
          main=title)
  }
  
  par(mfcol=c(modelNb*(modelNb-1)/2, 2), mex=.6)   
  for (m in 1:(modelNb-1)){for(mm in ((m+1):modelNb)){
    plot(Stabsel[[m]]$degree, Stabsel[[mm]]$degree, xlab=modelList[m], ylab=modelList[mm])
    abline(0, 1)}}
  for (m in 1:(modelNb-1)){for(mm in ((m+1):modelNb)){
    plot(Stabsel[[m]]$authority, Stabsel[[mm]]$authority, 
         xlab=modelList[m], ylab=modelList[mm])
    abline(0, 1)}}
  # par(mfrow=c(modelNb, modelNb))   
  # for (m in 1:(modelNb-1)){for(mm in ((m+1):modelNb)){
  #    plot(Stabsel[[m]]$betweenness, Stabsel[[mm]]$betweenness, main='betweenness', xlab=modelList[m], ylab=modelList[mm])
  #    abline(0, 1)}}
  # par(mfrow=c(modelNb, modelNb))   
  # for (m in 1:(modelNb-1)){for(mm in ((m+1):modelNb)){
  #    plot(Stabsel[[m]]$closeness, Stabsel[[mm]]$closeness, main='closeness', xlab=modelList[m], ylab=modelList[mm])
  #    abline(0, 1)}}
}


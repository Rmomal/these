######################################################################################
# Others
FitGlasso <- function(S){
   p = ncol(S)
   rholist = 10^seq(-3, 3, by=.1)
   GLpath = glassopath(S, rholist=rholist, penalize.diagonal=F, trace=0)
   rhoThresh = sapply(1:p, function(i){sapply(1:p, function(j){rholist[max(which(GLpath$wi[i, j, ] != 0))]})})
   return(rhoThresh)
}
CompGraph <- function(Gtrue, Ghat){
   p = ncol(Gtrue)
   Gtot = Gtrue + Ghat; Gtot[which(Gtot>1)] = 1
   TP = (Gtrue*Ghat); FN = (Gtot*(1-Ghat)); FP = ((1-Gtot)*Ghat)
   Col = 1*TP + 2*FN + 3*FP
   Position=gplot(Gtot, gmode='graph', edge.col=Col)
   gplot(TP, gmode='graph', coord=Position, edge.col=1)
   gplot(FN, gmode='graph', coord=Position, edge.col=2)
   gplot(FP, gmode='graph', coord=Position, edge.col=3)
}
EstimGraphDens <- function(Cov, n){
   p = ncol(Cov)
   Gtot = Gtrue + Ghat; Gtot[which(Gtot>1)] = 1
   TP = (Gtrue*Ghat); FN = (Gtot*(1-Ghat)); FP = ((1-Gtot)*Ghat)
   Col = 1*TP + 2*FN + 3*FP
   Position=gplot(Gtot, gmode='graph', edge.col=Col)
   gplot(TP, gmode='graph', coord=Position, edge.col=1)
   gplot(FN, gmode='graph', coord=Position, edge.col=2)
   gplot(FP, gmode='graph', coord=Position, edge.col=3)
}

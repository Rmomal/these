source("Functions.R")
#----------------
# Simulated example:
set.seed(1)
pal_nodes= c("#adc9e0","#e7bd42") ; pal_edges = "#31374f"
n=200; p=15

simu=data_from_scratch("erdos",p=p)
Y=simu$data # count data
G=1*(simu$omega!=0) ; diag(G) = 0 # corresponding original dependency graph
draw_network(G, layout="kk", btw_rank=4, nodes_size=c(3, 5), 
             pal_nodes= pal_nodes, pal_edges = pal_edges, nodes_label = 1:p)$G
# The nodes 6, 8 and 15 show the biggests betweenness centrality scores.
unlink=c(6,8,15)

# a classic fit of new_EMtree:
linkedFit=new_ResampleEMtree(counts=Y,unlinked = NULL, S=100, cores=3,maxIter = 50)
linked_network<-1*(freq_selec(linkedFit$Pmat,Pt=0.1)>0.4)
draw_network(linked_network, layout="kk", groupes =(1:p)%in%unlink, nodes_size=c(3, 5), 
             pal_nodes= pal_nodes, pal_edges = pal_edges,nodes_label = 1:p)$G

## now let's unlink 6, 8 and 15:
unlinkedFit=new_ResampleEMtree(counts=Y,unlinked = unlink, S=100, cores=3) 
unlinked_network<-1*(freq_selec(unlinkedFit$Pmat,Pt=0.1)>0.4)
draw_network(unlinked_network, layout="kk", groupes =(1:p)%in%unlink, nodes_size=c(3, 5), 
             pal_nodes= pal_nodes, pal_edges = pal_edges, nodes_label = 1:p)$G



####### test jia-hui data
test.this.Y.data <- read.csv("~/Downloads/test this Y data.csv")
testfit=new_ResampleEMtree(counts=test.this.Y.data,unlinked = 1:10, S=100, cores=3,maxIter = 50)
testfit$Pmat

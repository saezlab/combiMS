# Run a given graph and determine whether the source and target are connected
# Script kindly provided by Federica Eduati
# Modified by Marti Bernardo-Faura, final version october 2014

# INPUT:
# graph => optimized model converted to a graph (e.g. using graph<-sif2graph(model2sif(model=PKN, optimRes=list(bString=MyAverageBstring))))
# stimulus => the name of the stimulus
# readout => the name of the inhibitor
# link => the selected link (with the formalism of CNO model, e.g. A=B, !A=B, A+B=C)
#             
# OUTPUT:
# output => logic value: TRUE if there is a connection, FALSE if not

is.connected <-function(graph, stimulus, readout, link){
  link<-gsub("!", "", link)
  link<-unlist(strsplit(link, "=", fixed=TRUE))
  EndNode<-link[2]
  # can contain more than 1 element in case of AND nodes
  StartNodes<-unlist(strsplit(link[1], "+", fixed=TRUE))
  
  # 1. check if there is a path from the stimulus to the start node(s) of the selected link
  path1<-vector()
  for (i in 1:length(StartNodes)){
    # output is the lenght of the shortest path, or NA if there in no connection
    path1[i]<-sp.between(g=graph,start=stimulus,finish=StartNodes[i], detail=FALSE)[[1]]$length
    path1[i]<-!is.na(path1[i])
    # handle case in which the start node of the link is the stimulus
    if (StartNodes[i]==stimulus){
      path1[i]<-T
    }
  }
  path1<-all(path1==T)
  
  # 2. if the first connection is there, check if the selected link is in the network
  path2<-F
  if (path1==T){
    path2<-vector()
    for (i in 1:length(StartNodes)){
      # output is the lenght of the shortest path, or NA if there in no connection
      path2[i]<-sp.between(g=graph,start=StartNodes[i],finish=EndNode, detail=FALSE)[[1]]$length
      path2[i]<-!is.na(path2[i])
    }
    path2<-all(path2==T)
  }
  
  # 3. if both the first 2 path are present, check there is also the connection from the end node of the link to the readout
  path3<-F
  if (path1==T & path2==T){
    path3<-vector()
    path3<-sp.between(g=graph,start=EndNode,finish=readout, detail=FALSE)[[1]]$length
    path3<-!is.na(path3)
    if (EndNode==readout){
      path3<-T
    }
  }
  
  # if all 3 paths are present then the output is TRUE, other
  output<-F
  if (path1==T & path2==T & path3==T){
    output<-T
  }
  
  return(output)
}


# test
# library(RBGL)
# setwd("/Users/eduati/Dropbox/forMarti")
# load("exampleGraph.RData")
# plot(exampleGraph)
# is.connected(graph=exampleGraph, stimulus='IGF1', readout='TP53', link="AKT=MDM2")
# source("is.connected.R")
# is.connected(graph=exampleGraph, stimulus='IGF1', readout='TP53', link="IGF1=IGF1R")
# is.connected(graph=exampleGraph, stimulus='IGF1', readout='TP53', link="MDM2=TP53")
# 



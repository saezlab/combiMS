# Function to discard bias of final patient models due to technical factors such as
# treatment, center, disease subtipe.

#Final version October 2014, Marti Bernardo-Faura

#apclusterPatients<-function(Models){
  
  library(apcluster)
  library(gclus)
  library(ggplot2)
  library(reshape2)
  library(NMF)
  library(corrplot)
  
  
  load("/Users/marti/Documents/R/combiMS/modelling/AvgModelsInTol.RData")
  Models=AvgModelsInTol

  # ************************************ 
  # calculate patient correaltion as similarity matrix
  # ************************************
  #Cor=cor(t(Models))
  D=dist(Models)
  Cor=as.matrix(D)
  #aheatmap(Cor)
  #corrplot(Cor,order = "hclust", addrect = 3)
  annot=read.csv("/Users/marti/Documents/R/combiMS/modelling/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
  aheatmap(Cor, annColors=c("Paired","Dark2","Accent"),distfun="euclidean", annRow=annot[,c("Disease.Subtype","condition","Treatment")])
  
  # ************************************ 
  # ap cluster using similiartiy matrix
  # ************************************ 
  clustered=apcluster(Cor)
  #cat("give e.g. preference=-30 for a very reduced number of clusters \n")
  cat("affinity propogation optimal number of clusters:", length(clustered@clusters), "\n")
  
  # ************************************ 
  # save clusters in both annotation file and dataframe for different plots
  # ************************************ 
  # ********* create a molten data frame for plotting
  df=melt(Models)
  colnames(df)=c("patient","interaction","value")
  
  # ********* asign cluster names to patients
  df$cluster=vector(length=length(df$value)) #add an empty column to the df to contain the clusters
  df$exemplar=vector(length=length(df$value)) #add an empty column to the df to contain the exemplars
  annot$cluster=vector(length=length(annot$ID))
  
  for (index in 1:length(clustered@clusters)){  
    df$cluster[which(df$patient %in% names(clustered@clusters[[index]])==TRUE)]=index
    df$exemplar[which(df$patient %in% names(clustered@clusters[[index]])==TRUE)]=names(clustered@exemplars[index])    
    # asign cluster to annotation file
    annot$cluster[which(annot$ID %in% names(clustered@clusters[[index]])==TRUE)]=index
    
  }
  
  # ************************************ 
  # plot heatmap of patients annotated by cluster
  # ************************************ 
  
  ModelsOrderedClustered=Models[order(annot$cluster),]
  annotOrdered=annot[order(annot$cluster),]
  aheatmap(ModelsOrderedClustered, Rowv=NA, Colv=NA,annColors=c("Paired","Dark2","Accent","Dark2","Accent"),distfun="euclidean", annRow=annotOrdered[,c("cluster","Disease.Subtype","condition","Treatment","Center")])
  
  # ************************************ 
  # subset to test plotting small num of clusters
  # ***********************************
  # subset to plot a single cluster
  dfCluster16=df[which(df$patient %in% names(clustered@clusters[[16]])),]
  
  dfCluster16$cluster=factor(dfCluster16$cluster) # turn the cluster number into a factor so that plotting does not treat it as a continuos var
  ggplot(dfCluster16)+geom_point(aes(x=interaction,y=value,colour=patient))
  ggplot(dfCluster16)+geom_boxplot(aes(x=patient,y=value))
  
  # subset to plot three clusters
  dfThreeClusters=df[which(df$patient %in% c(names(clustered@clusters[[10]]),names(clustered@clusters[[11]]),names(clustered@clusters[[12]]))),]
  ThreeExemplars=dfThreeClusters[which(dfThreeClusters$patient %in% unique(dfThreeClusters$exemplar)),]
  
  dfThreeClusters$cluster=factor(dfThreeClusters$cluster) # turn the cluster number into a factor so that plotting does not treat it as a continuos var
  ggplot(dfThreeClusters)+geom_boxplot(aes(x=patient,y=value))+geom_boxplot(data=ThreeExemplars, aes(x=patient,y=value,colour=patient))+theme_bw()+facet_grid(cluster~.)
  
  
  # ************************************ 
  # ggplot the interaction  trajectories colored by cluster to which they belong
  # ************************************ 

  Exemplars=df[which(df$patient %in% unique(df$exemplar)),]
  df$cluster=factor(df$cluster) # turn the cluster number into a factor so that plotting does not treat it as a continuos var
  ggplot(df)+geom_point(aes(x=interaction,y=value))+geom_point(data=Exemplars, aes(x=interaction,y=value,colour=patient))+theme_bw()+facet_grid(cluster~.)
  ggplot(df)+geom_boxplot(aes(x=interaction,y=value))+geom_boxplot(data=Exemplars, aes(x=interaction,y=value,colour=patient))+theme_bw()+facet_grid(cluster~.)
  
  # ************************************ 
  # ggplot the patient intearction valuescolored by cluster to which they belong
  # ************************************ 
  ggplot(df)+geom_point(aes(x=patient,y=value))+geom_point(data=Exemplars, aes(x=patient,y=value,colour=patient))+facet_grid(cluster~.)
  
  
  # ************************************ 
  # calculate patients with most interactions on or off (but not in between)
  # ************************************ 
  
  numPat=dim(Models)[1]
  numInteractions=dim(Models)[2]
  confidentInteractions=vector(length=numPat)
  for (i in 1:numPat){
    
    confidentInteractions[i]=length(which(Models[i,]>.8)) + length(which(Models[i,]<0.2))
    
  }
  
  barplot(confidentInteractions)
  SuperconfidentInteractions=vector(length=numPat)
  for (i in 1:numPat){
    
    SuperconfidentInteractions[i]=length(which(Models[i,]>.9)) + length(which(Models[i,]<0.1))
    
  }
  
  barplot(SuperconfidentInteractions)
  
  MegaconfidentInteractions=vector(length=numPat)
  for (i in 1:numPat){
    
    MegaconfidentInteractions[i]=length(which(Models[i,]>.95)) + length(which(Models[i,]<0.05))
    
  }
  
  barplot(MegaconfidentInteractions)
  
  
  
  # ************************************ 
  # remove uncertain interactions, i.e. those that are for at least one patient between the superconfident treshold
  # reminder: this is a highly stringent filtering
  # ************************************ 
  
  #see clusterConfidentInteractions.R
  
  
  
  
  
  
  
  
  
  
  
  
  


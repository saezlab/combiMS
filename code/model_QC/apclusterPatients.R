# Function to discard bias of final patient models due to technical factors such as
# treatment, center, disease subtipe.

# Final version October 2014, Marti Bernardo-Faura
#
# Minor modifications regarding figure styling and Github repository directories, 2019, Melanie Rinas



#apclusterPatients<-function(Models){
  
  library(apcluster)
  library(gclus)
  library(ggplot2)
  library(reshape2)
  library(NMF)
  library(corrplot)
  
  # define directories / file locations  ------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  
  
  
  # ######################################################
  # get annotation data 
  # # ######################################################
  
  annot=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
  
  
  
  
  #simplify category
  annot2=annot
  annot2$Category[which(annot$Category=="dont know")]='PPMS'
  annot2$Category[which(annot2$Category=='PPMS slow')]='PPMS'
  annot2$Category[which(annot2$Category=='PPMS fast')]='PPMS'
  
  
  rm(annot)
  
  
  
  
  annot3 = annot2
  
  annot3$HealthStatus = annot2$Disease.Subtype
  #annot3$HealthStatus[which(annot2$Group=="Healthy")]='Healthy controls'
  annot3$HealthStatus[which(annot2$Group=="Healthy")]='Healthy'
  
  annot3$HealthStatus[which(annot2$Disease.Subtype=="PP")]='PPMS'
  annot3$HealthStatus[which(annot2$Disease.Subtype=="PR")]=annot2$Category[which(annot2$Disease.Subtype=="PR")]
  annot3$HealthStatus[which(annot2$Disease.Subtype=="RR")]='RRMS'
  
  annot3$HealthStatus[which(annot2$Disease.Subtype=="SP")]=annot2$Category[which(annot2$Disease.Subtype=="SP")]
  #annot3$HealthStatus[which(annot2$Disease.Subtype=="SP")]='SPMS'
  
  
  annot3$Condition = annot2$condition
  
  annot3$Treatment[which(annot2$Treatment=="Copaxone")]='GA'
  annot3$Treatment[which(annot2$Treatment=="Gilenya")]='FTY'
  annot3$Treatment[which(annot2$Treatment=="Tysabri")]='NTZ'
  annot3$Treatment[which(annot2$condition=="Untreated")]='Untreated MS'
  annot3$Treatment[which(annot2$condition=="Untreated")]='Untreated MS'
  annot3$Treatment[which(annot2$condition=="Healthy")]='Healthy'
  annot3$Treatment[which(annot2$condition=="Tysabri or Placebo")]='NTZ or Placebo'
  
  
  annot = annot3
  
  
  
  
  
  
  # *************************************************************************************************************************
  # ***********load predicted networks
  # *************************************************************************************************************************
  
  load(paste("../../files/median_models/allMedianModels.RData",sep=""))  # called allMedianNetworks               
  
  Models=allMedianNetworks

  
  
  
  
  
  # ************************************ 
  # calculate patient correlation as similarity matrix
  # ************************************
  #Cor=cor(t(Models))
  D=dist(Models)
  Cor=as.matrix(D)
  #aheatmap(Cor)
  #corrplot(Cor,order = "hclust", addrect = 3)

  #aheatmap(Cor, annColors=c("Paired","Dark2","Accent"),distfun="euclidean", annRow=annot[,c("Disease.Subtype","condition","Treatment")])
  
  
  
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
  annot$Cluster=vector(length=length(annot$ID))
  
  for (index in 1:length(clustered@clusters)){  
    df$cluster[which(df$patient %in% names(clustered@clusters[[index]])==TRUE)]=index
    df$exemplar[which(df$patient %in% names(clustered@clusters[[index]])==TRUE)]=names(clustered@exemplars[index])    
    # asign cluster to annotation file
    annot$Cluster[which(annot$ID %in% names(clustered@clusters[[index]])==TRUE)]=index
    
  }
  
  # ************************************ 
  # plot heatmap of patients annotated by cluster
  # ************************************ 
  
  ModelsOrderedClustered=Models[order(annot$Cluster),]
  annotOrdered=annot[order(annot$Cluster),]
  
  
  aheatmap(ModelsOrderedClustered, Rowv=NA, Colv=NA,annColors=c("Paired","Dark2","Accent","Dark2","Accent"),distfun="euclidean", annRow=annotOrdered[,c("Cluster","Disease.Subtype","condition","Treatment","Center")])
  
  
  
  # Define colores
  
  color_vec_clusters = c("#D53E4F","#F46D43","#FDAE61","#FEE08B" ,"#FFFFBF", "#E6F598", "#ABDDA4" ,"#66C2A5", "#5E4FA2")
  
  color_annotation = list(color_vec_clusters,                                                   # color for clusters
                          c('#E69679','#407FB7','#A8859E','#89CCCF'),                                           # color for disease subtype
                          c('#407FB7', '#8DC060', '#FABE50'),                        # color for condition
                          c("#0098A1","grey40","#57AB27","#BDCD00","#407FB7","#612158",'lightcyan2',"#7A6FAC",'grey80','#FABE50'),                                                  # color for treatment
                          c('#2D7F83', '#B65256', '#D0D95C', '#9B91C1'))                                                              # color for center
  
  
  
  
  pdf(file.path("../../figures/supplementary_merging/",paste("aheatmap_allMedianNetworks_ModelsCluster_QC.pdf",sep = "")))
  
  
  aheatmap(ModelsOrderedClustered, Rowv=NA, Colv=NA,
           color = c('grey90','black'),
           annColors = color_annotation,distfun="euclidean", annRow=annotOrdered[,c("Cluster","HealthStatus","Condition","Treatment","Center")])
  
  
  dev.off()
  
  
  
  
  
  
  
  
  
  # ************************************ 
  # subset to test plotting small num of clusters
  # ***********************************
  # subset to plot a single cluster
  # 
  test_cluster_ID = 3
  dfClusterTest=df[which(df$patient %in% names(clustered@clusters[[test_cluster_ID]])),]
  
  dfClusterTest$cluster=factor(dfClusterTest$cluster) # turn the cluster number into a factor so that plotting does not treat it as a continuos var
  ggplot(dfClusterTest)+geom_point(aes(x=interaction,y=value,colour=patient))
  ggplot(dfClusterTest)+geom_boxplot(aes(x=patient,y=value))
  
  # subset to plot three clusters
  dfThreeClusters=df[which(df$patient %in% c(names(clustered@clusters[[1]]),names(clustered@clusters[[2]]),names(clustered@clusters[[3]]))),]
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
  
  
  
  
  
  
  
  
  
  
  
  
  




# *************************
# this script loads a similarity matrix which was calculated for each patient in the cluster using calculateSimilarityBetweenMedianCluster.R
# the similarities were merged into one single structure using mergeSimilarityMatrix.R
# Marti July 2015
# *************************
#apclusterPatients<-function(Models){
  
  library(apcluster)
  library(gclus)
  library(ggplot2)
  library(reshape2)
  library(NMF)
  library(corrplot)
  library(gdata)
  
  
  
  
  
  # ************************************ 
  # load annotation file and similarity BETWEEN matrix
  # ************************************
  load("/Users/marti/Documents/R/combiMS/cluster/analysis/similarityMatrix.RData")
  annot=read.csv("/Users/marti/Documents/R/combiMS/modelling/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
#   aheatmap(similarityMatrix, annColors=c("Paired","Dark2","Accent","Paired"), Row=NA, Colv=NA, cexCol=2.5, labRow=NA, annCol=annot[,c("Disease.Subtype","condition","Treatment","Center")])
  # add some fields to annotation for ggplot  
  annot2=annot
  # correct that PPMS patients appear in group instead of by untreated and treatment
  annot2[which(annot2$Group=='PPMS' & annot2$Treatment=='no'),'Group']='Untreated'
  annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Group']=annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Treatment']
  annot2[which(annot2$Group=='Gilenya'),'Group']='Fingolimod'
  annot=annot2
  # ************************************ 
  # include (or not) similarities WITHIN
  # ************************************ 
  load('/Users/marti/Documents/R/combiMS/cluster/analysis/similarityWITHIN.RData')
  #colnames(similarityMatrix)==names(similarityVector)
  for (i in 1:dim(similarityMatrix)[1]){
    similarityMatrix[i,i]=similarityVector[i]
  }
 
  # ************************************ 
  # create dataframe of similarities to plot
  # ************************************
  #similarityMatrix[lower.tri(similarityMatrix,diag=F)]=NA
  sim_df=melt(similarityMatrix)
  sim_df$Group=rep(F,dim(similarityMatrix)[1])
  sim_df$GroupSim=rep(F,dim(similarityMatrix)[1])
  sim_df$Condition=rep(F,dim(similarityMatrix)[1])
  sim_df$ConditionSim=rep(F,dim(similarityMatrix)[1])
  sim_df$Self=rep(F,dim(similarityMatrix)[1])
  sim_df$RR_Untreated=rep(F,dim(similarityMatrix)[1])
  #*** add meta-fields for visualization
  
  #colnames(sim_df)=c('SimWith','Patient','Similarity','Condition','ConditionSim','Self')
  

  colnames(sim_df)=c('SimWith','Patient','Similarity','Group','GroupSim','Condition','ConditionSim','Self','RR_Untreated')  
  sim_df$GroupSim=rep(annot$Group,length(annot$ID))
  sim_df$ConditionSim=rep(annot$condition,length(annot$ID))
  sim_df$all="All"
  
  
  for (i in 1:length(annot$ID)){
    sim_df$Group[sim_df$Patient==annot[i,'ID']] = annot[i,'Group']
    sim_df$Condition[sim_df$Patient==annot[i,'ID']] = annot[i,'condition'] 
  }

  for (i in 1:length(sim_df$Patient)){
    if(as.character(sim_df$Patient[i])==as.character(sim_df$SimWith[i])){
      sim_df$Self[i]='Self'
    } else{sim_df$Self[i]='Another'}
    
    if(as.character(sim_df$Patient[i]) %in% annot[which(annot$Disease.Subtype=='RR' & annot$Treatment=='no'),'ID']){
      sim_df$RR_Untreated[i]='RR_Untreated'
      } else{sim_df$RR_Untreated[i]='N'}
    }
  
#   # ************************************ 
#   # Alternative dataframe based only on upper part of matrix
#   # ************************************
#   similarityMatrix[lower.tri(similarityMatrix,diag=F)]=NA
#   num_elements=(169*168)/2
#   sim_df=data.frame(Patient=vector('numeric'),SimWith=vector(mode='numeric'),Similarity=vector(mode='numeric'),Condition=vector(mode='numeric'),ConditionSim=vector(mode='numeric'),Self=vector(mode='numeric'))
#   sim_df$Patient=colnames(similarityMatrix)

  # ************************************ 
  # plot distribution of similarities
  # ************************************
  g=ggplot(sim_df) + geom_density(aes(x=Similarity, fill=Condition),alpha=0.2) +theme_bw()
  g+geom_density(aes(x=Similarity,fill=all),alpha=0.2)
  
  # ************************************ 
  # plot distribution of subsets where similarities are subsetting only between patients with same conditions
  # and against all
  # ************************************
  sub_same_group_similarity=subset(sim_df,sim_df$Group==sim_df$GroupSim)
  sub_self_similarity=subset(sub_same_group_similarity,sub_same_group_similarity$Self=='Self')
  sub_IFN=subset(sub_same_group_similarity,sub_same_group_similarity$RR_Untreated=='RR_Untreated')

  # boxplot by Group
  # In a notched box plot, the notches extend 1.58 * IQR / sqrt(n). This gives a roughly 95 interval for comparing medians
  # So non-overlapping notches strongly suggest that medians are significantly different
  # This is the case for RR untreated vis RR IFNb
  
  g=ggplot(sub_same_group_similarity) + geom_boxplot(notch=T,aes(x=Group, y=Similarity, fill=Group),alpha=0.2)
  g=g+geom_boxplot(notch=T,data=sim_df,aes(x=all,y=Similarity, fill=all),alpha=0.2)
  g=g+geom_boxplot(notch=T,data=sub_self_similarity,aes(x=Self,y=Similarity,fill=Self),alpha=0.2)
  g=g+geom_boxplot(notch=T,data=sub_IFN,aes(x=RR_Untreated,y=Similarity,fill=RR_Untreated),alpha=0.2)
  #g=g+geom_jitter(data=sub_self_similarity,aes(x=Self,y=Similarity,color=Self))
  g+theme_bw()
  
  
  # boxplot and points by Group
  sub_same_group_similarity=subset(sim_df,sim_df$Group==sim_df$GroupSim)
  g=ggplot(sub_same_group_similarity) + geom_boxplot(notch=T,aes(x=Group, y=Similarity, fill=Group),alpha=0.2)
  g=g+geom_point(aes(x=Group, y=Similarity, color=Group),position = position_jitter(width = 0.05))
  
  g=g+geom_boxplot(notch=T,data=sim_df,aes(x=all,y=Similarity, fill=all),alpha=0.2)
  g=g+geom_point(aes(x=all,y=Similarity,color=all),alpha=0.2)
  
  g=g+geom_boxplot(notch=T,data=sub_self_similarity,aes(x=Self,y=Similarity,fill=Self),alpha=0.2)
  g=g+geom_point(aes(x=Self,y=Similarity,color=Self))
  
  g=g+geom_boxplot(notch=T,data=sub_IFN,aes(x=RR_Untreated,y=Similarity,fill=RR_Untreated),alpha=0.2)
  g=g+geom_point(aes(x=RR_Untreated,y=Similarity,color=RR_Untreated))
  #g=g+geom_jitter(data=sub_self_similarity,aes(x=Self,y=Similarity,color=Self))
  g+theme_bw()
  
  #   # boxplot by condition
  #   sub_same_group_similarity=subset(sim_df,sim_df$Condition==sim_df$ConditionSim)
  #   g=ggplot(sub_same_group_similarity) + geom_boxplot(aes(x=Condition, y=Similarity, color=Condition),alpha=0.2)
  #   g=g+geom_boxplot(data=sim_df,aes(x=all,y=Similarity, color=all))
  #   g=g+geom_boxplot(data=sub_self_similarity,aes(x=Self,y=Similarity,color=Self))
  #   g+theme_bw()
  #   
  #   
  #   # violin plot by Group
  #   sub_same_group_similarity=subset(sim_df,sim_df$Group==sim_df$GroupSim)
  #   g=ggplot(sub_same_group_similarity) + geom_violin(aes(x=Group, y=Similarity, fill=Group),alpha=0.2)
  #   g=g+geom_violin(data=sim_df,aes(x=all,y=Similarity, fill=all),alpha=0.2)
  #   g=g+geom_violin(data=sub_self_similarity,aes(x=Self,y=Similarity,fill=Self),alpha=0.2)
  #   #g=g+geom_jitter(data=sub_self_similarity,aes(x=Self,y=Similarity,color=Self))
  #   g+theme_bw()
  
  # boxplot + jitter datapoints
  g=g+geom_jitter(data=sub_same_group_similarity,aes(Condition,Similarity, color=Condition)) 
  g=g+geom_jitter(data=sim_df,aes(all,Similarity,color=all),alpha=0.2)+ theme_bw()
  g
  #distributions
  g=ggplot(sub_Matrix) + geom_density(aes(x=Similarity, fill=Condition),alpha=0.2) +theme_bw()
  g+geom_density(aes(x=Similarity,fill=all),alpha=0.2)
  
  # ************************************ 
  # plot distribution of subsets where similarities are only between patients with same conditions
  # VS similarities between patients in different groups
  # ************************************
  sub_diff_group_similarity=subset(sim_df,sim_df$Condition!=sim_df$ConditionSim)
  sub_diff_group_similarity$Comparing=paste0(sub_diff_group_similarity$Condition,"_",sub_diff_group_similarity$ConditionSim)
  # boxplot
  g=ggplot(sub_diff_group_similarity) + geom_boxplot(aes(x=Comparing, y=Similarity, color=Comparing),alpha=0.2)
  g=g+geom_boxplot(data=sim_df,aes(x=all,y=Similarity))+theme_bw()
  g
  # histogram
  g=ggplot(sub_diff_group_similarity) + geom_histogram(aes(x=Similarity,fill=Comparing)) + facet_grid(Comparing~.)
  g
  
  # ************************************ 
  # load merged models
  # ************************************
  load(file=paste0("/Users/marti/Documents/R/combiMS/cluster/analysis/fivePerHundredThousTol/","allMedianModels.RData"))
  aheatmap(allMedianNetworks, annColors=c("Paired","Dark2","Accent","Paired"), Row=NA, Colv=NA, labRow=NA, annCol=annot[,c("Disease.Subtype","condition","Treatment","Center")])

  # ************************************ 
  # ap cluster using precalculated similiartiy matrix
  # ************************************ 
  clustered=apcluster(similarityMatrix, p=0.1)
  #cat("give e.g. preference=-30 for a very reduced number of clusters \n")
  cat("affinity propogation optimal number of clusters:", length(clustered@clusters), "\n")
  
  # ************************************ 
  # save clusters in both annotation file and dataframe for different plots
  # ************************************ 
  # ********* create a molten data frame for plotting
  df=melt(allMedianNetworks)
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
  
  ModelsOrderedClustered=allMedianNetworks[order(annot$cluster),]
  annotOrdered=annot[order(annot$cluster),]
  aheatmap(ModelsOrderedClustered, Rowv=NA, Colv=NA,annColors=c("Paired","Dark2","Accent","Dark2","Accent"),annRow=annotOrdered[,c("cluster","Disease.Subtype","condition","Treatment","Center")])
  
  
  # ************************************ 
  # test significance difference in similarites
  # ************************************ 
  sub_same_group_similarity=subset(sim_df,sim_df$Group==sim_df$GroupSim)
  sub_self_similarity=subset(sub_same_group_similarity,sub_same_group_similarity$Self=='Self')
  sub_IFN=subset(sub_same_group_similarity,sub_same_group_similarity$RR_Untreated=='RR_Untreated')
  sub_IFN[which(sub_IFN$RR_Untreated=='RR_Untreated'),'Patient']
  
  Yes=annot[which(annot$Group=='Interferon B'),'ID']
  No=annot[which(annot$Disease.Subtype=='RR' & annot$Treatment=='no'),'ID']
  
  # visualize distributions being tested
  hist(sim_df[which(sim_df$Patient %in% Yes),'Similarity'])
  hist(sim_df[which(sim_df$Patient %in% No),'Similarity'])
  prova=t.test(sim_df[which(sim_df$Patient %in% Yes),'Similarity'],sim_df[which(sim_df$Patient %in% No),'Similarity'])
  
  
  
  
  
  
  
  
  
  
  


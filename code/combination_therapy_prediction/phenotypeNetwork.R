# Funtion called by pathDrugTargetsFinalv3.R as well as calculateDefective.R to calculate the group network for a given set of patients model which belong to the same group
# E.g. patients treated with same drug or same disease subtype
# Final version April 2016
# Marti Bernardo-Faura

#this is the modification of phenotypeNetwork old. here i do not take the links above activity threshold
#this modification allows to calculate defective and only then remove those below threshold
#for the same reason, I cannot turn network into graph, this is now done in pathDrugTargetv3
#
#
#

# applylinkActivityQuantileThreshold ='no' can be used to produce a non-weighted network for plotting and graph search
# The option applylinkActivityQuantileThreshold is only used for grouping mode mean, as in the median model it is not necessary
# 
# applylinkActivityQuantileThreshold ='yes' leads to rounded mean values based on the choosen linkActivityQuantileThreshold

# Changes by Melanie Rinas, August 2018




# Copyright information ================================================================================================================================== -->
   
#  Copyright (c) 2018 - European Molecular Biology Laboratory, European Bioinformatics Institute, UK,
#                       Joint Research Center for Computational Biomedicine (JRC-COMBINE), RWTH-Aachen University, Faculty of Medicine, Aachen, Germany 
# 
#  File author(s): Marti Bernardo-Faura (marti.bernardo.faura@gmail.com), Melanie Rinas (melrinas@gmail.com) 
# 
#  Distributed under the GPLv3 License. 
#  See accompanying file LICENSE.txt or copy at
#  http://www.gnu.org/licenses/gpl-3.0.html 
   
#  ======================================================================================================================================================== -->
   

phenotypeNetwork<-function(IndexPatients, 
                           allMedianNetworks, 
                           model, 
                           mode='mean', 
                           linkActivityQuantileThreshold=0.75,
                           applylinkActivityQuantileThreshold ='no'){

  # *************************************************************************************************************************
  # ***********calculate grand median OR mean of median models for each phenotype
  # *************************************************************************************************************************
  if(mode=='median'){
  phenotypeNw=apply(allMedianNetworks[IndexPatients,],2,median)
  #Mistake: by convention, the median of an even number of multiple 1s is 0.5
  #I forcely correct this below
  phenotypeNw[which(phenotypeNw==0.5)]=1
  }else if (mode=='mean'){
    phenotypeNw=apply(allMedianNetworks[IndexPatients,],2,mean)
    
   
    if(applylinkActivityQuantileThreshold == 'yes'){
       
       
       phenotypeNw[which(phenotypeNw > quantile(phenotypeNw, probs = linkActivityQuantileThreshold))]=1
       phenotypeNw[which(phenotypeNw <= quantile(phenotypeNw, probs = linkActivityQuantileThreshold))]=0
       
    }
    #if mode=mean, select only active part based on the chosen linkActivityQuantileThreshold
   
  } else (stop("Wrong mode. Mean or median required"))


#   # *************************************************************************************************************************
#   # ***********create graph out of active disease network
#   # *************************************************************************************************************************
#   # create model with active part
#   modelPhenotype=cutModel(model,phenotypeNw)
#   plotModel(modelPhenotype,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))
#   # transform model into graph
#   tempSIF=model2sif(model,optimRes=list(bString=phenotypeNw))
#   # replace inhibition by activation, as the graph only accepts positive interactions
#   tempSIF[,2]=abs(as.numeric(tempSIF[,2]))
#   graphModel=sif2graph(tempSIF)
# 
#   ModelGraph=list("model"=modelPhenotype,"graph"=graphModel,"network"=phenotypeNw)
groupNetwork=list("network"=phenotypeNw)  
return(groupNetwork)

}


# *************************************************************************************************************************
# ***********compare mean vs median
# *************************************************************************************************************************
# IdxHealthy = annot$Group == 'Healthy'
# IdxMSuntreated = (annot$Treatment == "no" & annot$Disease.Subtype!='')
# 
# phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
# thisPhenotype=phenotypes[i]
# patientsThisPhenotype= annot$Treatment==thisPhenotype 



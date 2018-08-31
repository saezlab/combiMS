# Function called by pathDrugTargetsFinalv3.R to transform a graph after graph search into a network
# Marti Bernardo-Faura
# Final version April 2016
# 
# 
# 
# 
# 
# Changes by Melanie Rinas, July 2018
# 
# 
# 
# Copyright information ================================================================================================================================== -->

#  Copyright (c) 2018 - European Molecular Biology Laboratory, European Bioinformatics Institute, UK,
#                       Joint Research Center for Computational Biomedicine (JRC-COMBINE), RWTH-Aachen University, Faculty of Medicine, Aachen, Germany 
# 
#  File author(s): Marti Bernardo-Faura (marti.bernardo.faura@gmail.com), Jakob Wirbel, Melanie Rinas (melrinas@gmail.com) 
# 
#  Distributed under the GPLv3 License. 
#  See accompanying file LICENSE.txt or copy at
#  http://www.gnu.org/licenses/gpl-3.0.html 

#  ======================================================================================================================================================== -->

network2graph<-function(phenotypeNw,
                        mode='mean',
                        linkActivityQuantileThreshold=0.75){   # Note: The function network2graph requires a Boolean logic network (using not rounded mean causes several wrong results). 
                                                               # Thus the option applylinkActivityQuantileThreshold must be applylinkActivityQuantileThreshold = 'no' and can not be chosen by the user (to prevent wrong results).     
   

  # *************************************************************************************************************************
  # ***********calculate grand median OR mean of median models for each phenotype
  # *************************************************************************************************************************
  if(mode=='median'){
    
    warning(paste0('median model must be corrected for even number of 1. Hence, no 0.5 should appear\n'))
    
    #NOTHING TO DO. I ASSUME MEDIAN HAS BEEN CORRECTED BY PEHNOTEYPENETWORK.R AS DESCRIBED BELOW
    #phenotypeNw=apply(allMedianNetworks[IndexPatients,],2,median)
  #Mistake: by convention, the median of an even number of multiple 1s is 0.5
  #I forcely correct this below
  #phenotypeNw[which(phenotypeNw==0.5)]=1
  }else if (mode=='mean'){
    #phenotypeNw=apply(allMedianNetworks[IndexPatients,],2,mean)
    #if mode=mean, select only active part (apply 0.5 activation threshold)
    
    phenotypeNw[which(phenotypeNw > quantile(phenotypeNw, probs = linkActivityQuantileThreshold))]=1
    phenotypeNw[which(phenotypeNw <= quantile(phenotypeNw, probs = linkActivityQuantileThreshold))]=0
    
    
    warning(paste0('The results of the function network2graph are generated using the Boolean version of the network (rounded mean values) based on the chosen linkActivityQuantileThreshold.\n'))
    
    
  } else (stop("Wrong mode. Mean or median required"))


  # *************************************************************************************************************************
  # ***********create graph out of active disease network
  # *************************************************************************************************************************
  # create model with active part
  modelPhenotype=cutModel(model,phenotypeNw)
  #plotModel(modelPhenotype,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))
  # transform model into graph
  tempSIF=model2sif(model,optimRes=list(bString=phenotypeNw))
  cutNw=tempSIF
  # replace inhibition by activation, as the graph only accepts positive interactions
  tempSIF[,2]=abs(as.numeric(tempSIF[,2]))
  graphModel=sif2graph(tempSIF)

  ModelGraph=list("model"=modelPhenotype,
                  "graph"=graphModel,
                  "network"=phenotypeNw,
                  "cutNetwork"=cutNw)

return(ModelGraph)

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



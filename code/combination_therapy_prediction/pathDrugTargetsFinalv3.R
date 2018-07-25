#Script for combination therapy prediction
#Final version April 2016

# this is a modification of pathPredictedInts.R to identify paths that connect measurable nodes with stimulus used through interactions
# predicted to be targeted successfully, as described by their score in reverting the disease phenotype to a healthy edge weight
# which was calculated using substractNetworksByEdge.R 

# minor changed for github repository by Jakob Wirbel, 2017
# 
# 
# 
# Changes by Melanie Rinas, July 2018
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



library(CellNOptR)
library(reshape2)
library(NMF)
library(ggplot2)
library(corrplot)
library(beeswarm)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("./calculateDefective.R")
source("./phenotypeNetwork.R")

source("./makeMatrixPathStimSignal.R")
source("./network2graph.R")

figure_folder="../../files/cohort_signaling/"
drugScores_folder='../../files/drugScores/'







# *************************************************************************************************************************
# ***********select algorythm options
# *************************************************************************************************************************


groupingMode='mean'     # 2 groupingMode options are possible: 'mean' or 'median'

linkActivityThreshold=0.5
applyLinkActivityThreshold ='no'

thisDrugable="zero"     # 2 thisDrugable options are possible: "zero" or "negative"
                        # "zero": defectiveInts are interactions with a drugScore equal or smaller than 0 (drugScore <= 0)
                        # "negative": defectiveInts are only interactions with a negative drugScore (drugScore < 0)
                        
searchInactiveInts="yes"






# *************************************************************************************************************************
# ***********load anotation
# *************************************************************************************************************************

# ************load anotation to map patients to groups
data_folder="../../data/phosphos_processed/"
filenames=list.files(data_folder,pattern="*.csv",full.names=FALSE)
annot=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
numPat=length(filenames)
filenames2=filenames
for (j in 1:numPat){
  filenames2[j]=strsplit(filenames[j],"\\.")[[1]][1]
}

# ************load model and midas for annotation
patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
model_path="../../files/model/combiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED.sif"
fileName=patientData[1]
midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)  
sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
#model=preprocessing(midas,model,expansion=FALSE)
numInteractions=length(model$reacID)
sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),numInteractions)

#writeSIF(model,file="/Users/marti/Documents/R/combiMS/preprocessedModelCombi.sif")

# *************************************************************************************************************************
# ***********load predicted networks
# *************************************************************************************************************************

load("../../files/median_models/allMedianModels.RData")


# *************************************************************************************************************************
# *********** 1. Calculate Healthy, MS and drug phenotipic network and identify defective interactions
# *************************************************************************************************************************

defectiveScores=calculateDefective(thisMode=groupingMode,
                                   drugable=thisDrugable,
                                   linkActivityThreshold_used=linkActivityThreshold,
                                   applyLinkActivityThreshold_used = applyLinkActivityThreshold)
#save plots
# defectiveScores$scorePlot
# ggsave(paste0(figure_folder,"scoreCumulativePlot",groupingMode,thisDrugable,".pdf"))
# defectiveScores$numDefectivePlot
# ggsave(paste0(figure_folder,"numDefectivePlot",groupingMode,thisDrugable,".pdf"))
# defectiveScores$alwaysDefective

H=defectiveScores$healthyNW
MS=defectiveScores$MSuntreatedNw
allDrugNws=defectiveScores$allDrugNws
#plotModel(H$model,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))
#plotModel(MS$model,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))
#plotModel(MS$model,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))

# *************************************************************************************************************************
# *********** 1b. Add the calculation of MS subtypes networks, as they are not needed in calculateDefective
# *************************************************************************************************************************
phenotypeNws_folder='../../files/group_models/'


#simplify category
annot2=annot
annot2$Category[which(annot$Category=="dont know")]='PPMS'
annot2$Category[which(annot2$Category=='PPMS slow')]='PPMS'
annot2$Category[which(annot2$Category=='PPMS fast')]='PPMS'

#test if all 169 are correctly labeled
length(which(annot2$Category=='PPMS')) + length(which(annot2$Category=='RRMS')) + length(which(annot2$Category=='healthy'))
length(which(annot2$Category=='RRMS' & annot2$condition=='Untreated')) + length(which(annot2$Category=='PPMS' & annot2$condition=='Untreated')) + length(which(annot2$condition=='Treated')) + length(which(annot2$condition=='Healthy'))

#create PPMS and RRMS networks



applyLinkActivityThreshold__storing_text = ''

if(groupingMode =='mean' && applyLinkActivityThreshold =='no'){
   
   linkActivityThreshold_used_text = gsub('\\.', '_', linkActivityThreshold)
   applyLinkActivityThreshold__storing_text = paste('linkActivityThreshold_',linkActivityThreshold_used_text,'_NOTroundedRealNumber_',sep="")
   
} else if (groupingMode =='mean' && applyLinkActivityThreshold =='yes'){
   
   linkActivityThreshold_used_text = gsub('\\.', '_', linkActivityThreshold)
   applyLinkActivityThreshold__storing_text = paste('linkActivityThreshold_',linkActivityThreshold_used_text,'_rounded_',sep="")
   
}



thisPhenotype='RRMS'
Idx= which(annot2$Category==thisPhenotype & annot2$condition=='Untreated')
thisNw= phenotypeNetwork(Idx, allMedianNetworks,model,mode=groupingMode,linkActivityThreshold=linkActivityThreshold,applyLinkActivityThreshold = applyLinkActivityThreshold)  

write.table(thisNw$network,file=paste(phenotypeNws_folder,thisPhenotype,"__",applyLinkActivityThreshold__storing_text,groupingMode,"_",thisDrugable,".csv",sep=""),sep=",",row.names=T)

thisPhenotype='PPMS'
Idx= which(annot2$Category==thisPhenotype & annot2$condition=='Untreated')
thisNw= phenotypeNetwork(Idx, allMedianNetworks,model,mode=groupingMode,linkActivityThreshold=linkActivityThreshold,applyLinkActivityThreshold = applyLinkActivityThreshold)  
  
write.table(thisNw$network,paste(phenotypeNws_folder,thisPhenotype,"__",applyLinkActivityThreshold__storing_text,groupingMode,"_",thisDrugable,".csv",sep=""),sep=",",row.names=T)


# *************************************************************************************************************************
# *********** Plot Healthy, MS, and drug phenotypic networks -->now in createNwFigure.R
# *************************************************************************************************************************


drug=1     # 5 #**********select a drug here
#Gilenya=Fingolimod, Glatiramer=Copaxone
phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
numDrugs=5
# #plot(allDrugNws[[drug]]$graph)




Figure_folder_storage_name = paste("Figures_based_on__", phenotypes[drug] ,"__",applyLinkActivityThreshold__storing_text,groupingMode,"_",thisDrugable,sep="")        

ifelse(!dir.exists(file.path(figure_folder,Figure_folder_storage_name)), dir.create(file.path(figure_folder,Figure_folder_storage_name)), FALSE)              
figure_folder_of_current_script_settings = file.path(figure_folder,Figure_folder_storage_name)                                                                             
figure_folder_of_current_script_settings_ = paste(figure_folder_of_current_script_settings,"/",sep="")




# *************************************************************************************************************************
# *********** 2. Select only active parts of nw above activity threshold, then turn into graph
# *************************************************************************************************************************
# Transform network into graph. Booleanize network if mode is mean
drugNetwork=network2graph(allDrugNws[[drug]]$network,
                          mode=groupingMode,
                          linkActivityThreshold=linkActivityThreshold)  # Note: The function network2graph requires a Boolean logic network (using not rounded mean causes several wrong results). 
                                                                        # Thus the option applyLinkActivityThreshold must be applyLinkActivityThreshold = 'no' and can not be chosen by the user (to prevent wrong results).     


plotModel(drugNetwork$model,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))


# Save figure in figure_folder_of_current_script_settings_
FigureWidth_plotModel = 22    
FigureHeight_plotModel = 10   

pdf(file=file.path(figure_folder_of_current_script_settings_,paste("plotModel_drugNetwork_model__", phenotypes[drug] ,"__",applyLinkActivityThreshold__storing_text,groupingMode,"_",thisDrugable,".pdf",sep="")),width=FigureWidth_plotModel, height=FigureHeight_plotModel)

plotModel(drugNetwork$model,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))

dev.off()    


#Note the difference between activity in network (grouping), co-drugability score, and "situation not ok"score, i.e. difference between grouping in H and D.
#calculateDefective can work in "zero" scoring mode, identifying possibly inactive interactions with a non-positive co-drugability (groupingMode = median: score of 0 or -1; groupingMode = mean: score <= 0) where situation is not ok. 
# network2graph will cut such interactions out, thereby if we search that graph such interactions cannot be used. A rationale has been
# discussed where it is interesting to use them. Therefore, before prediction of combination therapy, force activation of these interactions.
# note:interactions with 0 scores and situation ok have been replaced by 1 in their score
if(searchInactiveInts=="no"){
   
  warning(paste0("searching inactive ints: ",searchInactiveInts))
   
} else if (searchInactiveInts=="yes"){
   
  phenotypeScores=read.table(file=paste(drugScores_folder,phenotypes[drug],"__",applyLinkActivityThreshold__storing_text,groupingMode,"_",thisDrugable,".csv",sep=""),sep=",")
  
   
  interactions0Score=drugNetwork$network[which(names(drugNetwork$network) %in% rownames(phenotypeScores)[which(phenotypeScores==0)])]
  interactions0ScoreToReplace=names(which(interactions0Score==0))
  
  interactionsNegativeScore=drugNetwork$network[which(names(drugNetwork$network) %in% rownames(phenotypeScores)[which(phenotypeScores<0)])]
  interactionsNegativeScoreToReplace=names(which(interactionsNegativeScore==0))
  
  
  network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts  =  drugNetwork$network   # create new list element to store the modified activity value
  
  network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts[which(names(drugNetwork$network) %in% interactions0ScoreToReplace)]=1 #replace activity by 1
  network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts[which(names(drugNetwork$network) %in% interactionsNegativeScoreToReplace)]=1 #replace activity by 1
  
  for(ii in 1:length(interactions0ScoreToReplace)){
     
     print(warning(paste0("(" , ii , ")  " , interactions0ScoreToReplace[ii]," defective inactive (OFF) interaction replaced with activity=1 (ON).\n")))
     
  }
  warning(paste0("Total of ",length(interactions0ScoreToReplace)," inactive (OFF) ints with 0 score and situation not ok \nout of all (ON and OFF) ",
                       length(rownames(phenotypeScores)[which(phenotypeScores==0)]) ," ints with 0 score and situation not ok.\n"))
  
  
  
  
  
  for(iii in 1:length(interactionsNegativeScoreToReplace)){
     
     print(warning(paste0("(" , iii , ")  " , interactionsNegativeScoreToReplace[iii]," defective inactive (OFF) interaction replaced with activity=1 (ON).\n")))

  }
  warning(paste0("Total of ",length(interactionsNegativeScoreToReplace)," inactive (OFF) ints with negative score and situation not ok \nout of all (ON and OFF) ",
                 length(rownames(phenotypeScores)[which(phenotypeScores<0)]) ," ints with negative score and situation not ok.\n"))
  

  
  
  
  # *************************************************************************************************************************
  # *********** Repeat calculations of network2graph.R using network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts and 
  #             store them as new list elements in drugNetwork
  #             
  #             In this way, drugNetwork will by a list containing following list elements:
  #             (1) "model"=modelPhenotype
  #             (2) "graph"=graphModel
  #             (3) "network"=phenotypeNw
  #             (4) "cutNetwork"=cutNw
  #             
  #             (5) "model_OFFdefectiveInts_ModifiedTo_ONdefectiveInts"
  #             (6) "graph_OFFdefectiveInts_ModifiedTo_ONdefectiveInts"
  #             (7) "network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts"
  #             (8) "cutNetwork_OFFdefectiveInts_ModifiedTo_ONdefectiveInts"
  #             
  # *************************************************************************************************************************
  
  
  # create model with active part
  modelPhenotype_OFFdefectiveInts_ModifiedTo_ONdefectiveInts = cutModel(model,network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts)
  drugNetwork$model_OFFdefectiveInts_ModifiedTo_ONdefectiveInts = modelPhenotype_OFFdefectiveInts_ModifiedTo_ONdefectiveInts
  # Note: For the number of interactions of drugNetwork$model_OFFdefectiveInts_ModifiedTo_ONdefectiveInts it holds:
  #       length(drugNetwork$model_OFFdefectiveInts_ModifiedTo_ONdefectiveInts$reacID) = length(drugNetwork$model$reacID) + length(interactions0ScoreToReplace) + length(interactionsNegativeScoreToReplace)
     
  plotModel(drugNetwork$model_OFFdefectiveInts_ModifiedTo_ONdefectiveInts,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))
  
  
  # Save figure in figure_folder_of_current_script_settings_
  FigureWidth_plotModel__OFFdefectiveInts_ModifiedTo_ONdefectiveInts = 32    
  FigureHeight_plotModel__OFFdefectiveInts_ModifiedTo_ONdefectiveInts = 15   
  
  pdf(file=file.path(figure_folder_of_current_script_settings_,paste("plotModel_drugNetwork_model_OFFdefectiveInts_ModifiedTo_ONdefectiveInts__", phenotypes[drug] ,"__",applyLinkActivityThreshold__storing_text,groupingMode,"_",thisDrugable,".pdf",sep="")),width=FigureWidth_plotModel, height=FigureHeight_plotModel)
  
  plotModel(drugNetwork$model_OFFdefectiveInts_ModifiedTo_ONdefectiveInts,midas,graphvizParams=list(fontsize=35,nodeWidth=4,nodeHeight=4))
  
  dev.off()    
  
  
  # transform model into graph
  tempSIF_OFFdefectiveInts_ModifiedTo_ONdefectiveInts = model2sif(model,optimRes=list(bString=network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts))
  cutNw__OFFdefectiveInts_ModifiedTo_ONdefectiveInts=tempSIF_OFFdefectiveInts_ModifiedTo_ONdefectiveInts
  # replace inhibition by activation, as the graph only accepts positive interactions
  tempSIF_OFFdefectiveInts_ModifiedTo_ONdefectiveInts[,2]=abs(as.numeric(tempSIF_OFFdefectiveInts_ModifiedTo_ONdefectiveInts[,2]))
  graphModel_OFFdefectiveInts_ModifiedTo_ONdefectiveInts=sif2graph(tempSIF_OFFdefectiveInts_ModifiedTo_ONdefectiveInts)
  

  drugNetwork$graph_OFFdefectiveInts_ModifiedTo_ONdefectiveInts = graphModel_OFFdefectiveInts_ModifiedTo_ONdefectiveInts
  drugNetwork$network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts = network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts
  drugNetwork$cutNetwork_OFFdefectiveInts_ModifiedTo_ONdefectiveInts = cutNw__OFFdefectiveInts_ModifiedTo_ONdefectiveInts
  
     
     #  ModelGraph=list("model"=modelPhenotype,"graph"=graphModel,"network"=phenotypeNw,"cutNetwork"=cutNw)

  # transform network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts into graph
  
  tempSIF__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts = model2sif(model,optimRes=list(bString=network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts))
  cut__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts = tempSIF__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts
  
  graphModel__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts=sif2graph(tempSIF__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts)                                                            
  
  # replace inhibition by activation, as the graph only accepts positive interactions
  tempSIF__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts[,2]=abs(as.numeric(tempSIF__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts[,2]))                                                
  graphModel__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts__inhibition_replaced_by_activation=sif2graph(tempSIF__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts)  
  
  
  drugNetwork$graph__codruggableInactiveInts__modified_to_active = graphModel__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts
  drugNetwork$graph__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation = graphModel__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts__inhibition_replaced_by_activation
  drugNetwork$cutnetwork_OFFdefectiveInts_ModifiedTo_ONdefectiveInts = cut__network_OFFdefectiveInts_ModifiedTo_ONdefectiveInts
  

  
} else {stop(warning("incorrect searching mode selected"))}






# *************************************************************************************************************************
# ***********3. Once a drug is selected, recreate drug network, then
#************determine if a stimulus connects to signals through Defective interactions
# *************************************************************************************************************************

# ***********3.1 Determine combinationExperiments using only active (ON) defective interactions # ****************************************************

Stimuli=colnames(midas@stimuli)
AllReadouts=colnames(midas@signals[[2]])
cat('**************Chosen drug:',phenotypes[drug],'\n')
topTargets=names(defectiveScores$defectiveInts[[drug]]) # these are defective interactions

# this allowed us to select only the interactions with a mean difference bigger than 0.15
# to select all interactions ranked by mean difference in targetting, use e.g. NatalizumabScore.csv
#AllInteractions=rownames(predictedInts)
presentTopInts=topTargets[which(topTargets %in% drugNetwork$model$reacID)]
presentStimuli=colnames(midas@stimuli)[which(colnames(midas@stimuli) %in% drugNetwork$graph@nodes)]
presentSignals=colnames(midas@signals[[2]])[which(colnames(midas@signals[[2]]) %in% drugNetwork$graph@nodes)]

cat("Defective interactions:",topTargets,"\n")
cat("Present (active) defective:",presentTopInts,"\n")
cat("Stimuli with activity:",presentStimuli,"\n")
cat("Readouts with activity",presentSignals,"\n")


#combinationExperiments=data.frame(Reaction=character(),Stimulus=numeric(),Readout=numeric(),drugScoreInteraction=numeric())
combinationExperiments=matrix(ncol=4,nrow=0)
colnames(combinationExperiments)=c('Reaction','Stimulus','Readout','drugScore')


for(indexInts in 1:length(presentTopInts)){
   
link=presentTopInts[indexInts]
matrixPaths=makeMatrixPathStimSignal(drugNetwork$graph,link)
numCombinations=length(AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]])
Stimuli[which(matrixPaths==1,arr.ind=T)[,2]]

experimentsThisInteraction=matrix(nrow=numCombinations,ncol=4)
experimentsThisInteraction[,1]=rep(link,numCombinations)
experimentsThisInteraction[,2]=Stimuli[which(matrixPaths==1,arr.ind=T)[,2]]
experimentsThisInteraction[,3]=AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]]
experimentsThisInteraction[,4]=rep(defectiveScores$defectiveInts[[drug]][link],numCombinations) 



# experimentsThisInteraction=data.frame(Reaction=rep(link,numCombinations),
#                                       Stimulus=Stimuli[which(matrixPaths==1,arr.ind=T)[,2]],
#                                       Readout=AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]],
#                                       drugScoreInteraction=rep(defectiveScores$defectiveInts[[drug]][indexInts],numCombinations) 
#                             )
combinationExperiments=rbind(combinationExperiments,experimentsThisInteraction)
#aheatmap(matrixPaths,Rowv=NA,Colv=NA)
}

drugName=phenotypes[drug]


#write.csv(combinationExperiments,file=paste("/Users/marti/Desktop/figuresCombiMS/combinations/",drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))
write.csv(combinationExperiments,file=paste(drugScores_folder,"combinationExperiments_",drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))
write.csv(combinationExperiments,file=paste(drugScores_folder,"combinationExperiments_",drugName,thisDrugable,searchInactiveInts,applyLinkActivityThreshold__storing_text,groupingMode,".csv",sep=""))





# *********** If searchInactiveInts = "yes" ********************************************************************************************************

# ***********3.2 Determine combinationExperiments using active (ON) and inactive (OFF) defective interactions # ****************************************************
#                by repeating 3.2 using 'OFFdefectiveInts_ModifiedTo_ONdefectiveInts' list elements of drugNetwork 


if(searchInactiveInts == "yes"){
   
      
   cat('**************Chosen drug:',phenotypes[drug],'\n')
   
   
   ONandOFF_TopInts=topTargets[which(topTargets %in% drugNetwork$model_OFFdefectiveInts_ModifiedTo_ONdefectiveInts$reacID)]
   ONandOFF_Stimuli=colnames(midas@stimuli)[which(colnames(midas@stimuli) %in% drugNetwork$graph_OFFdefectiveInts_ModifiedTo_ONdefectiveInts@nodes)]
   ONandOFF_Signals=colnames(midas@signals[[2]])[which(colnames(midas@signals[[2]]) %in% drugNetwork$graph_OFFdefectiveInts_ModifiedTo_ONdefectiveInts@nodes)]
   
   cat("Defective interactions:",topTargets,"\n")
   cat("ON and OFF (active and inactive) defective of model_OFFdefectiveInts_ModifiedTo_ONdefectiveInts :",ONandOFF_TopInts,"\n")
   cat("Stimuli of graph_OFFdefectiveInts_ModifiedTo_ONdefectiveInts:",ONandOFF_Stimuli,"\n")
   cat("Readouts of graph_OFFdefectiveInts_ModifiedTo_ONdefectiveInts",ONandOFF_Signals,"\n")
   
   
   #combinationExperiments__ONandOFF_DefectiveInts=data.frame(Reaction=character(),Stimulus=numeric(),Readout=numeric(),drugScoreInteraction=numeric())
   combinationExperiments__ONandOFF_DefectiveInts=matrix(ncol=6,nrow=0)
   colnames(combinationExperiments__ONandOFF_DefectiveInts)=c('Reaction','Stimulus','Readout','drugScore','Boolean network activity','Group network activity')
   
   
   for(indexInts in 1:length(ONandOFF_TopInts)){
      
      link=ONandOFF_TopInts[indexInts]
      matrixPaths=makeMatrixPathStimSignal(drugNetwork$graph_OFFdefectiveInts_ModifiedTo_ONdefectiveInts,link)
      numCombinations=length(AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]])
      Stimuli[which(matrixPaths==1,arr.ind=T)[,2]]
      
      experimentsThisInteraction=matrix(nrow=numCombinations,ncol=6)
      experimentsThisInteraction[,1]=rep(link,numCombinations)
      experimentsThisInteraction[,2]=Stimuli[which(matrixPaths==1,arr.ind=T)[,2]]
      experimentsThisInteraction[,3]=AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]]
      experimentsThisInteraction[,4]=rep(defectiveScores$defectiveInts[[drug]][link],numCombinations) 
      experimentsThisInteraction[,5]=rep(drugNetwork$network[link],numCombinations) 
      experimentsThisInteraction[,6]=rep(defectiveScores$allDrugNws[[drug]]$network[link],numCombinations) 
      
      
      
      # experimentsThisInteraction=data.frame(Reaction=rep(link,numCombinations),
      #                                       Stimulus=Stimuli[which(matrixPaths==1,arr.ind=T)[,2]],
      #                                       Readout=AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]],
      #                                       drugScoreInteraction=rep(defectiveScores$defectiveInts[[drug]][indexInts],numCombinations) 
      #                             )
      combinationExperiments__ONandOFF_DefectiveInts=rbind(combinationExperiments__ONandOFF_DefectiveInts,experimentsThisInteraction)
      #aheatmap(matrixPaths,Rowv=NA,Colv=NA)
}

   drugName=phenotypes[drug]
   
   
   #write.csv(combinationExperiments__ONandOFF_DefectiveInts,file=paste("/Users/marti/Desktop/figuresCombiMS/combinations/",drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))
   write.csv(combinationExperiments__ONandOFF_DefectiveInts,file=paste(drugScores_folder,"combinationExperiments__ONandOFF_DefectiveInts_",drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))
   write.csv(combinationExperiments__ONandOFF_DefectiveInts,file=paste(drugScores_folder,"combinationExperiments__ONandOFF_DefectiveInts_",drugName,thisDrugable,searchInactiveInts,applyLinkActivityThreshold__storing_text,groupingMode,".csv",sep=""))




}





write.table(drugNetwork$cutNetwork,file=paste0(phenotypeNws_folder,drugName,"CUT",thisDrugable,searchInactiveInts,groupingMode,".csv"),sep=",",row.names=T,quote=F)
write.table(drugNetwork$cutNetwork,file=paste0(phenotypeNws_folder,drugName,"CUT",thisDrugable,searchInactiveInts,applyLinkActivityThreshold__storing_text,groupingMode,".csv"),sep=",",row.names=T,quote=F)


write.table(drugNetwork$cutNetwork_OFFdefectiveInts_ModifiedTo_ONdefectiveInts,file=paste0(phenotypeNws_folder,drugName,"CUT_OFFdefectiveInts_ModifiedTo_ONdefectiveInts_",thisDrugable,searchInactiveInts,groupingMode,".csv"),sep=",",row.names=T,quote=F)
write.table(drugNetwork$cutNetwork_OFFdefectiveInts_ModifiedTo_ONdefectiveInts,file=paste0(phenotypeNws_folder,drugName,"CUT_OFFdefectiveInts_ModifiedTo_ONdefectiveInts_",thisDrugable,searchInactiveInts,applyLinkActivityThreshold__storing_text,groupingMode,".csv"),sep=",",row.names=T,quote=F)



save(drugNetwork,file=paste(phenotypeNws_folder,drugName,"_drugNetwork_",thisDrugable,searchInactiveInts,applyLinkActivityThreshold__storing_text,groupingMode,".RData",sep=""))


cat('**************num combinations by stimuli',length(unique(combinationExperiments[,2])),' based on active co-druggable interactions. \n')
cat('**************num combinations by reaction',length(unique(combinationExperiments[,1])),' based on active co-druggable interactions. \n')

if(searchInactiveInts == "yes"){
   
   cat('**************num combinations by stimuli',length(unique(combinationExperiments__ONandOFF_DefectiveInts[,2])),' based on active and inactive co-druggable interactions. \n')
   cat('**************num combinations by reaction',length(unique(combinationExperiments__ONandOFF_DefectiveInts[,1])),' based on active and inactive co-druggable interactions. \n')
}

# *************************************************************************************************************************
# *********** test each group against each other - careful with the interpretation, as it should be one group against all others
# *************************************************************************************************************************
#see pathDrugTargetsFinal.R

print("Script finished.")
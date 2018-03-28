
# pathDrugTargetsFinalv3_MR__using_Cluster_MR_Results_based_on_using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__gaBinaryT1_postCellNOptRupdate__OKpostCellNOptRupdate




# Some notes added
# 
# AND
# 
# paths adapted for CombiMS rerun
# 
# AND
# 
# different sources used
# 
# 
# AND
# 
# some modifications with no influence on the results
# 
# 
# AND
# 
# MR note:
# 
# option searchInactiveInts has no effect as
# the searchInactiveInts="yes" related modification of
# drugNetwork$network
# by setting
# drugNetwork$network[which(names(drugNetwork$network) %in% interactionsToReplace)]=1 #replace activity by 1
# 
# is afterwards not used 
# 
# 
# If inactive interactions should be considerd,
# it is necessary to change / modify 
# drugNetwork$graph
# 
# by replacing 
# (1) co-druggable interactions with score 0 and inactive status in MSdrug subgroup network
# AND (!!!)
# (2) co-druggable interactions with score -1 and inactive status in MSdrug subgroup network
# 
# please see below in code region marked with
# 
# **************************************************************************************************************************************# 
# 
# and please compare the new variables
# 
# drugNetwork$network__codruggableInactiveInts__modified_to_active
# 
# drugNetwork$graph__codruggableInactiveInts__modified_to_active
# drugNetwork$graph__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation
# 
# drugNetwork$cutNetwork__codruggableInactiveInts__modified_to_active
# 
# 
# 
# AND new RESULT!!!!!!!!!!!!!!
# 
# combinationExperiments__Active_AND_Inactive_codruggableInts
# 
# 
# 
# 
# 
# 
# 
# AND
# 
# ERROR correction !!!!!!!!!!!!!!!!!!!!!!!!!!!!
# concerning
# 
# 
# OLD/original:
# #   experimentsThisInteraction[,4]=rep(defectiveScores$defectiveInts[[drugID]][indexInts],numCombinations)   # ERROR!!!!!
#   
#                                                                                                              It holds that
#                                                                                                              indexInts in 1:length(presentTopInts) = number of the present defective interactions of the current drug
#                                                                                                              
#                                                                                                              AND NOT indexInts in 1:length(defectiveScores$defectiveInts[[drugID]]) = number of ALL defective interactions of the current drug
#                                                                                                              
#                                                                                                              with
#                                                                                                              length(presentTopInts) unequal length(defectiveScores$defectiveInts[[drugID]])
#
#
# NEW/corrected:
#   experimentsThisInteraction[,4]=rep(defectiveScores$defectiveInts[[drugID]][link],numCombinations)    # MR ERROR correction



# ERROR correction has an influence on 
# 
#    combinationExperiments=rbind(combinationExperiments,experimentsThisInteraction)
# stored at
#     write.csv(combinationExperiments,file=paste("/Users/marti/Desktop/figuresCombiMS/combinations/",drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))






# AND
# 
# ERROR correction with no influence on the result
# concerning
# 
# storing path adaption
# OLD/original:
#     write.csv(combinationExperiments,file=paste("/Users/marti/Desktop/figuresCombiMS/combinations/",drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))




# 
# by 
# Melanie Rinas
# December 2017 and January 2018




Results_storage_name = "OptCombiMSclusterFork10_MR__using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__gaBinaryT1_postCellNOptRupdate"        # MR inserted

data_used = "phosphos_processed__original_MB_JW"          # MR inserted


#model_used =   "combiMS_rerun_whole_pkn_preprocessed_MR_211117.sif"                      # MR inserted
model_used =   "combiMSplaneCUT.sif"                      # MR inserted





# pathDrugTargetsFinalv3.R

#Script for combination therapy prediction
#Final version April 2016

# this is a modification of pathPredictedInts.R to identify paths that connect measurable nodes with stimulus used through interactions
# predicted to be targeted successfully, as described by their score in reverting the disease phenotype to a healthy edge weight
# which was calculated using substractNetworksByEdge.R 

# minor changed for github repository by Jakob Wirbel, 2017

library(CellNOptR)
library(reshape2)
library(NMF)
library(ggplot2)
library(corrplot)
library(beeswarm)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#system("defaults write org.R-project.R force.LANG en_US.UTF-8")

#source('~/.active-rstudio-document', echo=TRUE)
# source.with.encoding('~/.active-rstudio-document', encoding='UTF-8', echo=TRUE) 





source("./calculateDefective.R")
#source("./calculateDefective_MR.R")     # MR modified

source("./phenotypeNetwork.R")

#source("./network2graph.R")
source("./network2graph_MR.R")         # MR modified

source("./makeMatrixPathStimSignal.R")
#source("./makeMatrixPathStimSignal_MR.R")           # MR modified







# figure_folder="../../files/cohort_signaling/"     # MR modified
# drugScores_folder='../../files/drugScores/'       # MR modified
figure_folder="../../files/cohort_signaling_MR"     # MR modified
drugScores_folder='../../files/drugScores_MR'       # MR modified

phenotypeNws_folder='../../files/group_models_MR'       # MR modified

ifelse(!dir.exists(file.path(figure_folder,Results_storage_name)), dir.create(file.path(figure_folder,Results_storage_name)), FALSE)              # MR inserted
figure_folder_of_current_script = file.path(figure_folder,Results_storage_name)                                                                              # MR inserted
figure_folder_of_current_script_ = paste(figure_folder_of_current_script,"/",sep="")


ifelse(!dir.exists(file.path(drugScores_folder,Results_storage_name)), dir.create(file.path(drugScores_folder,Results_storage_name)), FALSE)              # MR inserted
drugScores_folder_of_current_script = file.path(drugScores_folder,Results_storage_name)                                                                              # MR inserted
drugScores_folder_of_current_script_ = paste(drugScores_folder_of_current_script,"/",sep="")


ifelse(!dir.exists(file.path(phenotypeNws_folder,Results_storage_name)), dir.create(file.path(phenotypeNws_folder,Results_storage_name)), FALSE)              # MR inserted
phenotypeNws_folder_of_current_script = file.path(phenotypeNws_folder,Results_storage_name)                                                                              # MR inserted
phenotypeNws_folder_of_current_script_ = paste(phenotypeNws_folder_of_current_script,"/",sep="")


# *************************************************************************************************************************
# ***********select algorythm options
# *************************************************************************************************************************
groupingMode='median'
thisDrugable="zero"
searchInactiveInts="yes"                     # MR note:
# option searchInactiveInts has no effect as
# the searchInactiveInts="yes" related modification of
# drugNetwork$network
# by setting
# drugNetwork$network[which(names(drugNetwork$network) %in% interactionsToReplace)]=1 #replace activity by 1
# 
# is afterwards not used 
# 
# 
# If inactive interactions should be considerd,
# it is necessary to change / modify 
# drugNetwork$graph
# 
# by replacing 
# (1) co-druggable interactions with score 0 and inactive status in MSdrug subgroup network
# AND (!!!)
# (2) co-druggable interactions with score -1 and inactive status in MSdrug subgroup network
# 
# please see below in code region marked with
# 
# **************************************************************************************************************************************
# 
# and please compare the new variables
# 
# drugNetwork$network__codruggableInactiveInts__modified_to_active
# 
# drugNetwork$graph__codruggableInactiveInts__modified_to_active
# drugNetwork$graph__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation
# 
# drugNetwork$cutNetwork__codruggableInactiveInts__modified_to_active
# 
# 
# 
# AND new RESULT!!!!!!!!!!!!!!
# 
# combinationExperiments__Active_AND_Inactive_codruggableInts











# *************************************************************************************************************************
# ***********load anotation
# *************************************************************************************************************************

# ************load anotation to map patients to groups
#data_folder="../../data/phosphos_processed/"
data_folder=paste("../../data/", data_used , "/",sep="")   # MR modified

filenames=list.files(data_folder,pattern="*.csv",full.names=FALSE)


annot=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
numPat=length(filenames)
filenames2=filenames
for (j in 1:numPat){
  filenames2[j]=strsplit(filenames[j],"\\.")[[1]][1]
}




# ************load model and midas for annotation
patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)


#model_path="../../files/model/combiMSplaneCUT.sif"
model_path=paste("../../files/model/", model_used, sep="")   # MR modified

fileName=patientData[1]
midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)  

#sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))   MR modified
#model=preprocessing(midas,model,expansion=FALSE)
numInteractions=length(model$reacID)
sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),numInteractions)

#writeSIF(model,file="/Users/marti/Documents/R/combiMS/preprocessedModelCombi.sif")

# *************************************************************************************************************************
# ***********load predicted networks
# *************************************************************************************************************************

#load("../../files/median_models/allMedianModels.RData")
load(paste("../../files/median_models_MR/", Results_storage_name ,"/allMedianModels.RData",sep=""))  # called allMedianNetworks               # MR modified







# *************************************************************************************************************************
# *********** 1. Calculate Healthy, MS and drug phenotipic network and identify defective interactions
# *************************************************************************************************************************

defectiveScores=calculateDefective(thisMode=groupingMode,drugable=thisDrugable)
#defectiveScores=calculateDefective_MR__with_Results_storage_name_Input(thisMode=groupingMode,drugable=thisDrugable,Results_storage_name=Results_storage_name)   # MR modified
#defectiveScores=calculateDefective_MR(thisMode=groupingMode,drugable=thisDrugable)   # MR modified



save(defectiveScores,file=paste(drugScores_folder_of_current_script_,"defectiveScores_",thisDrugable,"_",searchInactiveInts,"_",groupingMode,".RData",sep=""))                                 # MR inserted


#save plots
# defectiveScores$scorePlot
# ggsave(paste0(figure_folder_of_current_script_,"scoreCumulativePlot",groupingMode,thisDrugable,".pdf"))
# defectiveScores$numDefectivePlot
# ggsave(paste0(figure_folder_of_current_script_,"numDefectivePlot",groupingMode,thisDrugable,".pdf"))
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

#phenotypeNws_folder='../../files/group_models/'
phenotypeNws_folder='../../files/group_models_MR/'     # MR modified

#simplify category
annot2=annot
annot2$Category[which(annot$Category=="dont know")]='PPMS'
annot2$Category[which(annot2$Category=='PPMS slow')]='PPMS'
annot2$Category[which(annot2$Category=='PPMS fast')]='PPMS'

#test if all 169 are correctly labeled

length(which(annot2$Category=='PPMS')) + length(which(annot2$Category=='RRMS')) + length(which(annot2$Category=='healthy'))

length(which(annot2$Category=='RRMS' & annot2$condition=='Untreated')) +
  length(which(annot2$Category=='PPMS' & annot2$condition=='Untreated')) + 
  length(which(annot2$condition=='Treated')) + 
  length(which(annot2$condition=='Healthy'))




numDrugs=5

num_patients_in_subgroups=vector(length=(2+numDrugs+2))                                                # MR inserted
names(num_patients_in_subgroups)=c("Healthy",
                                   'MS untreated',
                                   'Gilenya','IFNb','Copaxone','EGCG','Tysabri',
                                   'RRMS untreated','PPMS untreated')          # MR inserted


num_patients_in_subgroups[1]=  length(which(annot2$condition=='Healthy'))
num_patients_in_subgroups[2]=  length(which(annot2$Treatment == "no" & annot2$Disease.Subtype!=''))

num_patients_in_subgroups[3]=  length(which(annot2$Treatment == names(num_patients_in_subgroups)[3]))
num_patients_in_subgroups[4]=  length(which(annot2$Treatment == names(num_patients_in_subgroups)[4]))
num_patients_in_subgroups[5]=  length(which(annot2$Treatment == names(num_patients_in_subgroups)[5]))
num_patients_in_subgroups[6]=  length(which(annot2$Treatment == names(num_patients_in_subgroups)[6]))
num_patients_in_subgroups[7]=  length(which(annot2$Treatment == names(num_patients_in_subgroups)[7]))

num_patients_in_subgroups[8]=  length(which(annot2$Category=='RRMS' & annot2$condition=='Untreated'))
num_patients_in_subgroups[9]=  length(which(annot2$Category=='PPMS' & annot2$condition=='Untreated'))

num_patients_considerd_in_subgroups = sum(num_patients_in_subgroups[1:7])

#create PPMS and RRMS networks
thisPhenotype='RRMS'
Idx= which(annot2$Category==thisPhenotype & annot2$condition=='Untreated')
thisNw= phenotypeNetwork(Idx, allMedianNetworks,model,mode=groupingMode)  
write.table(thisNw$network,file=paste0(phenotypeNws_folder_of_current_script_,thisPhenotype,groupingMode,".csv"),sep=",",row.names=T)

thisPhenotype='PPMS'
Idx= which(annot2$Category==thisPhenotype & annot2$condition=='Untreated')
thisNw= phenotypeNetwork(Idx, allMedianNetworks,model,mode=groupingMode)  
write.table(thisNw$network,file=paste0(phenotypeNws_folder_of_current_script_,thisPhenotype,groupingMode,".csv"),sep=",",row.names=T)









# *************************************************************************************************************************
# *********** Plot Healthy, MS, and drug phenotypic networks -->now in createNwFigure.R
# *************************************************************************************************************************


#drug=5 #**********select a drug here                                                           # MR commented as using a loop to go through all drugs



#Gilenya=Fingolimod, Glatiramer=Copaxone
phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
numDrugs=5
# #plot(allDrugNws[[drug]]$graph)


# *************************************************************************************************************************
# *********** 2. Select only active parts of nw above activity threshold, then turn into graph
# *************************************************************************************************************************
# Transform network into graph. Booleanize network if mode is mean




for(drugID in 1:numDrugs){                                                          # MR inserted
  
  
  #drugNetwork=network2graph(allDrugNws[[drug]]$network,mode=groupingMode)             # MR modified
  #drugNetwork=network2graph(allDrugNws[[drug]]$network,mode="mean")    # MR inserted for testing
  #
  drugNetwork=network2graph_MR(allDrugNws[[drugID]]$network,mode=groupingMode)             # MR modified
  
  
  FigureWidth_plotModel = 22    # MR inserted
  FigureHeight_plotModel = 10   # MR inserted
  
  pdf(file=file.path(figure_folder_of_current_script_,paste("plotModel_drugNetwork_",phenotypes[drugID] ,".pdf", sep="")),width=FigureWidth_plotModel, height=FigureHeight_plotModel)    # MR inserted
  
  plotModel(drugNetwork$model,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))
  
  dev.off()    # MR inserted
  
  
  
  
  
  
  
  
  
  #Note the difference between activity in network (grouping), 
  #co-drugability score, and "situation not ok"score, i.e. difference between grouping in H and D.
  #
  #calculateDefective can work in "zero" scoring mode, 
  #identifying some interactions with a socre of 0 where situation is not ok. 
  #
  # network2graph will cut such interactions out, 
  # thereby if we search that graph such interactions cannot be used. 
  # A rationale has been discussed where it is interesting to use them. 
  # Therefore, before prediction of combination therapy, force activation of these interactions.
  
  # note:intereactions with 0 scores and situation ok have been replaced by 1 in their score, 
  # hence it is safe to just activate all remaining ints with 0 scores
  # 
  
  if(searchInactiveInts=="no"){
    
    warning(paste0("searching inactive ints: ",searchInactiveInts))
    
  } else if (searchInactiveInts=="yes"){
    
    
    
    # MR note:
    # 
    # option searchInactiveInts has no effect as
    # the searchInactiveInts="yes" related modification of
    # drugNetwork$network
    # by setting
    # drugNetwork$network[which(names(drugNetwork$network) %in% interactionsToReplace)]=1 #replace activity by 1
    # 
    # is afterwards not used 
    # 
    # 
    # If inactive interactions should be considerd,
    # it is necessary to change / modify 
    # drugNetwork$graph
    # 
    # by replacing 
    # (1) co-druggable interactions with score 0 and inactive status in MSdrug subgroup network
    # AND (!!!)
    # (2) co-druggable interactions with score -1 and inactive status in MSdrug subgroup network
    # 
    # please see below in code region marked with
    # 
    # **************************************************************************************************************************************
    # 
    # 
    # and please compare the new variables
    # 
    # drugNetwork$network__codruggableInactiveInts__modified_to_active
    # 
    # drugNetwork$graph__codruggableInactiveInts__modified_to_active
    # drugNetwork$graph__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation
    # 
    # drugNetwork$cutNetwork__codruggableInactiveInts__modified_to_active
    # 
    # 
    # 
    # AND new RESULT!!!!!!!!!!!!!!
    # 
    # combinationExperiments__Active_AND_Inactive_codruggableInts
    
    
    
    
    phenotypeScores=read.table(file=paste0(drugScores_folder_of_current_script_,phenotypes[drugID],groupingMode,thisDrugable,".csv"),sep=",")
    interactions0Score=drugNetwork$network[which(names(drugNetwork$network) %in% rownames(phenotypeScores)[which(phenotypeScores==0)])]
    interactionsToReplace=names(which(interactions0Score==0))
    
    
    drugNetwork$network__InactiveInts_with_0Score__modified_to_active  =  drugNetwork$network                                                                     # MR inserted
    
    
    
    #drugNetwork$network[which(names(drugNetwork$network) %in% interactionsToReplace)]=1 #replace activity by 1
    drugNetwork$network__InactiveInts_with_0Score__modified_to_active[which(names(drugNetwork$network) %in% interactionsToReplace)]=1 #replace activity by 1      # MR modified
    
    
    warning(paste0(interactionsToReplace," ints replaced with activity=1. Total of ",length(interactions0Score)," ints with 0 score and situation not ok"))
    
    
    
    
    
    # START MR inserted **************************************************************************************************************************************
    
    
    phenotypeScores=read.table(file=paste0(drugScores_folder_of_current_script_,phenotypes[drugID],groupingMode,thisDrugable,".csv"),sep=",")
    interactions0Score=drugNetwork$network[which(names(drugNetwork$network) %in% rownames(phenotypeScores)[which(phenotypeScores==0)])]
    interactions0ScoreToReplace=names(which(interactions0Score==0))
    
    interactionsMinus1Score=drugNetwork$network[which(names(drugNetwork$network) %in% rownames(phenotypeScores)[which(phenotypeScores==-1)])]
    interactionsMinus1ScoreToReplace=names(which(interactionsMinus1Score==0))
    
    
    drugNetwork$network__codruggableInactiveInts__modified_to_active  =  drugNetwork$network                                                                     # MR inserted
    
    
    drugNetwork$network__codruggableInactiveInts__modified_to_active[which(names(drugNetwork$network) %in% interactions0ScoreToReplace)]=1 #replace activity by 1      
    drugNetwork$network__codruggableInactiveInts__modified_to_active[which(names(drugNetwork$network) %in% interactionsMinus1ScoreToReplace)]=1 #replace activity by 1     
    
    
    warning(paste0(interactions0ScoreToReplace," ints replaced with activity=1. Total of ",length(interactions0ScoreToReplace)," ints with 0 score and situation not ok"))
    warning(paste0(interactionsMinus1ScoreToReplace," ints replaced with activity=1. Total of ",length(interactionsMinus1ScoreToReplace)," ints with -1 score and situation not ok"))
    
    
    
    
    # ____________________________________________________________________________________________________________________________________________________________________
    # The following lines are adaptations of some commands used in 
    # network2graph.R
    
    # transform network__codruggableInactiveInts__modified_to_active into graph
    tempSIF__network__codruggableInactiveInts__modified_to_active = model2sif(model,optimRes=list(bString=drugNetwork$network__codruggableInactiveInts__modified_to_active))
    cut__network__codruggableInactiveInts__modified_to_active = tempSIF__network__codruggableInactiveInts__modified_to_active
    
    graphModel__network__codruggableInactiveInts__modified_to_active=sif2graph(tempSIF__network__codruggableInactiveInts__modified_to_active)                                                            
    
    # replace inhibition by activation, as the graph only accepts positive interactions
    tempSIF__network__codruggableInactiveInts__modified_to_active[,2]=abs(as.numeric(tempSIF__network__codruggableInactiveInts__modified_to_active[,2]))                                                
    graphModel__network__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation=sif2graph(tempSIF__network__codruggableInactiveInts__modified_to_active)  
    
    
    drugNetwork$graph__codruggableInactiveInts__modified_to_active = graphModel__network__codruggableInactiveInts__modified_to_active
    drugNetwork$graph__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation = graphModel__network__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation
    drugNetwork$cutNetwork__codruggableInactiveInts__modified_to_active = cut__network__codruggableInactiveInts__modified_to_active
    
    # ____________________________________________________________________________________________________________________________________________________________________
    
    
    # END MR inserted **************************************************************************************************************************************
    
    
    
    
    
    
    
  } else {stop(warning("incorrect searching mode selected"))}
  
  
  
  
  
  
  
  # *************************************************************************************************************************
  # ***********3. Once a drug is selected, recreate drug network, then
  #************determine if a stimulus connects to signals through Defective interactions
  # *************************************************************************************************************************
  Stimuli=colnames(midas@stimuli)
  AllReadouts=colnames(midas@signals[[2]])
  cat('**************Chosen drug:',phenotypes[drugID],'\n')
  topTargets=names(defectiveScores$defectiveInts[[drugID]]) # these are defective interactions
  
  # this allowed us to select only the interactions with a mean difference bigger than 0.15
  # to select all interactions ranked by mean difference in targetting, use e.g. NatalizumabScore.csv
  #AllInteractions=rownames(predictedInts)
  #
  
  # Testing whether 
  # > identical(drugNetwork$graph_original@nodes,drugNetwork$graph_inhibition_replaced_by_activation@nodes)
  # [1] TRUE
  
  presentTopInts=topTargets[which(topTargets %in% drugNetwork$model$reacID)]
  presentStimuli=colnames(midas@stimuli)[which(colnames(midas@stimuli) %in% drugNetwork$graph_original@nodes)]
  presentSignals=colnames(midas@signals[[2]])[which(colnames(midas@signals[[2]]) %in% drugNetwork$graph_original@nodes)]
  
  cat("Defective intearctions:",topTargets,"\n")
  cat("Present defective:",presentTopInts,"\n")
  cat("Stimuli with activity:",presentStimuli,"\n")
  cat("Readouts with activity",presentSignals,"\n")
  
  
  #combinationExperiments=data.frame(Reaction=character(),Stimulus=numeric(),Readout=numeric(),drugScoreInteraction=numeric())
  combinationExperiments=matrix(ncol=4,nrow=0)
  colnames(combinationExperiments)=c('Reaction','Stimulus','Readout','drugScore')
  
  
  
  
  for(indexInts in 1:length(presentTopInts)){
    
    link=presentTopInts[indexInts]
    
    # MR inserted the following comment:
    # 
    # The function called
    # 
    # makeMatrixPathStimSignal
    # is using the function
    # 
    # is.connected             # 
    # and is.connected is using the function 
    # 
    # 'sp.between' 
    # and 'sp.between' requires that all edge weights are nonnegative!!!
    # 
    # Therfore use as input 
    # drugNetwork$graph_inhibition_replaced_by_activation
    # 
    # and NOT
    # drugNetwork$graph_original
    
    # matrixPaths=makeMatrixPathStimSignal_MR(drugNetwork$graph_inhibition_replaced_by_activation,link)       # MR modified
    matrixPaths=makeMatrixPathStimSignal(drugNetwork$graph_inhibition_replaced_by_activation,link)      
    
    numCombinations=length(AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]])
    Stimuli[which(matrixPaths==1,arr.ind=T)[,2]]
    
    experimentsThisInteraction=matrix(nrow=numCombinations,ncol=4)
    experimentsThisInteraction[,1]=rep(link,numCombinations)
    experimentsThisInteraction[,2]=Stimuli[which(matrixPaths==1,arr.ind=T)[,2]]
    experimentsThisInteraction[,3]=AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]]
    #   experimentsThisInteraction[,4]=rep(defectiveScores$defectiveInts[[drugID]][indexInts],numCombinations)   # ERROR!!!!!
    #   
    #                                                                                                              It holds that
    #                                                                                                              indexInts in 1:length(presentTopInts) = number of the present defective interactions of the current drug
    #                                                                                                              
    #                                                                                                              AND NOT indexInts in 1:length(defectiveScores$defectiveInts[[drugID]]) = number of ALL defective interactions of the current drug
    #                                                                                                              
    #                                                                                                              with
    #                                                                                                              length(presentTopInts) unequal length(defectiveScores$defectiveInts[[drugID]])
    
    experimentsThisInteraction[,4]=rep(defectiveScores$defectiveInts[[drugID]][link],numCombinations)    # MR ERROR correction
    
    
    # experimentsThisInteraction=data.frame(Reaction=rep(link,numCombinations),
    #                                       Stimulus=Stimuli[which(matrixPaths==1,arr.ind=T)[,2]],
    #                                       Readout=AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]],
    #                                       drugScoreInteraction=rep(defectiveScores$defectiveInts[[drugID]][indexInts],numCombinations) 
    #                             )
    combinationExperiments=rbind(combinationExperiments,experimentsThisInteraction)
    #aheatmap(matrixPaths,Rowv=NA,Colv=NA)
  }
  
  
  
  
  drugName=phenotypes[drugID]
  #write.csv(combinationExperiments,file=paste("/Users/marti/Desktop/figuresCombiMS/combinations/",drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))
  write.csv(combinationExperiments,file=paste(drugScores_folder_of_current_script_,drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))                                 # MR modified / ERROR correction
  write.csv(combinationExperiments,file=paste(drugScores_folder_of_current_script_,"combinationExperiments_",drugName,"_",thisDrugable,"_",searchInactiveInts,"_",groupingMode,".csv",sep=""))                                 # MR modified / ERROR correction
  
  write.table(drugNetwork$cutNetwork,file=paste0(phenotypeNws_folder_of_current_script_,drugName,"CUT",thisDrugable,searchInactiveInts,groupingMode,".csv"),sep=",",row.names=T,quote=F)
  write.table(drugNetwork$cutNetwork,file=paste0(phenotypeNws_folder_of_current_script_,"cutNetwork_",drugName,"_","CUT","_",thisDrugable,"_",searchInactiveInts,"_",groupingMode,".csv"),sep=",",row.names=T,quote=F)
  
  
  cat('**************num combinations by stimuli',length(unique(combinationExperiments[,2])),'\n')
  cat('**************num combinations by reaction',length(unique(combinationExperiments[,1])),'\n')
  # *************************************************************************************************************************
  # *********** test each group against each other - careful with the interpretation, as it should be one group against all others
  # *************************************************************************************************************************
  #see pathDrugTargetsFinal.R
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # START MR inserted 
  
  
  
  # *************************************************************************************************************************
  # *********** combinationExperiments__Active_AND_Inactive_codruggableInts
  # *************************************************************************************************************************
  
  
  combinationExperiments__Active_AND_Inactive_codruggableInts=matrix(ncol=5,nrow=0)
  colnames(combinationExperiments__Active_AND_Inactive_codruggableInts)=c('Reaction','Stimulus','Readout','drugScore','Active_Inactive_Status')
  
  
  
  
  for(indexAllInts in 1:length(topTargets)){
    
    link=topTargets[indexAllInts]
    
    # MR inserted the following comment:
    # 
    # The function called
    # 
    # makeMatrixPathStimSignal
    # is using the function
    # 
    # is.connected             # 
    # and is.connected is using the function 
    # 
    # 'sp.between' 
    # and 'sp.between' requires that all edge weights are nonnegative!!!
    # 
    # Therfore use as input 
    # drugNetwork$graph__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation
    # 
    # and NOT
    # drugNetwork$graph__codruggableInactiveInts__modified_to_active
    
    #matrixPaths_ALL=makeMatrixPathStimSignal_MR(drugNetwork$graph__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation,link)       # MR modified
    matrixPaths_ALL=makeMatrixPathStimSignal(drugNetwork$graph__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation,link)       # MR modified
    
    numCombinations=length(AllReadouts[which(matrixPaths_ALL==1,arr.ind=T)[,1]])
    Stimuli[which(matrixPaths_ALL==1,arr.ind=T)[,2]]
    
    experimentsThisInteraction=matrix(nrow=numCombinations,ncol=5)
    experimentsThisInteraction[,1]=rep(link,numCombinations)
    experimentsThisInteraction[,2]=Stimuli[which(matrixPaths_ALL==1,arr.ind=T)[,2]]
    experimentsThisInteraction[,3]=AllReadouts[which(matrixPaths_ALL==1,arr.ind=T)[,1]]
    
    experimentsThisInteraction[,4]=rep(defectiveScores$defectiveInts[[drugID]][link],numCombinations)    # MR ERROR correction
    
    experimentsThisInteraction[,5] = drugNetwork$network[link]
    
    
    
    
    # experimentsThisInteraction=data.frame(Reaction=rep(link,numCombinations),
    #                                       Stimulus=Stimuli[which(matrixPaths_ALL==1,arr.ind=T)[,2]],
    #                                       Readout=AllReadouts[which(matrixPaths_ALL==1,arr.ind=T)[,1]],
    #                                       drugScoreInteraction=rep(defectiveScores$defectiveInts[[drugID]][indexAllInts],numCombinations) 
    #                             )
    combinationExperiments__Active_AND_Inactive_codruggableInts=rbind(combinationExperiments__Active_AND_Inactive_codruggableInts,experimentsThisInteraction)
    #aheatmap(matrixPaths_ALL,Rowv=NA,Colv=NA)
  }
  
  
  
  
  drugName=phenotypes[drugID]
  #write.csv(combinationExperiments__Active_AND_Inactive_codruggableInts,file=paste("/Users/marti/Desktop/figuresCombiMS/combinations/",drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))
  #write.csv(combinationExperiments__Active_AND_Inactive_codruggableInts,file=paste(drugScores_folder_of_current_script_,drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))                                 # MR modified / ERROR correction
  write.csv(combinationExperiments__Active_AND_Inactive_codruggableInts,file=paste(drugScores_folder_of_current_script_,"combinationExperiments__Active_AND_Inactive_codruggableInts_",drugName,"_",thisDrugable,"_",searchInactiveInts,"OK_",groupingMode,".csv",sep=""))                                 # MR modified / ERROR correction
  
  #write.table(drugNetwork$cutNetwork,file=paste0(phenotypeNws_folder_of_current_script_,drugName,"CUT",thisDrugable,searchInactiveInts,groupingMode,".csv"),sep=",",row.names=T,quote=F)
  write.table(drugNetwork$network__codruggableInactiveInts__modified_to_active,file=paste0(phenotypeNws_folder_of_current_script_,"cutNetwork__codruggableInactiveInts__modified_to_active_",drugName,"_","CUT","_",thisDrugable,"_",searchInactiveInts,"_",groupingMode,".csv"),sep=",",row.names=T,quote=F)
  
  
  cat('**************num combinationExperiments__Active_AND_Inactive_codruggableInts by stimuli',length(unique(combinationExperiments__Active_AND_Inactive_codruggableInts[,2])),'\n')
  cat('**************num combinationExperiments__Active_AND_Inactive_codruggableInts by reaction',length(unique(combinationExperiments__Active_AND_Inactive_codruggableInts[,1])),'\n')
  
  
  
  
  
  
  # *************************************************************************************************************************
  # *********** test each group against each other - careful with the interpretation, as it should be one group against all others
  # *************************************************************************************************************************
  #see pathDrugTargetsFinal.R
  
  
  
  
  
  
  # END MR inserted 
  
  
  
  
  
  
  
  
  
  
  
} # MR inserted end of for loop drugID









print("Script finished!")              # MR inserted

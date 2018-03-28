
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






#source("../utils/newJaccard_MR.R")     # MR inserted
source("../utils/newJaccard.R")     # MR inserted





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
    
    #matrixPaths=makeMatrixPathStimSignal_MR(drugNetwork$graph_inhibition_replaced_by_activation,link)       # MR modified
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
    matrixPaths_ALL=makeMatrixPathStimSignal(drugNetwork$graph__codruggableInactiveInts__modified_to_active__inhibition_replaced_by_activation,link)       
    
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







# START MR inserted 


# *************************************************************************************************************************
# *********** Some analysis and plotting of the results
# *************************************************************************************************************************






names(num_patients_in_subgroups)=c("Healthy",
                                   'MS untreated',
                                   'Gilenya','IFNb','Copaxone','EGCG','Tysabri',
                                   'RRMS untreated','PPMS untreated')          # MR inserted


num_patients_in_drug_subgroups = num_patients_in_subgroups[3:7]

num_codruggableInteractions_of_drug_subgroups = c(length(defectiveScores$defectiveInts[[1]]),
                                                  length(defectiveScores$defectiveInts[[2]]),
                                                  length(defectiveScores$defectiveInts[[3]]),
                                                  length(defectiveScores$defectiveInts[[4]]),
                                                  length(defectiveScores$defectiveInts[[5]]))


names(num_codruggableInteractions_of_drug_subgroups) = phenotypes


correlation__numPatients__numCodruggableInteractions = cor(num_patients_in_drug_subgroups,num_codruggableInteractions_of_drug_subgroups)
correlation__numPatients__numCodruggableInteractions


DF_1drugs_2_numPatients_3_numCodruggableInteractions = data.frame(Treatment=phenotypes,
                                                                  NumPatients = num_patients_in_drug_subgroups,
                                                                  NumCodruggableInteractions = num_codruggableInteractions_of_drug_subgroups)




phenotypes_color_vec__ggplot = c('#BDCD00', '#0098A1', '#57AB27', '#612158', '#7A6FAC')


FigureWidth__geom_point = 5.5   # MR inserted
FigureHeight__geom_point = 3.5   # MR inserted

gg_DF_drugs_subgroups = ggplot(data=DF_1drugs_2_numPatients_3_numCodruggableInteractions,aes(NumPatients,NumCodruggableInteractions)) +           # MR modified
  geom_smooth(method=lm) + 
  geom_point(size=2,aes(colour=Treatment)) + 
  xlab('Number of treated patients') +           
  ylab('Number of co-druggable interactions') + 
  labs(title = paste('Correlation coefficient = ',round(correlation__numPatients__numCodruggableInteractions,digits = 2))) +
  #title(paste('Correlation coefficient = ',correlation__numPatients__numCodruggableInteractions))
  theme_bw()+
  scale_color_manual(values=phenotypes_color_vec__ggplot)
ggsave(file.path(figure_folder_of_current_script_,paste("geom_point__x__numPatients__y__numCodruggableInteractions.pdf",sep = "")), width = FigureWidth__geom_point, height = FigureHeight__geom_point)             # MR inserted
















#***********************************************************************
# *************** calculate similarity between median/mean subgroup networks and median/mean healthy network
#***********************************************************************

# H=defectiveScores$healthyNW
# MS=defectiveScores$MSuntreatedNw
# allDrugNws=defectiveScores$allDrugNws
# 
# 


IdxHealthy = which(annot2$Group == 'Healthy')
num_patientsHealthy=length(IdxHealthy)    # MR inserted

IdxMSuntreated = which(annot2$Treatment == "no" & annot2$Disease.Subtype!='')
num_patientsMSuntreated=length(IdxMSuntreated)    # MR inserted

Network__subgroup__H = phenotypeNetwork(IdxHealthy, allMedianNetworks,model,mode=groupingMode)  
# identical(Network__subgroup__H,H)
# [1] TRUE



Idx__H_untreatedMS_treatedMS= which(annot2$Treatment=='no' | annot2$Treatment == phenotypes[1] | annot2$Treatment == phenotypes[2] | annot2$Treatment == phenotypes[3] | annot2$Treatment == phenotypes[4] | annot2$Treatment == phenotypes[5])
Network__subgroup__H_untreatedMS_treatedMS= phenotypeNetwork(Idx__H_untreatedMS_treatedMS, allMedianNetworks,model,mode=groupingMode)  

Idx__treatedMS= which(annot2$Treatment == phenotypes[1] | annot2$Treatment == phenotypes[2] | annot2$Treatment == phenotypes[3] | annot2$Treatment == phenotypes[4] | annot2$Treatment == phenotypes[5])
Network__subgroup__treatedMS= phenotypeNetwork(Idx__treatedMS, allMedianNetworks,model,mode=groupingMode)  

Idx__untreatedMS= which(annot2$Treatment == "no" & annot2$Disease.Subtype!='')
#IdxMSuntreated = (annot$Treatment == "no" & annot$Disease.Subtype!='')
Network__subgroup__untreatedMS= phenotypeNetwork(Idx__untreatedMS, allMedianNetworks,model,mode=groupingMode)  


Idx__untreatedMS_treatedMS= c(Idx__untreatedMS,Idx__treatedMS)
Network__subgroup__untreatedMS_treatedMS= phenotypeNetwork(Idx__untreatedMS_treatedMS, allMedianNetworks,model,mode=groupingMode)  



Jaccard__H_treatedMS=jaccNello(H$network,Network__subgroup__treatedMS$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
# > identical(H$network,Network__subgroup__treatedMS$network)
# [1] TRUE


Jaccard__untreatedMS_treatedMS=jaccNello(Network__subgroup__untreatedMS$network,Network__subgroup__treatedMS$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)


Jaccard__all_versus_H=jaccNello(H$network,Network__subgroup__H_untreatedMS_treatedMS$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
Jaccard__all_versus_Gilenya=jaccNello(Network__subgroup__H_untreatedMS_treatedMS$network,allDrugNws[[1]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
Jaccard__all_versus_IFNb=jaccNello(Network__subgroup__H_untreatedMS_treatedMS$network,allDrugNws[[2]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
Jaccard__all_versus_Copaxone=jaccNello(Network__subgroup__H_untreatedMS_treatedMS$network,allDrugNws[[3]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
Jaccard__all_versus_EGCG=jaccNello(Network__subgroup__H_untreatedMS_treatedMS$network,allDrugNws[[4]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
Jaccard__all_versus_Tysabri=jaccNello(Network__subgroup__H_untreatedMS_treatedMS$network,allDrugNws[[5]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)




Jaccard__H__versus__untreatedMS_treatedMS=jaccNello(H$network,Network__subgroup__untreatedMS_treatedMS$network) 





Jaccard__H_untreatedMS=jaccNello(H$network,MS$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)

Jaccard__H_Gilenya=jaccNello(H$network,allDrugNws[[1]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
Jaccard__H_IFNb=jaccNello(H$network,allDrugNws[[2]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
Jaccard__H_Copaxone=jaccNello(H$network,allDrugNws[[3]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
Jaccard__H_EGCG=jaccNello(H$network,allDrugNws[[4]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
Jaccard__H_Tysabri=jaccNello(H$network,allDrugNws[[5]]$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)


Jaccard_vector___H__versus__treatedMS = c(Jaccard__H_Gilenya,
                                          Jaccard__H_IFNb,
                                          Jaccard__H_Copaxone,
                                          Jaccard__H_EGCG,
                                          Jaccard__H_Tysabri)


names(Jaccard_vector___H__versus__treatedMS) = phenotypes


DF_1drugs_2_numPatients_3_Jaccard = data.frame(Treatment=phenotypes,
                                               NumPatients = num_patients_in_drug_subgroups,
                                               JaccardSimilarity = Jaccard_vector___H__versus__treatedMS)



correlation__numPatients__Jaccard = cor(num_patients_in_drug_subgroups,Jaccard_vector___H__versus__treatedMS)
correlation__numPatients__Jaccard



FigureWidth__geom_point = 5.5   # MR inserted
FigureHeight__geom_point = 3.5   # MR inserted

#phenotypes_color_vec__ggplot = c('#612158', '#7A6FAC', '#57AB27', '#BDCD00', '#0098A1')


phenotypes_color_vec__ggplot = c('#BDCD00', '#0098A1', '#57AB27', '#612158', '#7A6FAC')


gg_DF_drugs_Jaccard = ggplot(data=DF_1drugs_2_numPatients_3_Jaccard,aes(NumPatients,JaccardSimilarity)) +           # MR modified
  geom_smooth(method=lm) + 
  geom_point(size=2,aes(colour=Treatment)) + 
  xlab('Number of treated patients') +           
  ylab('Jaccard similarity index') + 
  ylim(0.5,1.1)+
  labs(title = paste('Correlation coefficient = ',round(correlation__numPatients__Jaccard,digits = 2))) +
  #title(paste('Correlation coefficient = ',correlation__numPatients__numCodruggableInteractions))
  theme_bw()+
  scale_color_manual(values=phenotypes_color_vec__ggplot)
ggsave(file.path(figure_folder_of_current_script_,paste("geom_point__ylim__x__numPatients__y__JaccardSimilarity.pdf",sep = "")), width = FigureWidth__geom_point, height = FigureHeight__geom_point)             # MR inserted




numDefectiveDF=data.frame(Drug=character(length=numDrugs),
                          NumDefective=vector(length=numDrugs))

numDefectiveDF$Drug = phenotypes
numDefectiveDF$NumDefective = num_codruggableInteractions_of_drug_subgroups



FigureWidth__geom_bar = 3
FigureHeight__geom_bar = 3




#To prevent ggplot2 from reordering alphabetically the labels according to the factors, which are alphabetical
#We specify that the factors need to be ordered as they already are
# numDefectiveDF$Drug=factor(numDefectiveDF$Drug,levels=numDefectiveDF$Drug)
# manual re-naming
# numDefectiveDF$drug_name = c('FTY', 'IFNb', 'GA', 'EGCG', 'NTZ')
# numDefectiveDF$drug_name = factor(numDefectiveDF$drug_name, levels=c('IFNb', 'NTZ', 'FTY', 'GA', 'EGCG'))
# Barplot with the number of defective interaction for each drug, serves as legend for the main part of the figure
a = ggplot(data=numDefectiveDF,
           aes(Drug,NumDefective,fill=Drug)) + 
  geom_bar(stat='identity') + 
  theme_classic() + theme(legend.position = 'none') + 
  xlab('') + ylab('Number of \n co-druggable reactions') +
  theme(axis.title.y = element_text(size=8)) +
  scale_fill_manual(values=phenotypes_color_vec__ggplot, guide=FALSE) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(figure_folder_of_current_script_,paste("geom_bar__x__drug_name__y__NumDefective.pdf",sep = "")), width = FigureWidth__geom_bar, height = FigureHeight__geom_bar)            # MR inserted





#To prevent ggplot2 from reordering alphabetically the labels according to the factors, which are alphabetical
#We specify that the factors need to be ordered as they already are
# numDefectiveDF$Drug=factor(numDefectiveDF$Drug,levels=numDefectiveDF$Drug)
# manual re-naming
numDefectiveDF$drug_name = c('FTY', 'IFNb', 'GA', 'EGCG', 'NTZ')
numDefectiveDF$drug_name = factor(numDefectiveDF$drug_name, levels=c('IFNb', 'NTZ', 'FTY', 'GA', 'EGCG'))
# Barplot with the number of defective interaction for each drug, serves as legend for the main part of the figure
a2 = ggplot(data=numDefectiveDF,
            aes(drug_name,NumDefective,fill=drug_name)) + 
  geom_bar(stat='identity') + 
  theme_classic() + theme(legend.position = 'none') + 
  xlab('') + ylab('Number of \n co-druggable reactions') +
  theme(axis.title.y = element_text(size=8)) +
  scale_fill_manual(values=c('#612158', '#7A6FAC', '#57AB27', '#BDCD00', '#0098A1'), guide=FALSE) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(figure_folder_of_current_script_,paste("geom_bar__x__drug_name__y__NumDefective__2.pdf",sep = "")), width = FigureWidth__geom_bar, height = FigureHeight__geom_bar)









numTreatedPatientsDF=data.frame(Drug=character(length=numDrugs),
                                NumTreatedPatients=vector(length=numDrugs))

numTreatedPatientsDF$Drug = phenotypes
numTreatedPatientsDF$NumTreatedPatients = num_patients_in_subgroups[3:7]

#To prevent ggplot2 from reordering alphabetically the labels according to the factors, which are alphabetical
#We specify that the factors need to be ordered as they already are
# numTreatedPatientsDF$Drug=factor(numTreatedPatientsDF$Drug,levels=numTreatedPatientsDF$Drug)
# manual re-naming
numTreatedPatientsDF$drug_name = c('FTY', 'IFNb', 'GA', 'EGCG', 'NTZ')
numTreatedPatientsDF$drug_name = factor(numTreatedPatientsDF$drug_name, levels=c('IFNb', 'NTZ', 'FTY', 'GA', 'EGCG'))
# Barplot with the number of defective interaction for each drug, serves as legend for the main part of the figure
numTreatedPatientsDF_2 = ggplot(data=numTreatedPatientsDF,
                                aes(drug_name,NumTreatedPatients,fill=drug_name)) + 
  geom_bar(stat='identity') + 
  theme_classic() + theme(legend.position = 'none') + 
  xlab('') + ylab('Number of \n treated MS patients') +
  theme(axis.title.y = element_text(size=8)) +
  scale_fill_manual(values=c('#612158', '#7A6FAC', '#57AB27', '#BDCD00', '#0098A1'), guide=FALSE) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(figure_folder_of_current_script_,paste("geom_bar__x__drug_name__y__numTreatedPatients__2.pdf",sep = "")), width = FigureWidth__geom_bar, height = FigureHeight__geom_bar)













# Testing Jaccard similarity between ______________________________________________________________________________________________________________________
# 
# (1) Network__subgroup__H
#  VERSUS
# (2) randomly selected treatedMS patients with different subgroup size




number_of_random_samples = 10000



#matrix__ROWS_random_samples__COLS_Idx__random_subgroup_treatedMS = matrix(data=NA,nrow = as.numeric(number_of_random_samples),ncol = as.numeric(random_subgroup_size))
matrix__ROWS_random_samples__COLS_Jaccard__H_RandomTreatedMS_with_same_random_subgroup_size_as_drug_subgroups= matrix(data=NA,nrow = as.numeric(number_of_random_samples),ncol = numDrugs)
colnames(matrix__ROWS_random_samples__COLS_Jaccard__H_RandomTreatedMS_with_same_random_subgroup_size_as_drug_subgroups) = phenotypes

for(drug_ID in 1:numDrugs){
  
  random_subgroup_size = num_patients_in_drug_subgroups[drug_ID]
  
  
  for(index_random_sample in 1:number_of_random_samples){
    
    
    
    # Testing
    #sort(sample(1:11, 10, replace=FALSE))
    #sort(sample(1:11, 10, replace=T))
    Idx__random_subgroup_treatedMS = sample(Idx__treatedMS, random_subgroup_size, replace=FALSE)
    Network__random_subgroup_treatedMS= phenotypeNetwork(Idx__random_subgroup_treatedMS, allMedianNetworks,model,mode=groupingMode)  
    
    
    # matrix__ROWS_random_samples__COLS_Idx__random_subgroup_treatedMS[index_random_sample,]=sample(Idx__treatedMS, random_subgroup_size, replace=FALSE)
    # Network__random_subgroup_treatedMS= phenotypeNetwork(matrix__ROWS_random_samples__COLS_Idx__random_subgroup_treatedMS[index_random_sample,], allMedianNetworks,model,mode=groupingMode)  
    
    matrix__ROWS_random_samples__COLS_Jaccard__H_RandomTreatedMS_with_same_random_subgroup_size_as_drug_subgroups[index_random_sample,drug_ID]=jaccNello(H$network,Network__random_subgroup_treatedMS$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
    
    
  }
  
}

FigureWidth_hist = 5
FigureHeight_hist = 4
line_strength = 4

# drug_name_vec =c('IFNb', 'Tysabri', 'Gilenya', 'Copaxone', 'EGCG')
# drug_color_vec=c('#612158', '#7A6FAC', '#57AB27', '#BDCD00', '#0098A1')



phenotypes_color_vec = c('#57AB27','#612158','#BDCD00','#0098A1','#7A6FAC')


for(drug_ID in 1:numDrugs){
  
  
  
  matrix_used_for_plotting = as.numeric(matrix__ROWS_random_samples__COLS_Jaccard__H_RandomTreatedMS_with_same_random_subgroup_size_as_drug_subgroups[,colnames(matrix__ROWS_random_samples__COLS_Jaccard__H_RandomTreatedMS_with_same_random_subgroup_size_as_drug_subgroups)==phenotypes[drug_ID]])
  
  
  mean_to_plot = mean( matrix_used_for_plotting)
  std_to_plot = sd( matrix_used_for_plotting )
  mean_minus_1std_to_plot = mean_to_plot - std_to_plot
  mean_plus_1std_to_plot = mean_to_plot + std_to_plot
  
  
  pdf(file.path(figure_folder_of_current_script_,paste("hist__Jaccard__H_RandomTreatedMS_with_same_random_subgroup_size_as_",phenotypes[drug_ID] ,"_drug_subgroups.pdf",sep = "")),width = FigureWidth_hist, height = FigureHeight_hist)
  
  hist(x=matrix_used_for_plotting,
       main=paste("Random subgroup size as ",phenotypes[drug_ID] , "=",num_patients_in_drug_subgroups[drug_ID] ,"\n (mean = ", round(mean_to_plot,digits = 2),", std = ", round(std_to_plot,digits = 2), ")",sep = ""),
       xlab=paste("Jaccard similarity index\n between healthy subgroup and random treatedMS subgroup",sep = ""),
       xlim = c(0,1),
       col=phenotypes_color_vec[drug_ID])
  # ,
  #      breaks=seq(0,72,2),xaxt="n",
  #      main=paste("HC all \n (#",number_of_HC,", mean = ", round(mean_to_plot,digits = 1),", std = ", round(std_to_plot,digits = 1), ")",sep=""),
  #      xlim = c(0,72.1),ylim = c(0,ylim_hist__age),
  #      col=color_HC_all)
  # 
  abline(v = mean_to_plot,
         col = 'black',
         lwd = line_strength)
  
  abline(v = mean_minus_1std_to_plot,
         col = 'black',
         lwd = line_strength,
         lty="dashed")
  
  abline(v = mean_plus_1std_to_plot,
         col = 'black',
         lwd = line_strength,
         lty="dashed")
  # 
  # axis(side=1, at=seq(0,72,4))
  
  dev.off()
  
  
}








































number_of_random_samples = 10000
different_random_subgroup_size_vec = c(37,50,60,69)
matrix__ROWS_random_samples__COLS_Jaccard__H_RandomTreatedMS_with_different_random_subgroup_size= matrix(data=NA,nrow = as.numeric(number_of_random_samples),ncol = length(different_random_subgroup_size_vec))
colnames(matrix__ROWS_random_samples__COLS_Jaccard__H_RandomTreatedMS_with_different_random_subgroup_size) = as.character(different_random_subgroup_size_vec)

for(ID in 1:length(different_random_subgroup_size_vec)){
  
  random_subgroup_size = different_random_subgroup_size_vec[ID]
  
  
  for(index_random_sample in 1:number_of_random_samples){
    
    
    
    # Testing
    #sort(sample(1:11, 10, replace=FALSE))
    #sort(sample(1:11, 10, replace=T))
    Idx__random_subgroup_treatedMS = sample(Idx__treatedMS, random_subgroup_size, replace=FALSE)
    Network__random_subgroup_treatedMS= phenotypeNetwork(Idx__random_subgroup_treatedMS, allMedianNetworks,model,mode=groupingMode)  
    
    
    # matrix__ROWS_random_samples__COLS_Idx__random_subgroup_treatedMS[index_random_sample,]=sample(Idx__treatedMS, random_subgroup_size, replace=FALSE)
    # Network__random_subgroup_treatedMS= phenotypeNetwork(matrix__ROWS_random_samples__COLS_Idx__random_subgroup_treatedMS[index_random_sample,], allMedianNetworks,model,mode=groupingMode)  
    
    matrix__ROWS_random_samples__COLS_Jaccard__H_RandomTreatedMS_with_different_random_subgroup_size[index_random_sample,ID]=jaccNello(H$network,Network__random_subgroup_treatedMS$network)     # MR note: Be careful, jaccNello calculates the SMC = simple matching coefficientn (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
    
    
  }
  
}



FigureWidth_hist = 5
FigureHeight_hist = 4
line_strength = 4

# drug_name_vec =c('IFNb', 'Tysabri', 'Gilenya', 'Copaxone', 'EGCG')
# drug_color_vec=c('#612158', '#7A6FAC', '#57AB27', '#BDCD00', '#0098A1')



#phenotypes_color_vec = c('#57AB27','#612158','#BDCD00','#0098A1','#7A6FAC')


for(ID in 1:length(different_random_subgroup_size_vec)){
  
  
  
  matrix_used_for_plotting = as.numeric(matrix__ROWS_random_samples__COLS_Jaccard__H_RandomTreatedMS_with_different_random_subgroup_size[,ID])
  
  
  mean_to_plot = mean( matrix_used_for_plotting)
  std_to_plot = sd( matrix_used_for_plotting )
  mean_minus_1std_to_plot = mean_to_plot - std_to_plot
  mean_plus_1std_to_plot = mean_to_plot + std_to_plot
  
  
  pdf(file.path(figure_folder_of_current_script_,paste("hist__Jaccard__H_RandomTreatedMS_with_random_subgroup_size_",different_random_subgroup_size_vec[ID] ,".pdf",sep = "")),width = FigureWidth_hist, height = FigureHeight_hist)
  
  hist(x=matrix_used_for_plotting,
       main=paste("Random subgroup size ",different_random_subgroup_size_vec[ID] ,"\n (mean = ", round(mean_to_plot,digits = 2),", std = ", round(std_to_plot,digits = 2), ")",sep = ""),
       xlab=paste("Jaccard similarity index\n between healthy subgroup and random treatedMS subgroup",sep = ""),
       xlim = c(0,1))
  # ,
  #      col=phenotypes_color_vec[ID])
  # ,
  #      breaks=seq(0,72,2),xaxt="n",
  #      main=paste("HC all \n (#",number_of_HC,", mean = ", round(mean_to_plot,digits = 1),", std = ", round(std_to_plot,digits = 1), ")",sep=""),
  #      xlim = c(0,72.1),ylim = c(0,ylim_hist__age),
  #      col=color_HC_all)
  # 
  abline(v = mean_to_plot,
         col = 'black',
         lwd = line_strength)
  
  abline(v = mean_minus_1std_to_plot,
         col = 'black',
         lwd = line_strength,
         lty="dashed")
  
  abline(v = mean_plus_1std_to_plot,
         col = 'black',
         lwd = line_strength,
         lty="dashed")
  # 
  # axis(side=1, at=seq(0,72,4))
  
  dev.off()
  
  
}











# __________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# Analysis of the probability that an interaction is inactive or active  ____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# __________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


annot_OK = annot2[Idx__H_untreatedMS_treatedMS,]
allMedianNetworks_OK = allMedianNetworks[Idx__H_untreatedMS_treatedMS,]

colSums_allMedianNetworks_OK = colSums(allMedianNetworks_OK)

count_interactions_active_OK = colSums_allMedianNetworks_OK
count_interactions_inactive_OK = nrow(allMedianNetworks_OK) - colSums_allMedianNetworks_OK

sum__count_interactions_active_inactive_OK = count_interactions_active_OK+count_interactions_inactive_OK
min(sum__count_interactions_active_inactive_OK)
max(sum__count_interactions_active_inactive_OK)



matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK = rbind(count_interactions_active_OK,count_interactions_inactive_OK)


FigureWidth_barplot = 20
FigureHeight_barplot = 10
line_strength = 2

pdf(file.path(figure_folder_of_current_script_,paste("barplot__matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK.pdf",sep = "")),width = FigureWidth_barplot, height = FigureHeight_barplot)


num_all_interactions = ncol(allMedianNetworks_OK)
plot_row_number = 3

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(3, 0, 0, 0),
    mfrow=c(plot_row_number,1)) 


#par(mfrow=c(plot_row_number,1)) 
barplot(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,(1:(1*(num_all_interactions/plot_row_number)))],
        main = "Counting active and inactive interactions of allMedianNetworks_OK",
        ylab = "Number of donor specific networks",
        #xlab = "Interactions",
        col = c("green","orange"),
        beside = TRUE,las=2
)

# lablist<-as.vector(colnames(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,(1:(1*(num_all_interactions/plot_row_number)))]))
# axis(1, labels = FALSE)
# text(labels = lablist, srt = 45,xpd = TRUE)

abline(h = 37,
       col = 'black',
       lwd = line_strength)

barplot(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
        #main = "Counting active and inactive interactions of allMedianNetworks_OK",
        ylab = "Number of donor specific networks",
        #xlab = "Interactions",
        col = c("green","orange"),
        beside = TRUE,las=2
)

abline(h = 37,
       col = 'black',
       lwd = line_strength)

barplot(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
        #main = "Counting active and inactive interactions of allMedianNetworks_OK",
        ylab = "Number of donor specific networks",
        #xlab = "Interactions",
        col = c("green","orange"),
        beside = TRUE,las=2
)


abline(h = 37,
       col = 'black',
       lwd = line_strength)


legend("right",
       c("Active","Inactive"),
       fill = c("green","orange")
)

#text(srt=45, adj=1, xpd=TRUE)

dev.off()









































# __________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# Analysis of the probability that an interaction is inactive or active  
# based on drug subgroups
# ____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# __________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

# > phenotypes
# [1] "Gilenya"  "IFNb"     "Copaxone" "EGCG"     "Tysabri" 




IdxHealthy = which(annot2$Group == 'Healthy')
numHealthy=length(IdxHealthy)    # MR inserted

Healthy_DonorSpecificNetworks = allMedianNetworks[IdxHealthy,]


colSums_allMedianNetworks_of_HealthyDonors = colSums(Healthy_DonorSpecificNetworks)

count_interactions_active__HealthyDonors = colSums_allMedianNetworks_of_HealthyDonors
count_interactions_inactive__HealthyDonors = nrow(Healthy_DonorSpecificNetworks) - colSums_allMedianNetworks_of_HealthyDonors


matrix__ROW1_CountbInteractionsON__ROW2_CountInteractionsOFF__HealthyDonors = rbind(count_interactions_active__HealthyDonors,count_interactions_inactive__HealthyDonors)
matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__HealthyDonors = matrix__ROW1_CountbInteractionsON__ROW2_CountInteractionsOFF__HealthyDonors / nrow(Healthy_DonorSpecificNetworks)


phenotypes_color_vec = c('#57AB27','#612158','#BDCD00','#0098A1','#7A6FAC')

FigureWidth_barplot = 20
FigureHeight_barplot = 10
line_strength = 2






# drug_no_i = 1
# drug_name_i = phenotypes[drug_no_i]
for(drug_no_i in 1:numDrugs){
  
  
  drug_name_i = phenotypes[drug_no_i]
  
  
  
  
  Idx__MS_treated_with_drug_no_i= which(annot2$Treatment == phenotypes[drug_no_i])
  Network__subgroup__MS_treated_with_drug_no_i= phenotypeNetwork(Idx__MS_treated_with_drug_no_i, allMedianNetworks,model,mode=groupingMode)  
  
  
  allMedianNetworks_of_MS_treated_with_drug_no_i = allMedianNetworks[Idx__MS_treated_with_drug_no_i,]
  
  
  
  colSums_allMedianNetworks_of_MS_treated_with_drug_no_i = colSums(allMedianNetworks_of_MS_treated_with_drug_no_i)
  
  count_interactions_active__MS_treated_with_drug_no_i = colSums_allMedianNetworks_of_MS_treated_with_drug_no_i
  count_interactions_inactive__MS_treated_with_drug_no_i = nrow(allMedianNetworks_of_MS_treated_with_drug_no_i) - colSums_allMedianNetworks_of_MS_treated_with_drug_no_i
  
  # for testing
  # sum__count_interactions_active_inactive__MS_treated_with_drug_no_i = count_interactions_active__MS_treated_with_drug_no_i + count_interactions_inactive__MS_treated_with_drug_no_i
  # min(sum__count_interactions_active_inactive__MS_treated_with_drug_no_i)
  # max(sum__count_interactions_active_inactive__MS_treated_with_drug_no_i)
  
  
  
  matrix__ROW1_count_interactions_active__ROW2_count_interactions_inactive___of__MS_treated_with_drug_no_i = rbind(count_interactions_active__MS_treated_with_drug_no_i,count_interactions_inactive__MS_treated_with_drug_no_i)
  matrix__ROW1_probability_interactions_active__ROW2_probability_interactions_inactive___of__MS_treated_with_drug_no_i = matrix__ROW1_count_interactions_active__ROW2_count_interactions_inactive___of__MS_treated_with_drug_no_i / nrow(allMedianNetworks_of_MS_treated_with_drug_no_i)
  
  
  
  
  
  
  
  FigureWidth_barplot = 20
  FigureHeight_barplot = 10
  line_strength = 2
  
  pdf(file.path(figure_folder_of_current_script_,paste("barplot__matrix__ROW1_count_interactions_active__ROW2_count_interactions_inactive___of__MS_treated_with_drug_", phenotypes[drug_no_i],".pdf",sep = "")),width = FigureWidth_barplot, height = FigureHeight_barplot)
  
  
  num_all_interactions = ncol(allMedianNetworks_of_MS_treated_with_drug_no_i)
  plot_row_number = 3
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(3, 0, 0, 0),
      mfrow=c(plot_row_number,1)) 
  
  
  #par(mfrow=c(plot_row_number,1)) 
  barplot(matrix__ROW1_count_interactions_active__ROW2_count_interactions_inactive___of__MS_treated_with_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
          main = paste("Counting active and inactive interactions of \nallMedianNetworks of MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
          ylab = "Number of donor specific networks",
          #xlab = "Interactions",
          col = c("green","orange"),
          beside = TRUE,las=2
  )
  
  # lablist<-as.vector(colnames(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,(1:(1*(num_all_interactions/plot_row_number)))]))
  # axis(1, labels = FALSE)
  # text(labels = lablist, srt = 45,xpd = TRUE)
  
  abline(h = sum(matrix__ROW1_count_interactions_active__ROW2_count_interactions_inactive___of__MS_treated_with_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  
  barplot(matrix__ROW1_count_interactions_active__ROW2_count_interactions_inactive___of__MS_treated_with_drug_no_i[,((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
          #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
          ylab = "Number of donor specific networks",
          #xlab = "Interactions",
          col = c("green","orange"),
          beside = TRUE,las=2
  )
  
  abline(h = sum(matrix__ROW1_count_interactions_active__ROW2_count_interactions_inactive___of__MS_treated_with_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  
  barplot(matrix__ROW1_count_interactions_active__ROW2_count_interactions_inactive___of__MS_treated_with_drug_no_i[,((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
          #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
          ylab = "Number of donor specific networks",
          #xlab = "Interactions",
          col = c("green","orange"),
          beside = TRUE,las=2
  )
  
  
  abline(h = sum(matrix__ROW1_count_interactions_active__ROW2_count_interactions_inactive___of__MS_treated_with_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  
  legend("right",
         c("Active","Inactive"),
         fill = c("green","orange")
  )
  
  #text(srt=45, adj=1, xpd=TRUE)
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  FigureWidth_barplot = 20
  FigureHeight_barplot = 10
  line_strength = 2
  
  pdf(file.path(figure_folder_of_current_script_,paste("barplot__matrix__ROW1_probability_interactions_active__ROW2_probability_interactions_inactive___of__MS_treated_with_drug_", phenotypes[drug_no_i],".pdf",sep = "")),width = FigureWidth_barplot, height = FigureHeight_barplot)
  
  
  num_all_interactions = ncol(allMedianNetworks_of_MS_treated_with_drug_no_i)
  plot_row_number = 3
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(3, 0, 0, 0),
      mfrow=c(plot_row_number,1)) 
  
  
  #par(mfrow=c(plot_row_number,1)) 
  barplot(matrix__ROW1_probability_interactions_active__ROW2_probability_interactions_inactive___of__MS_treated_with_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
          main = paste("Probability of active and inactive interactions based on \nallMedianNetworks of MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
          ylab = "Drug-subgroup probability",
          #xlab = "Interactions",
          col = c("green","orange"),
          beside = TRUE,las=2
  )
  
  # lablist<-as.vector(colnames(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,(1:(1*(num_all_interactions/plot_row_number)))]))
  # axis(1, labels = FALSE)
  # text(labels = lablist, srt = 45,xpd = TRUE)
  
  abline(h = sum(matrix__ROW1_probability_interactions_active__ROW2_probability_interactions_inactive___of__MS_treated_with_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  barplot(matrix__ROW1_probability_interactions_active__ROW2_probability_interactions_inactive___of__MS_treated_with_drug_no_i[,((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
          #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
          ylab = "Drug-subgroup probability",
          #xlab = "Interactions",
          col = c("green","orange"),
          beside = TRUE,las=2
  )
  
  abline(h = sum(matrix__ROW1_probability_interactions_active__ROW2_probability_interactions_inactive___of__MS_treated_with_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  barplot(matrix__ROW1_probability_interactions_active__ROW2_probability_interactions_inactive___of__MS_treated_with_drug_no_i[,((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
          #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
          ylab = "Drug-subgroup probability",
          #xlab = "Interactions",
          col = c("green","orange"),
          beside = TRUE,las=2
  )
  
  
  abline(h = sum(matrix__ROW1_probability_interactions_active__ROW2_probability_interactions_inactive___of__MS_treated_with_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  legend("right",
         c("Active","Inactive"),
         fill = c("green","orange")
  )
  
  #text(srt=45, adj=1, xpd=TRUE)
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  drug_name_i = phenotypes[drug_no_i]
  
  
  
  
  Idx__MS_treated_with_drug_no_i= which(annot2$Treatment == phenotypes[drug_no_i])
  #Network__subgroup__MS_treated_with_drug_no_i= phenotypeNetwork(Idx__MS_treated_with_drug_no_i, allMedianNetworks,model,mode=groupingMode)  
  
  
  allMedianNetworks_of_MS_treated_with_drug_no_i = allMedianNetworks[Idx__MS_treated_with_drug_no_i,]
  
  
  
  colSums_allMedianNetworks_of_MS_treated_with_drug_no_i = colSums(allMedianNetworks_of_MS_treated_with_drug_no_i)
  
  count_interactions_active__MS_treated_with_drug_no_i = colSums_allMedianNetworks_of_MS_treated_with_drug_no_i
  count_interactions_inactive__MS_treated_with_drug_no_i = nrow(allMedianNetworks_of_MS_treated_with_drug_no_i) - colSums_allMedianNetworks_of_MS_treated_with_drug_no_i
  
  # for testing
  # sum__count_interactions_active_inactive__MS_treated_with_drug_no_i = count_interactions_active__MS_treated_with_drug_no_i + count_interactions_inactive__MS_treated_with_drug_no_i
  # min(sum__count_interactions_active_inactive__MS_treated_with_drug_no_i)
  # max(sum__count_interactions_active_inactive__MS_treated_with_drug_no_i)
  
  
  
  matrix__ROW1_CountInteractionsON__ROW2_CountInteractionsOFF__MStreated_drug_no_i = rbind(count_interactions_active__MS_treated_with_drug_no_i,count_interactions_inactive__MS_treated_with_drug_no_i)
  matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__MStreated_drug_no_i = matrix__ROW1_CountInteractionsON__ROW2_CountInteractionsOFF__MStreated_drug_no_i / nrow(allMedianNetworks_of_MS_treated_with_drug_no_i)
  
  matrix__ROW1_CountInteractionsOFF__ROW2_CountInteractionsON__MStreated_drug_no_i = rbind(count_interactions_inactive__MS_treated_with_drug_no_i,count_interactions_active__MS_treated_with_drug_no_i)
  matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i = matrix__ROW1_CountInteractionsOFF__ROW2_CountInteractionsON__MStreated_drug_no_i / nrow(allMedianNetworks_of_MS_treated_with_drug_no_i)
  
  
  matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i = rbind(matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__HealthyDonors[1,],
                                                                                                       matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__MStreated_drug_no_i[1,])
  
  
  
  
  CoDruggableInteractions_predicted_for_MStreated_drug_no_i = names(defectiveScores$defectiveInts[[drug_no_i]])
  vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i = NA*matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i[1,]
  vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[CoDruggableInteractions_predicted_for_MStreated_drug_no_i]=-0.1
  #Indices_CoDruggableInteractions_predicted_for_MStreated_drug_no_i = which(colnames(matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i))
  
  
  
  
  # 
  # pdf(file.path(figure_folder_of_current_script_,paste("barplot__matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__MStreated_drug_", phenotypes[drug_no_i],".pdf",sep = "")),width = FigureWidth_barplot, height = FigureHeight_barplot)
  # 
  # 
  # num_all_interactions = ncol(allMedianNetworks_of_MS_treated_with_drug_no_i)
  # plot_row_number = 3
  # 
  # mar.default <- c(5,4,4,2) + 0.1
  # par(mar = mar.default + c(3, 0, 0, 0),
  #     mfrow=c(plot_row_number,1)) 
  # 
  # 
  # #par(mfrow=c(plot_row_number,1)) 
  # #
  # 
  # 
  # # barplot(vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
  # #         main = paste("Probability of active and inactive interactions based on \nallMedianNetworks of MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
  # #         ylab = "Drug-subgroup probability",
  # #         ylim = c(0,1.1),
  # #         #xlab = "Interactions",
  # #         col = c("green","orange"),
  # #         beside = TRUE,las=2
  # # )
  # 
  # barplot(matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__MStreated_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
  #         main = paste("Probability of active and inactive interactions based on \nallMedianNetworks of MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
  #         ylab = "Drug-subgroup probability",
  #         ylim = c(0,1.1),
  #         #xlab = "Interactions",
  #         col = c("green","orange"),
  #         beside = TRUE,las=2
  # )
  # 
  # # lablist<-as.vector(colnames(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,(1:(1*(num_all_interactions/plot_row_number)))]))
  # # axis(1, labels = FALSE)
  # # text(labels = lablist, srt = 45,xpd = TRUE)
  # 
  # abline(h = sum(matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__MStreated_drug_no_i[,1]),
  #        col = 'black',
  #        lwd = line_strength)
  # 
  # abline(h = 0.5,
  #        col = 'black',
  #        lwd = line_strength)
  # 
  # barplot(matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__MStreated_drug_no_i[,((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
  #         #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
  #         ylab = "Drug-subgroup probability",
  #         ylim = c(0,1.1),
  #         #xlab = "Interactions",
  #         col = c("green","orange"),
  #         beside = TRUE,las=2
  # )
  # 
  # abline(h = sum(matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__MStreated_drug_no_i[,1]),
  #        col = 'black',
  #        lwd = line_strength)
  # 
  # abline(h = 0.5,
  #        col = 'black',
  #        lwd = line_strength)
  # 
  # barplot(matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__MStreated_drug_no_i[,((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
  #         #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
  #         ylab = "Drug-subgroup probability",
  #         ylim = c(0,1.1),
  #         #xlab = "Interactions",
  #         col = c("green","orange"),
  #         beside = TRUE,las=2
  # )
  # 
  # 
  # abline(h = sum(matrix__ROW1_ProbInteractionsON__ROW2_ProbInteractionsOFF__MStreated_drug_no_i[,1]),
  #        col = 'black',
  #        lwd = line_strength)
  # 
  # abline(h = 0.5,
  #        col = 'black',
  #        lwd = line_strength)
  # 
  # legend("topright",
  #        c("ON-probability","OFF-probability"),
  #        fill = c("green","orange")
  # )
  # 
  # #text(srt=45, adj=1, xpd=TRUE)
  # 
  # dev.off()
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # pdf(file.path(figure_folder_of_current_script_,paste("barplot__matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_", phenotypes[drug_no_i],".pdf",sep = "")),width = FigureWidth_barplot, height = FigureHeight_barplot)
  # 
  # 
  # num_all_interactions = ncol(allMedianNetworks_of_MS_treated_with_drug_no_i)
  # plot_row_number = 3
  # 
  # 
  # 
  # 
  # #par(mfrow=c(plot_row_number,1)) 
  # #
  # #
  # mar.default <- c(5,4,4,2) + 0.1
  # par(mar = mar.default + c(3, 0, 0, 0),
  #     mfrow=c(plot_row_number,1))
  # #col.axis = colorlist[labs %in% redlabs +1 ]) 
  # #col.axis = rep("red",186))
  # 
  # 
  # bp1 = barplot(matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
  #               main = paste("Probability of active interactions based on \nallMedianNetworks of healthy donors versus MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
  #               ylab = "ON-probability",
  #               ylim = c(-0.15,1.05),
  #               #xlab = "Interactions",
  #               col = c("grey",phenotypes_color_vec[drug_no_i]),
  #               beside = TRUE,las=2
  # )
  # 
  # points(x = (bp1[1,]+bp1[2,])/2, y = vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[(1:(1*(num_all_interactions/plot_row_number)))],
  #        col = "blue", cex = 2,pch = 23,
  #        bg="blue", lwd=1)
  # 
  # 
  # # text(bp1, 0.9, 1:62,cex=1,pos=3) 
  # # # draw asterics
  # # text(bp1, x[1]+((x[2]-x[1])/2),y+offset,"**")
  # # 
  # # labs = colnames(matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))])
  # # redlabs = CoDruggableInteractions_predicted_for_MStreated_drug_no_i
  # # colorlist = c("black","red")
  # # # one of many ways to generate the color labels
  # # #axiscolor = colorlist[labs %in% redlabs +1 ]
  # # axis(1, col.axis = colorlist[labs %in% redlabs +1 ])
  # # 
  # # lablist<-as.vector(colnames(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,(1:(1*(num_all_interactions/plot_row_number)))]))
  # # axis(1, labels = FALSE)
  # # text(labels = lablist, srt = 45,xpd = TRUE)
  # 
  # 
  # 
  # abline(h = 0.5,
  #        col = 'black',
  #        lwd = line_strength)
  # 
  # bp2 = barplot(matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i[,((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
  #               #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
  #               ylab = "ON-probability",
  #               ylim = c(-0.15,1.05),
  #               #xlab = "Interactions",
  #               col = c("grey",phenotypes_color_vec[drug_no_i]),
  #               beside = TRUE,las=2
  # )
  # 
  # points(x = (bp2[1,]+bp2[2,])/2, y = vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
  #        col = "blue", cex = 2,pch = 23,
  #        bg="blue", lwd=1)
  # 
  # 
  # 
  # abline(h = 0.5,
  #        col = 'black',
  #        lwd = line_strength)
  # 
  # bp3 = barplot(matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i[,((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
  #               #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
  #               ylab = "ON-probability",
  #               ylim = c(-0.15,1.05),
  #               #xlab = "Interactions",
  #               col = c("grey",phenotypes_color_vec[drug_no_i]),
  #               beside = TRUE,las=2
  # )
  # 
  # points(x = (bp3[1,]+bp3[2,])/2, y = vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
  #        col = "blue", cex = 2,pch = 23,
  #        bg="blue", lwd=1)
  # 
  # 
  # 
  # abline(h = 0.5,
  #        col = 'black',
  #        lwd = line_strength)
  # 
  # 
  # 
  # 
  # 
  # legend("topright",
  #        c("ON-probability of healthy donors",paste("ON-probability of ", phenotypes[drug_no_i] , " treated MS patients",sep="")),
  #        fill = c("grey",phenotypes_color_vec[drug_no_i])
  # )
  # 
  # #text(srt=45, adj=1, xpd=TRUE)
  # 
  # dev.off()
  # 
  # 
  # 
  
  
  
  
  
  
  pdf(file.path(figure_folder_of_current_script_,paste("barplot__matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_", phenotypes[drug_no_i],".pdf",sep = "")),width = FigureWidth_barplot, height = FigureHeight_barplot)
  
  
  num_all_interactions = ncol(allMedianNetworks_of_MS_treated_with_drug_no_i)
  plot_row_number = 3
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(3, 0, 0, 0),
      mfrow=c(plot_row_number,1)) 
  
  
  #par(mfrow=c(plot_row_number,1)) 
  #
  
  
  # barplot(vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
  #         main = paste("Probability of active and inactive interactions based on \nallMedianNetworks of MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
  #         ylab = "Drug-subgroup probability",
  #         ylim = c(0,1.1),
  #         #xlab = "Interactions",
  #         col = c("orange",phenotypes_color_vec[drug_no_i]),
  #         beside = TRUE,las=2
  # )
  
  barplot(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
          main = paste("Probability of active and inactive interactions based on \nallMedianNetworks of MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
          ylab = "Drug-subgroup probability",
          ylim = c(0,1.1),
          #xlab = "Interactions",
          col = c("orange",phenotypes_color_vec[drug_no_i]),
          beside = TRUE,las=2
  )
  
  # lablist<-as.vector(colnames(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,(1:(1*(num_all_interactions/plot_row_number)))]))
  # axis(1, labels = FALSE)
  # text(labels = lablist, srt = 45,xpd = TRUE)
  
  abline(h = sum(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  barplot(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
          #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
          ylab = "Drug-subgroup probability",
          ylim = c(0,1.1),
          #xlab = "Interactions",
          col = c("orange",phenotypes_color_vec[drug_no_i]),
          beside = TRUE,las=2
  )
  
  abline(h = sum(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  barplot(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
          #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
          ylab = "Drug-subgroup probability",
          ylim = c(0,1.1),
          #xlab = "Interactions",
          col = c("orange",phenotypes_color_vec[drug_no_i]),
          beside = TRUE,las=2
  )
  
  
  abline(h = sum(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  legend("left",
         c("OFF","ON"),
         fill = c("orange",phenotypes_color_vec[drug_no_i])
  )
  
  #text(srt=45, adj=1, xpd=TRUE)
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  pdf(file.path(figure_folder_of_current_script_,paste("barplot_HighlightCoDruggablePrediction__matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_", phenotypes[drug_no_i],".pdf",sep = "")),width = FigureWidth_barplot, height = FigureHeight_barplot)
  
  
  num_all_interactions = ncol(allMedianNetworks_of_MS_treated_with_drug_no_i)
  plot_row_number = 3
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(3, 0, 0, 0),
      mfrow=c(plot_row_number,1)) 
  
  
  #par(mfrow=c(plot_row_number,1)) 
  #
  
  
  # barplot(vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
  #         main = paste("Probability of active and inactive interactions based on \nallMedianNetworks of MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
  #         ylab = "Drug-subgroup probability",
  #         ylim = c(0,1.1),
  #         #xlab = "Interactions",
  #         col = c("orange",phenotypes_color_vec[drug_no_i]),
  #         beside = TRUE,las=2
  # )
  
  b1 = barplot(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
               main = paste("Probability of active and inactive interactions based on \nallMedianNetworks of MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
               ylab = "Drug-subgroup probability",
               ylim = c(-0.15,1.05),
               #xlab = "Interactions",
               col = c("orange",phenotypes_color_vec[drug_no_i]),
               beside = TRUE,las=2
  )
  
  
  points(x = (b1[1,]+b1[2,])/2, y = vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[(1:(1*(num_all_interactions/plot_row_number)))],
         col = "blue", cex = 2,pch = 23,
         bg="blue", lwd=1)
  
  # lablist<-as.vector(colnames(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,(1:(1*(num_all_interactions/plot_row_number)))]))
  # axis(1, labels = FALSE)
  # text(labels = lablist, srt = 45,xpd = TRUE)
  
  abline(h = sum(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  b2 = barplot(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
               #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
               ylab = "Drug-subgroup probability",
               ylim = c(-0.15,1.05),
               #xlab = "Interactions",
               col = c("orange",phenotypes_color_vec[drug_no_i]),
               beside = TRUE,las=2
  )
  
  points(x = (b2[1,]+b2[2,])/2, y = vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
         col = "blue", cex = 2,pch = 23,
         bg="blue", lwd=1)
  
  abline(h = sum(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  b3 = barplot(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
               #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
               ylab = "Drug-subgroup probability",
               ylim = c(-0.15,1.05),
               #xlab = "Interactions",
               col = c("orange",phenotypes_color_vec[drug_no_i]),
               beside = TRUE,las=2
  )
  
  points(x = (b3[1,]+b3[2,])/2, y = vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
         col = "blue", cex = 2,pch = 23,
         bg="blue", lwd=1)
  
  abline(h = sum(matrix__ROW1_ProbInteractionsOFF__ROW2_ProbInteractionsON__MStreated_drug_no_i[,1]),
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  legend("left",
         c("OFF","ON"),
         fill = c("orange",phenotypes_color_vec[drug_no_i])
  )
  
  #text(srt=45, adj=1, xpd=TRUE)
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  pdf(file.path(figure_folder_of_current_script_,paste("barplot__matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_", phenotypes[drug_no_i],".pdf",sep = "")),width = FigureWidth_barplot, height = FigureHeight_barplot)
  
  
  num_all_interactions = ncol(allMedianNetworks_of_MS_treated_with_drug_no_i)
  plot_row_number = 3
  
  
  
  #labs = colnames(matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i)
  
  
  #par(mfrow=c(plot_row_number,1)) 
  #
  #
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(3, 0, 0, 0),
      mfrow=c(plot_row_number,1))
  #col.axis = colorlist[labs %in% redlabs +1 ]) 
  #col.axis = rep("red",186))
  
  
  bp1 = barplot(matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i[,(1:(1*(num_all_interactions/plot_row_number)))],
                main = paste("Probability of active interactions based on \nallMedianNetworks of healthy donors versus MS patients treated with drug ", phenotypes[drug_no_i] , sep=""),
                ylab = "ON-probability",
                ylim = c(-0.15,1.05),
                #xlab = "Interactions",
                col = c("grey",phenotypes_color_vec[drug_no_i]),
                beside = TRUE,las=2
  )
  
  points(x = (bp1[1,]+bp1[2,])/2, y = vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[(1:(1*(num_all_interactions/plot_row_number)))],
         col = "blue", cex = 2,pch = 23,
         bg="blue", lwd=1)
  
  
  # text(bp1, 0.9, 1:62,cex=1,pos=3) 
  # # draw asterics
  # text(bp1, x[1]+((x[2]-x[1])/2),y+offset,"**")
  # 
  
  # colorlist = c("black","red")
  # # one of many ways to generate the color labels
  # #axiscolor = colorlist[labs %in% redlabs +1 ]
  # axis(1, col.axis = colorlist[labs %in% redlabs +1 ])
  # 
  # lablist<-as.vector(colnames(matrix__ROW1_count_interactions_active_OK__ROW2_count_interactions_inactive_OK[,(1:(1*(num_all_interactions/plot_row_number)))]))
  # axis(1, labels = FALSE)
  # text(labels = lablist, srt = 45,xpd = TRUE)
  
  abline(h = 1,
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  bp2 = barplot(matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i[,((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
                #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
                ylab = "ON-probability",
                ylim = c(-0.15,1.05),
                #xlab = "Interactions",
                col = c("grey",phenotypes_color_vec[drug_no_i]),
                beside = TRUE,las=2
  )
  
  points(x = (bp2[1,]+bp2[2,])/2, y = vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[((1*(num_all_interactions/plot_row_number)+1) : (2*(num_all_interactions/plot_row_number)))],
         col = "blue", cex = 2,pch = 23,
         bg="blue", lwd=1)
  
  abline(h = 1,
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  bp3 = barplot(matrix__ROW1_ProbInteractionsON__HealthyDonors__ROW2_ProbInteractionsON__MStreated_drug_no_i[,((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
                #main = "Counting active and inactive interactions of allMedianNetworks_of_MS_treated_with_drug_no_i",
                ylab = "ON-probability",
                ylim = c(-0.15,1.05),
                #xlab = "Interactions",
                col = c("grey",phenotypes_color_vec[drug_no_i]),
                beside = TRUE,las=2
  )
  
  points(x = (bp3[1,]+bp3[2,])/2, y = vec_logical_CoDruggableInteractions_predicted_for_MStreated_drug_no_i[((2*(num_all_interactions/plot_row_number)+1) : num_all_interactions )],
         col = "blue", cex = 2,pch = 23,
         bg="blue", lwd=1)
  
  abline(h = 1,
         col = 'black',
         lwd = line_strength)
  
  abline(h = 0.5,
         col = 'black',
         lwd = line_strength)
  
  
  
  
  
  legend("topright",
         c("ON-probability of healthy donors",paste("ON-probability of ", phenotypes[drug_no_i] , " treated MS patients",sep="")),
         fill = c("grey",phenotypes_color_vec[drug_no_i])
  )
  
  #text(srt=45, adj=1, xpd=TRUE)
  
  dev.off()
  
  
  
  
  
  
  
}















prob_interactions_active_OK = count_interactions_active_OK / nrow(allMedianNetworks_OK)
sort_prob_interactions_active_OK = sort(prob_interactions_active_OK)

unlist__defectiveScores_defectiveInts = unlist(defectiveScores$defectiveInts)

prob_interactions_active_OK__of__defectiveInts = prob_interactions_active_OK[names(unlist__defectiveScores_defectiveInts)]
sort_prob_interactions_active_OK__of__defectiveInts = sort(prob_interactions_active_OK__of__defectiveInts)

save(prob_interactions_active_OK__of__defectiveInts,file=paste(drugScores_folder_of_current_script_,"prob_interactions_active_OK__of__defectiveInts_",thisDrugable,"_",searchInactiveInts,"_",groupingMode,".RData",sep=""))                                 # MR inserted


max_prob_interactions_active_OK__of__defectiveInts = max(prob_interactions_active_OK__of__defectiveInts)
min_prob_interactions_active_OK__of__defectiveInts = min(prob_interactions_active_OK__of__defectiveInts)

diff__0_5_max_prob_interactions_active_OK__of__defectiveInts = max_prob_interactions_active_OK__of__defectiveInts - 0.5
diff__0_5_min_prob_interactions_active_OK__of__defectiveInts = 0.5 -min_prob_interactions_active_OK__of__defectiveInts

percentage_max_num_patients_in_drug_subgroups = max(num_patients_in_drug_subgroups)/nrow(allMedianNetworks_OK)

percentage_patientsHealthy = num_patientsHealthy/nrow(allMedianNetworks_OK)
































predicted_defectiveInts = defectiveScores$defectiveInts

# > num_codruggableInteractions_of_drug_subgroups
# Gilenya     IFNb Copaxone     EGCG  Tysabri 
# 13        3       14       23        4 

matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore = matrix(data = NA,nrow=sum(num_codruggableInteractions_of_drug_subgroups),ncol = 3)
matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[1:num_codruggableInteractions_of_drug_subgroups[1]] = 
  
  matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[,1] = c(rep(phenotypes[1],num_codruggableInteractions_of_drug_subgroups[1]),
                                                                                                         rep(phenotypes[2],num_codruggableInteractions_of_drug_subgroups[2]),
                                                                                                         rep(phenotypes[3],num_codruggableInteractions_of_drug_subgroups[3]),
                                                                                                         rep(phenotypes[4],num_codruggableInteractions_of_drug_subgroups[4]),
                                                                                                         rep(phenotypes[5],num_codruggableInteractions_of_drug_subgroups[5]))

matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[1:cumsum(num_codruggableInteractions_of_drug_subgroups)[1],2] = names(predicted_defectiveInts[[1]])
matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[1:cumsum(num_codruggableInteractions_of_drug_subgroups)[1],3] = unname(predicted_defectiveInts[[1]])


matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[(cumsum(num_codruggableInteractions_of_drug_subgroups)[1]+1) :(cumsum(num_codruggableInteractions_of_drug_subgroups)[2]) ,2] = names(predicted_defectiveInts[[2]])
matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[(cumsum(num_codruggableInteractions_of_drug_subgroups)[1]+1) :(cumsum(num_codruggableInteractions_of_drug_subgroups)[2]) ,3] = unname(predicted_defectiveInts[[2]])

matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[(cumsum(num_codruggableInteractions_of_drug_subgroups)[2]+1) :(cumsum(num_codruggableInteractions_of_drug_subgroups)[3]) ,2] = names(predicted_defectiveInts[[3]])
matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[(cumsum(num_codruggableInteractions_of_drug_subgroups)[2]+1) :(cumsum(num_codruggableInteractions_of_drug_subgroups)[3]) ,3] = unname(predicted_defectiveInts[[3]])

matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[(cumsum(num_codruggableInteractions_of_drug_subgroups)[3]+1) :(cumsum(num_codruggableInteractions_of_drug_subgroups)[4]) ,2] = names(predicted_defectiveInts[[4]])
matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[(cumsum(num_codruggableInteractions_of_drug_subgroups)[3]+1) :(cumsum(num_codruggableInteractions_of_drug_subgroups)[4]) ,3] = unname(predicted_defectiveInts[[4]])

matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[(cumsum(num_codruggableInteractions_of_drug_subgroups)[4]+1) :(cumsum(num_codruggableInteractions_of_drug_subgroups)[5]) ,2] = names(predicted_defectiveInts[[5]])
matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore[(cumsum(num_codruggableInteractions_of_drug_subgroups)[4]+1) :(cumsum(num_codruggableInteractions_of_drug_subgroups)[5]) ,3] = unname(predicted_defectiveInts[[5]])

colnames(matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore) = c("Drug","Interaction","DrugScore")

save(matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore,file=paste(drugScores_folder_of_current_script_,"matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore_",thisDrugable,"_",searchInactiveInts,"_",groupingMode,".RData",sep=""))                                 # MR inserted
write.table(matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore,file=paste(drugScores_folder_of_current_script_,"matrix_CoDruggableInteraction_Healthy_VS_DrugTreatedMS__COL1Drug_COL2Interaction_COL3DrugScore_",thisDrugable,"_",searchInactiveInts,"_",groupingMode,".csv"),sep=",",row.names=F)

# END MR inserted 



print("Script finished!")              # MR inserted

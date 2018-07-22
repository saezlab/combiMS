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

# option searchInactiveInts has no effect in this script as
# the searchInactiveInts="yes" related modification of
# drugNetwork$network
# by setting
# drugNetwork$network[which(names(drugNetwork$network) %in% interactionsToReplace)]=1 #replace activity by 1
# 
# is afterwards not used and not saved






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
   
   applyLinkActivityThreshold__storing_text = 'NOTroundedRealNumber_'
   
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


drug=5 #**********select a drug here
#Gilenya=Fingolimod, Glatiramer=Copaxone
phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
numDrugs=5
# #plot(allDrugNws[[drug]]$graph)


# *************************************************************************************************************************
# *********** 2. Select only active parts of nw above activity threshold, then turn into graph
# *************************************************************************************************************************
# Transform network into graph. Booleanize network if mode is mean
drugNetwork=network2graph(allDrugNws[[drug]]$network,
                          mode=groupingMode,
                          linkActivityThreshold=linkActivityThreshold)  # Note: The function network2graph requires a Boolean logic network (using not rounded mean causes several wrong results). 
                                                                        # Thus the option applyLinkActivityThreshold must be applyLinkActivityThreshold = 'no' and can not be chosen by the user (to prevent wrong results).     


plotModel(drugNetwork$model,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))

#Note the difference between activity in network (grouping), co-drugability score, and "situation not ok"score, i.e. difference between grouping in H and D.
#calculateDefective can work in "zero" scoring mode, identifying some interactions with a socre of 0 where situation is not ok. 
# network2graph will cut such interactions out, thereby if we search that graph such interactions cannot be used. A rationale has been
# discussed where it is interesting to use them. Therefore, before prediction of combination therapy, force activation of these interactions.
# note:intereactions with 0 scores and situation ok have been replaced by 1 in their score, hence it is safe to just activate all remaining ints with 0 scores
if(searchInactiveInts=="no"){
   
  warning(paste0("searching inactive ints: ",searchInactiveInts))
   
} else if (searchInactiveInts=="yes"){
   
  phenotypeScores=read.table(file=paste(drugScores_folder,phenotypes[drug],"__",applyLinkActivityThreshold__storing_text,groupingMode,"_",thisDrugable,".csv",sep=""),sep=",")
  
   
  interactions0Score=drugNetwork$network[which(names(drugNetwork$network) %in% rownames(phenotypeScores)[which(phenotypeScores==0)])]
  interactionsToReplace=names(which(interactions0Score==0))
  drugNetwork$network[which(names(drugNetwork$network) %in% interactionsToReplace)]=1 #replace activity by 1
  warning(paste0(interactionsToReplace," ints replaced with activity=1. Total of ",length(interactions0Score)," ints with 0 score and situation not ok"))
} else {stop(warning("incorrect searching mode selected"))}

# *************************************************************************************************************************
# ***********3. Once a drug is selected, recreate drug network, then
#************determine if a stimulus connects to signals through Defective interactions
# *************************************************************************************************************************
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

cat("Defective intearctions:",topTargets,"\n")
cat("Present defective:",presentTopInts,"\n")
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
experimentsThisInteraction[,4]=rep(defectiveScores$defectiveInts[[drug]][indexInts],numCombinations) 



# experimentsThisInteraction=data.frame(Reaction=rep(link,numCombinations),
#                                       Stimulus=Stimuli[which(matrixPaths==1,arr.ind=T)[,2]],
#                                       Readout=AllReadouts[which(matrixPaths==1,arr.ind=T)[,1]],
#                                       drugScoreInteraction=rep(defectiveScores$defectiveInts[[drug]][indexInts],numCombinations) 
#                             )
combinationExperiments=rbind(combinationExperiments,experimentsThisInteraction)
#aheatmap(matrixPaths,Rowv=NA,Colv=NA)
}

drugName=phenotypes[drug]
write.csv(combinationExperiments,file=paste("/Users/marti/Desktop/figuresCombiMS/combinations/",drugName,thisDrugable,searchInactiveInts,groupingMode,".csv",sep=""))
write.table(drugNetwork$cutNetwork,file=paste0(phenotypeNws_folder,drugName,"CUT",thisDrugable,searchInactiveInts,groupingMode,".csv"),sep=",",row.names=T,quote=F)

cat('**************num combinations by stimuli',length(unique(combinationExperiments[,2])),'\n')
cat('**************num combinations by reaction',length(unique(combinationExperiments[,1])),'\n')
# *************************************************************************************************************************
# *********** test each group against each other - careful with the interpretation, as it should be one group against all others
# *************************************************************************************************************************
#see pathDrugTargetsFinal.R


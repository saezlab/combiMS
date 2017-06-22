#***********************************************************************
# This scripts loads the results of loadAllModelsPlaneNW2.R,
# which is the concatenation of all networks solutions within 0.0000005% reltol for 
# the runs that did not crush out of 10 attempts per patient
# Juny 2015
# Marti Bernardo-Faura

# Minor changes to adjust the paths to the combiMS Github project
# Jakob Wirbel, June 2017
#***********************************************************************

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#library(reshape2)
#library(ggplot2)
library(CellNOptR)
#library(corrplot)
#library(vegan)
source("./calculateScore.R")
source("./edgeContribution2.R") # is the 2 intentional?
source("./loadPatientRuns.R")
#******************************************************************
#*************** define directory depending on appropriate cluster run
#***********************************************************************
run_folder="../../files/cluster/7th10hReltol005/results/"
networks_folder="../../files/cluster/7th10hReltol005/networks/"

similarity_folder="../../files/cluster/7th10hReltol005/similarity/"
completed_patients_folder="../../files/cluster/completedPatients/"
networks_for_Completion="../../files/cluster/8th10hReltol005/results/"
modelsForCompletion=list.files(networks_for_Completion,pattern="*.RData",full.names=FALSE)


solutionsOneRunPatient=list.files(run_folder,pattern="*.RData",full.names=FALSE) 
numModels=length(solutionsOneRunPatient)
solutionsAllRunsPatient=list.files(networks_folder,pattern="*.RData",full.names=FALSE)
#***********************************************************************
# *************** preprocess model
#***********************************************************************
data_folder="/Users/marti/Documents/ebi/combiMS/data/phosphosMergedAbsMax/processed/normalized/secondRoundProcessedMidas/"
patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
model_path="/Users/marti/Documents/R/combiMS/combiMSplane.sif"
fileName=patientData[1]

midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)  
sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))

model=preprocessing(midas,model,expansion=FALSE)
numInteractions=length(model$reacID)
sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))


#***********************************************************************
# *************** identify patients for which some solutions crashed, which need to be completed
#***********************************************************************
numRuns=vector()
for ( i in 1:length(patientData)){
  load(paste0(networks_folder,solutionsAllRunsPatient[i])) # this loads PatientResults
  cat("--patient",PatientResults$name,"\n")
  numRuns[i]=length(PatientResults$NwsInRuns)
}
numPatientsToComplete=length(numRuns[which(numRuns<10)])
cat("number of patients that need to be completed",numPatientsToComplete)
AllPatientsToComplete=vector(length=numPatientsToComplete)
AllPatientsToComplete=10-numRuns[which(numRuns<10)]
names(AllPatientsToComplete)=solutionsAllRunsPatient[which(numRuns<10)]

save(AllPatientsToComplete,file=paste0(completed_patients_folder,"runsToCompleteAllPatients",".RData"))
#***********************************************************************
# *************** generate a merged structure with the names of all files for each patient that will be used for completion
#***********************************************************************
# *************** not needed anymore, this is now in loadPatientRuns.R

#***********************************************************************
# *************** load data to complete patients
#***********************************************************************
networks_for_Completion="/Users/marti/Documents/R/combiMS/cluster/8th10hReltol005/results/"
for(i in 1:numPatientsToComplete){
  patientToComplete=names(AllPatientsToComplete[i])
  runsNeededPatient=AllPatientsToComplete[i]  
  networksToAdd=loadPatientRuns(patientToComplete,runsNeededPatient,numInteractions,networks_for_Completion)
  save(networksToAdd,file=paste0(completed_patients_folder,patientToComplete))
}


#***********************************************************************
# *************** bind data to complete patients from 7th run with the above from 8th
#***********************************************************************
# double check number of runs to complete
networks_folder="/Users/marti/Documents/R/combiMS/cluster/7th10hReltol005/networks/"
totalPatients=list.files(networks_folder,pattern="*.RData",full.names=FALSE)
incompletePatientSummary=list()
indexIncomplete=0
for(indexBoth in 1:length(totalPatients)){
  
  load(paste0(networks_folder,totalPatients[indexBoth]))
  cat("----------completing patient",PatientResults$name,"\n")
  if(length(PatientResults$NwsInRuns)==10){
    cat("patient complete","\n")
  } else {
    indexIncomplete=indexIncomplete+1
    incompletePatientSummary[[indexIncomplete]]<-list()
    incompletePatientSummary[[indexIncomplete]]$runsCompleted=length(PatientResults$NwsInRuns)
    incompletePatientSummary[[indexIncomplete]]$name=PatientResults$name
    cat("incomplete patient",PatientResults$name,"total incomplete",indexIncomplete,"\n")
  }
}
# bind to networks previously loaded from 8th
completed_patients_folder="/Users/marti/Documents/R/combiMS/cluster/completedPatients/"
patientsToBind=list.files(completed_patients_folder,pattern="*.RData",full.names=FALSE)
patientsToBind=patientsToBind[-c(10,11)] #remove two directories

for (i in 1:length(incompletePatientSummary)){
  cat(paste0(incompletePatientSummary[[i]]$name,".RData") %in% patientsToBind,incompletePatientSummary[[i]]$name,"\n")
}
#so...
incomplete_patients=c('UZ022.RData', 'UZ027.RData', 'IB068.RData', 'IB030.RData')
# now bind the others, which have been prepared previously
for (i in 1:length(patientsToBind)){
  load(paste0(completed_patients_folder,patientsToBind[i]))#this loads networksToAdd
  
  networks_folder="/Users/marti/Documents/R/combiMS/cluster/7th10hReltol005/networks/"
  load(paste0(networks_folder,patientsToBind[i]))#this loads PatientResults
  cat("completing",PatientResults$name,"with",networksToAdd$name,"\n")
  cat("completing",length(PatientResults$NwsInRuns),"with",length(networksToAdd$NwsInRuns),"\n")
  completed=list()
  completed$name=PatientResults$name
  completed$AllNwsPatient=rbind(PatientResults$AllNwsPatient,networksToAdd$AllNwsPatient)
  completed$AllScores=c(PatientResults$AllScores,networksToAdd$AllScores)
  completed$NwsInRuns=c(PatientResults$NwsInRuns,networksToAdd$NwsInRuns)
  save(completed,file=paste0(completed_patients_folder,"completed/",patientsToBind[i]))
  
}
#***********************************************************************
# *************** 15 patients were not complete in 7th cluster run.
#After completing above using 8th run, 4 are still not complete.
#Here we finsih them by searching in 9th cluster run
#***********************************************************************
incomplete_patients=c('UZ022.RData', 'IB068.RData', 'IB030.RData', 'UZ027.RData')
i=4
#load PatientResults from 7th run
load(paste0(networks_folder,incomplete_patients[i]))

#load networksToAdd from 8th run
load(paste0(completed_patients_folder,"Lacking2runs/",incomplete_patients[i]))

# print to be sure
cat("completing",PatientResults$name,"with",networksToAdd$name,"\n")
cat("completing",length(PatientResults$NwsInRuns),"with",length(networksToAdd$NwsInRuns),"\n")

#..and in the darkness bind them
re_completed=list()
re_completed$name=PatientResults$name
re_completed$AllNwsPatient=rbind(PatientResults$AllNwsPatient,networksToAdd$AllNwsPatient)
re_completed$AllScores=c(PatientResults$AllScores,networksToAdd$AllScores)
re_completed$NwsInRuns=c(PatientResults$NwsInRuns,networksToAdd$NwsInRuns)

#calculate how many runs are still missing
runsNeeded=10-length(re_completed$NwsInRuns)

#load missing runs from 9th run
networks_for_Re_Completion="/Users/marti/Documents/R/combiMS/cluster/9th10hReltol005/results/"
NetworksToAdd=loadPatientRuns(incomplete_patients[i],runsNeeded,numInteractions,networks_for_Re_Completion)

# finally bind from 9th (in thrice_completed) to result of binding from 7th (in networksToAdd) and 8th (in re_completed)
thrice_completed=list()
thrice_completed$name=NetworksToAdd$name
thrice_completed$AllNwsPatient=rbind(re_completed$AllNwsPatient,NetworksToAdd$AllNwsPatients)
thrice_completed$AllScores=c(re_completed$AllScores,NetworksToAdd$AllScores)
thrice_completed$NwsInRuns=c(re_completed$NwsInRuns,NetworksToAdd$NwsInRuns)
save(thrice_completed,file=paste0(completed_patients_folder,"completed/",incomplete_patients[i]))
#save(thrice_completed,file=paste0(completed_patients_folder,"stillLacking1RunAfter9th/",incomplete_patients[i]))

#***********************************************************************
# *************** IB068 is still incomplete by one run
#***********************************************************************
incomplete_patient='IB068.RData'
networks_folder="/Users/marti/Documents/R/combiMS/cluster/completedPatients/stillLacking1RunAfter9th/"

# load networks completed from 7th, 8th, and 9th run
load(paste0(networks_folder,incomplete_patient))

#load networksToAdd from 9bis run
networks_for_Completion="/Users/marti/Documents/R/combiMS/cluster/completedPatients/OnePatient/"
runsNeededPatient=10-length(thrice_completed$NwsInRuns) 
NetworksToAdd=loadPatientRuns(incomplete_patient,runsNeededPatient,numInteractions,networks_for_Completion)

# bind them
completed_patients_folder="/Users/marti/Documents/R/combiMS/cluster/completedPatients/completed/"

four_completed=list()
four_completed$name=NetworksToAdd$name
four_completed$AllNwsPatient=rbind(thrice_completed$AllNwsPatient,NetworksToAdd$AllNwsPatients)
four_completed$AllScores=c(thrice_completed$AllScores,NetworksToAdd$AllScores)
four_completed$NwsInRuns=c(thrice_completed$NwsInRuns,NetworksToAdd$NwsInRuns)
save(four_completed,file=paste0(completed_patients_folder,incomplete_patient))





#***********************************************************************
# *************** clean network solutions
#***********************************************************************
i=1
load(paste0(networks_folder,solutionsAllRunsPatient[i])) # this loads PatientResults

# calculate networks within 1% rel tol to shrink the original 5%
newRelTol=0.001
scoreTol=min(PatientResults$AllScores)+min(PatientResults$AllScores)*newRelTol
ModelsInTol=PatientResults$AllNwsPatient[which(PatientResults$AllScores<scoreTol),]
distanceModels=newJaccard(ModelsInTol,ModelsInTol)
#subset to shrink even more and test dissimilarity
subsetModels=ModelsInTol[1:10,] 
#calculate dissimilarity across runs in single patient

distanceModels=as.matrix(vegdist(subsetModels, method="jaccard"))
write.table(distanceModels,file=paste0(similarity_folder,patientData[i]),sep=",")

corrplot.mixed(distanceModels,is.corr=F)

# 
# #prune networks to remove spurious interactions
# prova=edgeContribution2(ModelsInTol[i,],model,midas)
# # are all networks unique? this seems to indicate too large solution space, as they should duplicate through runs
# i=1
# dim(ModelsInTol)[1]==dim(unique(ModelsInTol))[1]
# 
# #distanceModelsInTol=vegdist(ModelsInTol, method="jaccard")
# 
# #transform intor string, find out which ones are the same
# stringSolutions=apply(PatientResults$AllNwsPatient, 1, function(x) paste(x, collapse=","))
# 

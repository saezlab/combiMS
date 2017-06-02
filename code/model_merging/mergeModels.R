# Function to merge all cluster run solutions to allow each single model to be calculated from successful 10 cluster runs.
# If less than 10 jobs converged, the job was submitted again. When 10  solutions were obtained, all networks found within 0.00005 relative tolerance of the best networkd were concatenated.
# This allowed fair merging of all networks for subsequent single model calculation and was analysed in depth (other functions used)
# I investigated the reason for non-convergence, which was due to memory usage linked to larger number of networks withi rel tol (see code for model quality control)
# Created by Marti Bernardo-Faura, final version July 2015

#library(reshape2)
#library(ggplot2)
library(CellNOptR)

# *************** set working directory for relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


source("../utils/calculateScore.R")
# source("/Users/marti/Documents/r/combiMS/edgeContribution2.R")
# source("/Users/marti/Documents/r/combiMS/loadPatientRuns.R")

#***********************************************************************
# *************** preprocess model
#***********************************************************************
data_folder="../../data/phosphos_processed/"
patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
model_path="../../files/model/combiMSplaneCUT.sif"
fileName=patientData[1]

midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)  
sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))

# model=preprocessing(midas,model,expansion=FALSE)
numInteractions=length(model$reacID)
# sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))

#***********************************************************************
# *************** test that score is calculated properly
#***********************************************************************
#keptEdges=which(bString==1)
# initBstring=rep(1,times=numInteractions)
# Opt=gaBinaryT1(CNOlist=midas,model=model,initBstring=initBstring,stallGenMax=600,maxTime=20, maxGens=100000, verbose=FALSE,popSize=100,elitism=2)
# myscore$score=calculateScore(model,midas,bString=Opt$bString)
# myscore$score==Opt$bScore



#***********************************************************************
# *************** # calculate models in new tol. calculate median nw. 
# save as two structures:
#       create structures with  (i) median network and (ii) its score for subsequent network analysis
#       create structures with  (iii) mean network -score can not be calculated on weighted nw
#       create structure with (iv) stats for sanity check, i.e. are num networks in runs and error related
#***********************************************************************

networks_folder="../../files/all_solutions_models/"
totalPatients=list.files(networks_folder,pattern="*.RData",full.names=FALSE)
totalPatientsSimple=sapply(totalPatients, function(x){
  strsplit(x,"\\.")[[1]][1]
})

num_patients=length(totalPatients)
#num_patients=3
allMedianNetworks=matrix(nrow=num_patients,ncol=numInteractions)                         
medianNetworkErrors=list()
allMeanNetworks=matrix(nrow=num_patients,ncol=numInteractions)                         
total_patientDF=data.frame(name=character(),num_nws=vector(mode="numeric"),avg_score=vector(mode="numeric"),total_nws=vector(),best_score=vector(),merged_model_score=vector(),nws_new_reltol=vector(),convergence=character())


# convergence_2=c("CH024.RData", "CH055.RData", "IB027.RData", "IB041.RData", "IB044.RData", "IB054.RData", "IB066.RData", "IB067.RData", "IB069.RData", "UZ002.RData","UZ004.RData")
# convergence_3=c('UZ022.RData', 'IB030.RData', 'UZ027.RData')
# convergence_4='IB068.RData'

# create annotation vector for convergence
conv=rep("1",169)
convergence_2=c("CH024", "CH055", "IB027", "IB041", "IB044", "IB054", "IB066", "IB067", "IB069", "UZ002","UZ004")
convergence_3=c("UZ022", "IB030", "UZ027")
convergence_4="IB068"
names(conv)=totalPatientsSimple
conv[convergence_2]="2"
conv[convergence_3]="3"
conv[convergence_4]="4"

#for(index in 1:length(totalPatients)){
for(index in 1:num_patients){
  
  cat("----------loading patient",totalPatients[index],"\n")
  load(paste0(networks_folder,totalPatients[index])) # this loads PatientResults
  
  # this is not necessary since i use completePatientsFinal.R
  # if patient was completed the name of the object needs to be changed
#   if (totalPatientsSimple[index] %in% convergence_2){
#     cat("patient completed twice","\n")
#     PatientResults=completed
#   }else if (totalPatientsSimple[index] %in% convergence_3){
#     cat("patient completed thrice","\n")
#     PatientResults=thrice_completed
#   }else if (totalPatientsSimple[index] %in% convergence_4){
#     cat("patient completed four x","\n")
#     PatientResults=four_completed
#   }
#   
  # check that the num networks matches the scores
  if (! dim(PatientResults$AllNwsPatient)[1] == length(PatientResults$AllScores)){
    cat("***************************************************************************","\n")
    cat("dim networks does not macth num scores!!","\n")
  }
  
  
  #extract info from this patient
  num_success_runs=length(PatientResults$NwsInRuns)
  
  if (num_success_runs<10){
    cat("patient incomplete",PatientResults$name,"!!","\n")
  } else {
    
    # because networks come from different runs, new rel tol must be calculated
    newRelTol=0.00005 # 5 per hundred thowsend
    cat("recalculating models in global tolerance","\n")
    scoreTol=min(PatientResults$AllScores)+min(PatientResults$AllScores)*newRelTol
    ModelsInTol=PatientResults$AllNwsPatient[which(PatientResults$AllScores<scoreTol),]
    
    # clean repeated networks to keep only unique
    unique_nws=unique(ModelsInTol)
    if(!dim(ModelsInTol)[1]==dim(unique_nws)[1]){
      cat(dim(unique_nws)[1],"were unique out of", dim(ModelsInTol)[1],"\n")
      ModelsInTol=unique_nws
    }
    nws_new_reltol=dim(ModelsInTol)[1]
    
    # create 1st structure: 1 median network 
    cat("calculating median of interactions","\n")
    medianNw=apply(ModelsInTol,2,median)
    allMedianNetworks[index,]=medianNw
    
    # 2nd structure: error for medians
    #midas MUST be that of right patient
    fileName=patientData[index]
    midas=CNOlist(paste(data_folder,fileName,sep=""))
    cat("calculating score of median network, data is",fileName,"\n")    
    medianNetworkErrors[[index]]=calculateScore(model,midas,bString=medianNw)
    
    # create 3rd structure: 1 mean network 
    cat("calculating average of interactions","\n")
    meanNw=apply(ModelsInTol,2,mean)
    allMeanNetworks[index,]=meanNw

#     # 4th structure:error for mean ns
#     cat("calculating score of median network","\n")
#     meanNetworkErrors[[index]]=calculateScore(model,midas,bString=meanNw)
    
    # create 4th structure: stats for sanity check, e.g. are num networks (in each run) and error related?
    cat("calculating number of networks per run","\n")
    one_patientDF=data.frame(name=character(length=10),num_nws=vector(mode="numeric",length=10),avg_score=vector(mode="numeric",length=10),total_nws=vector(mode="numeric",length=10),best_score=vector(mode="numeric",length=10),merged_model_score=vector(mode="numeric",length=10),nws_new_reltol=vector(mode="numeric",length=10),convergence=character(length=10))
    one_patientDF$name=rep(PatientResults$name,num_success_runs)
    one_patientDF$num_nws=PatientResults$NwsInRuns
    one_patientDF$convergence=rep(conv[PatientResults$name],10)
    one_patientDF$total_nws=sum(PatientResults$NwsInRuns)
    one_patientDF$merged_model_score=rep(medianNetworkErrors[[index]]$score,10)
    one_patientDF$nws_new_reltol=rep(nws_new_reltol,10)

    #for all networks in one run calculate average error and best error
    run_start_index=1
    run_end_index=PatientResults$NwsInRuns[1]
    i=1
    cat("calculating avg score for each run","\n")
    cat("calculating best score for each run","\n")

    while (i <=num_success_runs){
      scores_one_run=PatientResults$AllScores[run_start_index:run_end_index]
      one_patientDF$avg_score[i]=mean(scores_one_run,na.rm=T)
      one_patientDF$best_score[i]=min(scores_one_run,na.rm=T)
      
      #move indexes
      run_start_index=run_start_index+PatientResults$NwsInRuns[i]
      run_end_index=run_end_index+PatientResults$NwsInRuns[i]
      i=i+1
      
    }
    total_patientDF=rbind(total_patientDF,one_patientDF) #concatenate this with total
    
  }
}

if(dim(total_patientDF)[1]!=1690){stop("Hey!! Not all patients were added, is IB068 missing?")}

#***********************************************************************
# *************** save all models and metadata structures
#***********************************************************************
rownames(allMedianNetworks)=totalPatientsSimple
colnames(allMedianNetworks)=model$reacID
save(allMedianNetworks,file="../../files/median_models/allMedianModels.RData")

rownames(allMeanNetworks)=totalPatientsSimple
colnames(allMeanNetworks)=model$reacID
save(allMeanNetworks,file="../../files/median_models/allMeanModels.RData")

names(medianNetworkErrors)=totalPatientsSimple
save(medianNetworkErrors,file="../../files/median_models/allErrors.RData")

# *************** calculate best score across runs
patient_names=unique(total_patientDF$name)
for (i in 1:length(patient_names)){
  all_best_scores=vector(mode="numeric",length=10)
  all_best_scores=total_patientDF$best_score[which(total_patientDF$name==patient_names[i])]
  total_patientDF$best_of_bests[which(total_patientDF$name==patient_names[i])]=rep(min(all_best_scores),10)
}
# ************** add an identifier for runs
run_id=c("1","2","3","4","5","6","7","8","9","10")
for (i in 1:length(patient_names)){
  total_patientDF$run_id[which(total_patientDF$name==patient_names[i])]=run_id
}
save(total_patientDF,file="../../files/model_merging/statsModels.RData")

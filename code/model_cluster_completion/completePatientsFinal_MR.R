# Function to merge all cluster run solutions to allow each single model to be calculated from successful 10 cluster runs.
# If less than 10 jobs converged, the job was submitted again. When 10  solutions were obtained, all networks found within 0.00005 relative tolerance of the best networkd were concatenated.
# This allowed fair merging of all networks for subsequent single model calculation and was analysed in depth (other functions used)
# I investigated the reason for non-convergence, which was due to memory usage linked to larger number of networks withi rel tol (see code for model quality control)
# Created by Marti Bernardo-Faura, final version July 2015

# Minor changes to adjust the paths to the combiMS Github project
# Jakob Wirbel, June 2017


# with some notes added by 
# Melanie Rinas
# November 2017


# One small "error " in the original script detected:
# numInteractions=length(model$reacID)     needs to be uncommented as needed 


# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(CellNOptR)
source("./loadPatientRuns.R")
source("./calculateScore.R")
source("./bindTwoRuns.R")
source("./checkNwsAndScores.R")
#******************************************************************
#*************** define directories with optimisation results
#***********************************************************************
first_folder="../../files/cluster/7th10hReltol005/results/"
second_folder="../../files/cluster/8th10hReltol005/results/"
third_folder="../../files/cluster/9th10hReltol005/results/"
fourth_folder="../../files/cluster/9thB/OnePatient/"
all_folders=c(first_folder,second_folder,third_folder,fourth_folder)
#***********************************************************************
# *************** preprocess model
#***********************************************************************
data_folder="../../data/phosphos_processed/"
patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
model_path="../../files/model/combiMSplaneCUT.sif"
fileName=patientData[grep("IB068",patientData)]

midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)  
sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))

# model=preprocessing(midas,model,expansion=FALSE)
numInteractions=length(model$reacID)                    # MR modified
# sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))

# initBstring=rep(1,length(model$reacID))
# Opt=gaBinaryT1(CNOlist=midas,model=model,initBstring=initBstring,stallGenMax=600,maxTime=20, maxGens=100000, verbose=FALSE,popSize=100,elitism=2)
# prova=calculateScore(model,midas,bString=initBstring)

#***********************************************************************
# *************** generate lists of patients that need to be completed
#***********************************************************************
load("../../files/completedPatients/runsToCompleteAllPatients.RData")
#this loads AllPatientsToComplete, which was calculated (TBC) in studyModelDistribution.R

#***********************************************************************
# *************** complete
#***********************************************************************

total_pat_complete=length(AllPatientsToComplete)
problem_finding_files=vector(mode="character")
for (i in 1:total_pat_complete){
  patientToComplete=names(AllPatientsToComplete)[i]
  runs=1
  # load midas for error calculation
  midasName=strsplit(patientToComplete,"\\.")[[1]][1]
  midasName=paste0(midasName,".csv")
  cat("****************loading midas for patient",midasName,"\n")
  midas=CNOlist(paste(data_folder,midasName,sep=""))
  
  # load data from 7th run
  runsNeededPatient=10
  networks_folder=all_folders[runs]
  cat("-------------------------------------------------loading patient",patientToComplete,"***run ",runs,"\n")
  networks_first_run=loadPatientRuns(patientToComplete,runsNeededPatient,numInteractions,networks_folder)
  
  if(is.na(networks_first_run$AllNwsPatient[1,1])){
    cat("there is a weird NA line, cause unknown. Remove!","\n")
    networks_first_run$AllNwsPatient=networks_first_run$AllNwsPatient[-1,]
  }
  
  #Check if removing first NA line messed up
  recalculated_score=vector()
  for (j in 1:2){
    recalculated_score[j]=calculateScore(model,midas,bString=networks_first_run$AllNwsPatient[j,])$score
  }
  checkNwsAndScores(networks_first_run,recalculated_score)
  
  # load data from as many runs as necessary, max is 4
  runsNeededPatient=runsNeededPatient-length(networks_first_run$NwsInRuns)
  PatientResults=networks_first_run
  awesomness=F
  runs=2
  while (awesomness==F){
    
    networks_folder=all_folders[runs]
    cat("*****************loading patient",patientToComplete,"***run ",runs,"\n")
    networks_additional=loadPatientRuns(patient_name=patientToComplete,
                                        runsNeededPatient=runsNeededPatient,
                                        numInteractions=numInteractions,
                                        networks_for_Re_Completion=networks_folder)
    #there is an indentified bug in the code: 
    #for IB030 and IB068, loadPatientRuns finds the folder of results in 3rd runs empty, hence returns empty networks.
    #however, if repeated it is then fine
    
    if(dim(networks_additional$AllNwsPatient)[1]==0){
      networks_additional=loadPatientRuns(patient_name=patientToComplete,
                                          runsNeededPatient=runsNeededPatient,
                                          numInteractions=numInteractions,
                                          networks_for_Re_Completion=networks_folder)
      problem_finding_files=c(problem_finding_files,patientToComplete)
    }
    
    
    if(is.na(networks_additional$AllNwsPatient[1,1])){
      cat("there is a weird NA line, cause unknown. Remove!","\n")
      networks_additional$AllNwsPatient=networks_additional$AllNwsPatient[-1,]
    }
    
    #Check if removing first NA line messed up 
    recalculated_score=vector()
    for (j in 1:2){
      recalculated_score[j]=calculateScore(model,midas,bString=networks_additional$AllNwsPatient[j,])$score
    }
    checkNwsAndScores(networks_additional,recalculated_score)
    
    # *************** bind results from two rounds
    cat("****binding results from first run and run",runs,"\n")
    PatientResults=bindTwoRuns(PatientResults,networks_additional)
    
    # update vars to check if we are done
    runsNeededPatient=runsNeededPatient-length(networks_additional$NwsInRuns)
    runs=runs+1
    if(runsNeededPatient==0){awesomness=T}    
  }
  
  # *************** save
  cat("saving patient","\n")
  save(PatientResults,file=paste0("../files/completedPatients/completed/",patientToComplete))
  cat("-------------------------------------------------patient completed","\n")
  
}  

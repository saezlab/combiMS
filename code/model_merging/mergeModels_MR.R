
# Function to create a single model from all solutions found within relative tolerance of the best model for a given patient
# (-) For a given patient, take a subset of models within reltol. Then clean the non unique.
# (-) For each patient, calculate a single network (i.e. the mean and the median), and keep stats such as error, num of networks within reltol
# (-) Merge all single models (each one for a patient) into one structure of networks
# Created by Marti Bernardo-Faura, final version July 2015


# Some modifications regarding avg_score and DonorBestNetworks calculation and storing
# Adaption of paths, SIF file, initializations of objects
# by Melanie Rinas, November 2017



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
model_path="../../files/combiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED.sif"

fileName=patientData[1]

midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)  
#sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
sprintf("********* The used model has  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))                                       


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





sub_dir__Results_model_merging = "../../files/model_merging"

sub_dir__Results_median_models = "../../files/median_models"




networks_folder="../../files/all_solutions_models/"
totalPatients=list.files(networks_folder,pattern="*.RData",full.names=FALSE)
totalPatientsSimple=sapply(totalPatients, function(x){
  strsplit(x,"\\.")[[1]][1]
})

num_patients=length(totalPatients)
#num_patients=3
allMedianNetworks=matrix(nrow=num_patients,ncol=numInteractions)          
DonorBestNetworks=list()           

medianNetworkErrors=list()

allMeanNetworks=matrix(nrow=num_patients,ncol=numInteractions)                         
total_patientDF=data.frame(name=character(),
                           num_nws=vector(mode="numeric"),
                           avg_score=vector(mode="numeric"),
                           total_nws=vector(),
                           best_score=vector(),
                           merged_model_score=vector(),
                           nws_new_reltol=vector(),
                           convergence=character())




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
  # load(paste0(networks_folder,totalPatients[index])) # this loads PatientResults
  load(paste(networks_folder,totalPatients[index],sep="/")) # this loads PatientResults        
  
  
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
    # 
    # 5/100 = 0.05 = 5 %
    # 5/100000 = 5e-05 = 0.00005
    # 
    newRelTol=0.00005 # 5 per hundred thowsend
    cat("recalculating models in global tolerance","\n")
    scoreTol=min(PatientResults$AllScores)+min(PatientResults$AllScores)*newRelTol
    ModelsInTol=PatientResults$AllNwsPatient[which(PatientResults$AllScores<scoreTol),]
    
    

    if(length(which(PatientResults$AllScores<scoreTol)) == 1){
      
      ModelsInTol = as.matrix(t(ModelsInTol))              # tranform vector to matrix 
      
    }
    
    
    
   
    
    unique_nws=unique(ModelsInTol)
    if(!dim(ModelsInTol)[1]==dim(unique_nws)[1]){           # command is the same as: if(!(dim(ModelsInTol)[1]==dim(unique_nws)[1])){
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
    # medianNetworkErrors[[index]]=calculateScore(model,midas,bString=medianNw)      
    medianNetworkErrors[[index]]=calculateScore(model,midas,bString=medianNw)     
    
    # create 3rd structure: 1 mean network 
    cat("calculating average of interactions","\n")
    meanNw=apply(ModelsInTol,2,mean)
    allMeanNetworks[index,]=meanNw
    
    #     # 4th structure:error for mean ns
    #cat("calculating score of median network","\n")

    # create 4th structure: stats for sanity check, e.g. are num networks (in each run) and error related?
    cat("calculating number of networks per run","\n")
    one_patientDF=data.frame(name=character(length=10),
                             num_nws=vector(mode="numeric",length=10),
                             avg_score=vector(mode="numeric",length=10),
                             total_nws=vector(mode="numeric",length=10),
                             best_score=vector(mode="numeric",length=10),
                             merged_model_score=vector(mode="numeric",length=10),
                             nws_new_reltol=vector(mode="numeric",length=10),
                             convergence=character(length=10))
    
    
    one_patientDF$name=rep(PatientResults$name,num_success_runs)       # 1 of one_patientDF
    one_patientDF$num_nws=PatientResults$NwsInRuns                     # 2 of one_patientDF
    one_patientDF$convergence=rep(conv[PatientResults$name],10)        # 8 of one_patientDF
    one_patientDF$total_nws=sum(PatientResults$NwsInRuns)              # 4 of one_patientDF
    one_patientDF$merged_model_score=rep(medianNetworkErrors[[index]]$score,10)  # 6 of one_patientDF
    one_patientDF$nws_new_reltol=rep(nws_new_reltol,10)                          # 7 of one_patientDF
    
    

    
    

    
    
    
    #for all networks in one run calculate average error and best error
    
    
    cat("calculating avg score for each run","\n")
    cat("calculating best score for each run","\n")
    
    
    run_end_index=0
    i=1
    
    while (i <=num_success_runs){
      
      
      run_start_index= run_end_index + 1
      run_end_index= run_end_index + PatientResults$NwsInRuns[i]
      
      scores_one_run=PatientResults$AllScores[run_start_index:run_end_index]
      one_patientDF$avg_score[i]=mean(scores_one_run,na.rm=T)          # 3 of one_patientDF
      one_patientDF$best_score[i]=min(scores_one_run,na.rm=T)          # 5 of one_patientDF
      
      i=i+1
      
    }
    
    
    
    total_patientDF=rbind(total_patientDF,one_patientDF) #concatenate this with total

  }
  
 
  
  
  # create 5th structure: best networks (best of best with min score)
  
  
  if(length(which(PatientResults$AllScores==min(PatientResults$AllScores)))==1){
    
    DonorBestNetworks[[index]] = PatientResults$AllNwsPatient[which(PatientResults$AllScores==min(PatientResults$AllScores)),]   
    
  }else{
    DonorBestNetworks[[index]] = unique(PatientResults$AllNwsPatient[which(PatientResults$AllScores==min(PatientResults$AllScores)),])    
  }
  
  
  
  
 
  
  
}

if(dim(total_patientDF)[1]!=1690){stop("Hey!! Not all patients were added, is IB068 missing?")}




  #***********************************************************************
  # *************** save all models and metadata structures
  #***********************************************************************
  rownames(allMedianNetworks)=totalPatientsSimple
  colnames(allMedianNetworks)=model$reacID
  #save(allMedianNetworks,file="../../files/median_models/allMedianModels.RData")
  save(allMedianNetworks,file=file.path(sub_dir__Results_median_models,"allMedianModels.RData"))            
  
  
  rownames(allMeanNetworks)=totalPatientsSimple
  colnames(allMeanNetworks)=model$reacID
  #save(allMeanNetworks,file="../../files/median_models/allMeanModels.RData")
  save(allMeanNetworks,file=file.path(sub_dir__Results_median_models,"allMeanModels.RData"))          
  
  
  names(medianNetworkErrors)=totalPatientsSimple
  #save(medianNetworkErrors,file="../../files/median_models/allErrors.RData")
  save(medianNetworkErrors,file=file.path(sub_dir__Results_median_models,"allErrors.RData"))         
  
  
  
  names(DonorBestNetworks)=totalPatientsSimple
  #save(DonorBestNetworks,file="../../files/median_models/DonorBestNetworks.RData")
  save(DonorBestNetworks,file=file.path(sub_dir__Results_median_models,"DonorBestNetworks.RData"))          
  
  
  
  
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
  
  
  
  
  
 
  
  
  
  
  # save(total_patientDF,file="../../files/model_merging/statsModels.RData")
  save(total_patientDF,file=file.path(sub_dir__Results_model_merging,"statsModels.RData"))           

  
  
  


print("Script finished!")   


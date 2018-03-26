

# mergeModels_MR


Results_storage_name = "OptCombiMSclusterFork10_MR__using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__postCellNOptRupdate090318"  # MR inserted

save_results_of_current_script=T

# with some notes
# AND
# important modifications 
# AND
# ERROR correction
#                  regarding the indices shifts for the calculation of
#                           one_patientDF$avg_score[i]=mean(scores_one_run,na.rm=T)          # 3 of one_patientDF
#                           one_patientDF$best_score[i]=min(scores_one_run,na.rm=T)          # 5 of one_patientDF
# AND
# DonorBestNetworks calculation and storing


# 
# 
# by 
# Melanie Rinas
# November 2017



# Function to create a single model from all solutions found within relative tolerance of the best model for a given patient
# (-) For a given patient, take a subset of models within reltol. Then clean the non unique.
# (-) For each patient, calculate a single network (i.e. the mean and the median), and keep stats such as error, num of networks within reltol
# (-) Merge all single models (each one for a patient) into one structure of networks
# Created by Marti Bernardo-Faura, final version July 2015

#library(reshape2)
#library(ggplot2)
library(CellNOptR)

# *************** set working directory for relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#source("../utils/calculateScore.R")
source("../utils/calculateScore_MR.R")                                  # MR modified


# source("/Users/marti/Documents/r/combiMS/edgeContribution2.R")
# source("/Users/marti/Documents/r/combiMS/loadPatientRuns.R")

#***********************************************************************
# *************** preprocess model
#***********************************************************************
data_folder="../../data/phosphos_processed__original_MB_JW/"             # MR modified


patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
model_path="../../files/model/combiMSplaneCUT.sif"
#model_path="../../files/model/combiMS_rerun_whole_pkn_preprocessed_MR_211117.sif"                        # MR modified

fileName=patientData[1]

midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)  
#sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
sprintf("********* The used model has  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))                                        # MR modified


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





# START MR modified

files_Results_folder = "../../files/Results"

sub_dir__Results_single_model_optimization = file.path(files_Results_folder,"Results_single_model_optimization")

sub_dir__Results_single_model_optimization__completePatientsFinal_MR = file.path(sub_dir__Results_single_model_optimization,Results_storage_name)




files_Results_model_merging_MR_folder = "../../files/model_merging_MR"

ifelse(!dir.exists(file.path(files_Results_model_merging_MR_folder,Results_storage_name)), dir.create(file.path(files_Results_model_merging_MR_folder,Results_storage_name)), FALSE)
sub_dir__Results_model_merging_MR__current_script = file.path(files_Results_model_merging_MR_folder,Results_storage_name)



files_Results_median_models_MR_folder = "../../files/median_models_MR"

ifelse(!dir.exists(file.path(files_Results_median_models_MR_folder,Results_storage_name)), dir.create(file.path(files_Results_median_models_MR_folder,Results_storage_name)), FALSE)
sub_dir__Results_median_models_MR__current_script = file.path(files_Results_median_models_MR_folder,Results_storage_name)


# END MR modified












# networks_folder="../../files/all_solutions_models/"
networks_folder= sub_dir__Results_single_model_optimization__completePatientsFinal_MR  # MR modified

totalPatients=list.files(networks_folder,pattern="*.RData",full.names=FALSE)
totalPatientsSimple=sapply(totalPatients, function(x){
  strsplit(x,"\\.")[[1]][1]
})

num_patients=length(totalPatients)
#num_patients=3
allMedianNetworks=matrix(nrow=num_patients,ncol=numInteractions)          
DonorBestNetworks=list()           # MR modified

medianNetworkErrors=list()
#meanNetworkErrors=list()                                                     # MR modified for testing

allMeanNetworks=matrix(nrow=num_patients,ncol=numInteractions)                         
total_patientDF=data.frame(name=character(),
                           num_nws=vector(mode="numeric"),
                           avg_score=vector(mode="numeric"),
                           total_nws=vector(),
                           best_score=vector(),
                           merged_model_score=vector(),
                           nws_new_reltol=vector(),
                           convergence=character())



total_patientDF__with_wrong_results=data.frame(name=character(),   # MR inserted
                                               num_nws=vector(mode="numeric"),
                                               avg_score=vector(mode="numeric"),
                                               wrong_avg_score=vector(mode="numeric"),
                                               total_nws=vector(),
                                               best_score=vector(),
                                               wrong_best_score=vector(),
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
  load(paste(networks_folder,totalPatients[index],sep="/")) # this loads PatientResults         # MR modified
  
  
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
    #cat("dim networks does not macth num scores!!","\n")
    cat("dim networks does not match num scores!!","\n")   # MR modified
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
    
    
    # START MR modified
    if(length(which(PatientResults$AllScores<scoreTol)) == 1){
      
      ModelsInTol = as.matrix(t(ModelsInTol))              # tranform vector to matrix 
      
    }
    
    
    
    # clean repeated networks to keep only unique
    # 
    # test that unique function is keeping unique rows and not columns
    # 
    # > Mat
    # [,1] [,2] [,3]
    # [1,]    1    2    3
    # [2,]    1    2    3
    # [3,]   10    2    3
    # > unique(Mat)
    # [,1] [,2] [,3]
    # [1,]    1    2    3
    # [2,]   10    2    3
    # 
    # 
    
    # END MR modified
    
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
    medianNetworkErrors[[index]]=calculateScore_MR(model,midas,bString=medianNw)      # MR modified
    
    # create 3rd structure: 1 mean network 
    cat("calculating average of interactions","\n")
    meanNw=apply(ModelsInTol,2,mean)
    allMeanNetworks[index,]=meanNw
    
    #     # 4th structure:error for mean ns
    #cat("calculating score of median network","\n")
    #meanNetworkErrors[[index]]=calculateScore_MR(model,midas,bString=meanNw)     # MR modified for testing
    
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
    
    
    # START MR modified
    
    
    ##############################
    ## START MB original version:
    ## ##############################
    # #for all networks in one run calculate average error and best error
    # run_start_index=1
    # run_end_index=PatientResults$NwsInRuns[1]
    # i=1
    # cat("calculating avg score for each run","\n")
    # cat("calculating best score for each run","\n")
    # 
    # while (i <=num_success_runs){
    #    scores_one_run=PatientResults$AllScores[run_start_index:run_end_index]
    #    one_patientDF$avg_score[i]=mean(scores_one_run,na.rm=T)          # 3 of one_patientDF
    #    one_patientDF$best_score[i]=min(scores_one_run,na.rm=T)          # 5 of one_patientDF
    #    
    #    #move indexes
    #    run_start_index=run_start_index+PatientResults$NwsInRuns[i]
    #    run_end_index=run_end_index+PatientResults$NwsInRuns[i]
    #    i=i+1
    #    
    # }
    # 
    ##############################
    ## END MB original version:
    ## ##############################
    
    
    
    
    
    
    
    
    
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
    
    
    
    
    
    
    
    
    
    # #for all networks in one run calculate average error and best error
    # 
    one_patientDF__with_wrong_avg_score_AND_wrong_best_score = one_patientDF
    
    run_start_index__MB=1
    run_end_index__MB=PatientResults$NwsInRuns[1]
    ii=1
    cat("calculating avg score for each run using MB's method","\n")
    cat("calculating best score for each run using MB's method","\n")
    
    while (ii <=num_success_runs){
      scores_one_run__MB=PatientResults$AllScores[run_start_index__MB:run_end_index__MB]
      one_patientDF__with_wrong_avg_score_AND_wrong_best_score$wrong_avg_score[ii]=mean(scores_one_run__MB,na.rm=T)          # 3 of one_patientDF
      one_patientDF__with_wrong_avg_score_AND_wrong_best_score$wrong_best_score[ii]=min(scores_one_run__MB,na.rm=T)          # 5 of one_patientDF
      
      #move indexes
      run_start_index__MB=run_start_index__MB+PatientResults$NwsInRuns[ii]
      run_end_index__MB=run_end_index__MB+PatientResults$NwsInRuns[ii]
      ii=ii+1
      
    }
    
    
    
    
    
    
    # END MR modified
    
    
    total_patientDF=rbind(total_patientDF,one_patientDF) #concatenate this with total
    total_patientDF__with_wrong_results=rbind(total_patientDF__with_wrong_results,one_patientDF__with_wrong_avg_score_AND_wrong_best_score) # MR inserted
    
  }
  
  # START MR modified
  
  
  # create 5th structure: best networks (best of best with min score)
  
  
  if(length(which(PatientResults$AllScores==min(PatientResults$AllScores)))==1){
    
    DonorBestNetworks[[index]] = PatientResults$AllNwsPatient[which(PatientResults$AllScores==min(PatientResults$AllScores)),]    # MR modified
    
  }else{
    DonorBestNetworks[[index]] = unique(PatientResults$AllNwsPatient[which(PatientResults$AllScores==min(PatientResults$AllScores)),])    # MR modified
  }
  
  
  
  
  # END MR modified
  
  
}

if(dim(total_patientDF)[1]!=1690){stop("Hey!! Not all patients were added, is IB068 missing?")}



if(save_results_of_current_script ==T){
  
  
  
  #***********************************************************************
  # *************** save all models and metadata structures
  #***********************************************************************
  rownames(allMedianNetworks)=totalPatientsSimple
  colnames(allMedianNetworks)=model$reacID
  #save(allMedianNetworks,file="../../files/median_models_MR/allMedianModels.RData")
  save(allMedianNetworks,file=file.path(sub_dir__Results_median_models_MR__current_script,"allMedianModels.RData"))            # MR modified
  
  
  rownames(allMeanNetworks)=totalPatientsSimple
  colnames(allMeanNetworks)=model$reacID
  #save(allMeanNetworks,file="../../files/median_models_MR/allMeanModels.RData")
  save(allMeanNetworks,file=file.path(sub_dir__Results_median_models_MR__current_script,"allMeanModels.RData"))            # MR modified
  
  
  names(medianNetworkErrors)=totalPatientsSimple
  #save(medianNetworkErrors,file="../../files/median_models_MR/allErrors.RData")
  save(medianNetworkErrors,file=file.path(sub_dir__Results_median_models_MR__current_script,"allErrors.RData"))            # MR modified
  
  
  
  names(DonorBestNetworks)=totalPatientsSimple
  #save(DonorBestNetworks,file="../../files/median_models_MR/DonorBestNetworks.RData")
  save(DonorBestNetworks,file=file.path(sub_dir__Results_median_models_MR__current_script,"DonorBestNetworks.RData"))            # MR modified
  
  
  
  
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
  
  
  
  
  
  
  # START MR inserted
  
  # *************** calculate best score across runs
  #patient_names=unique(total_patientDF$name)
  for (ii in 1:length(patient_names)){
    all_best_scores__MB=vector(mode="numeric",length=10)
    all_best_scores__MB=total_patientDF__with_wrong_results$wrong_best_score[which(total_patientDF__with_wrong_results$name==patient_names[ii])]
    total_patientDF__with_wrong_results$wrong_best_of_bests[which(total_patientDF__with_wrong_results$name==patient_names[ii])]=rep(min(all_best_scores__MB),10)
  }
  # ************** add an identifier for runs
  #run_id=c("1","2","3","4","5","6","7","8","9","10")
  #for (ii in 1:length(patient_names)){
  # total_patientDF__with_wrong_results$run_id[which(total_patientDF__with_wrong_results$name==patient_names[ii])]=run_id
  #}
  
  total_patientDF__with_wrong_results$best_of_bests = total_patientDF$best_of_bests     # MR inserted
  total_patientDF__with_wrong_results$run_id = total_patientDF$run_id                    # MR inserted
  
  
  # END MR inserted
  
  
  
  
  # save(total_patientDF,file="../../files/model_merging_MR/statsModels.RData")
  save(total_patientDF,file=file.path(sub_dir__Results_model_merging_MR__current_script,"statsModels.RData"))            # MR modified
  save(total_patientDF__with_wrong_results,file=file.path(sub_dir__Results_model_merging_MR__current_script,"statsModels__with_wrong_results.RData"))            # MR modified
  
  
  
  
}

print("Script finished!")   # MR inserted



# plotMedianModels_MR__using_Cluster_MR_Results_based_on_using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__postCellNOptRupdate090318



calculate_data_frames_again__instead_of_loading_stored_ones = T #TRUE #FALSE         # MR inserted




# Some notes added
# 
# AND
# 
# paths adapted for CombiMS rerun
# 
# AND
# 
# some modifications with no influence on the results
# 
# AND
# 
# ERROR correction concerning the MSEs calculation 
# compare MSEs_MR...
# 
# 
# ERROR correction concerning the calculation of the Jaccard index 
# by using jaccNello_2 instead of jaccNello, 
# as jaccNello calculates the old / uncorrect Jaccard index
# But this ERROR correction has no influence on the result as it is just a test for donor CH003 and CH004 without storing




# 
# by 
# Melanie Rinas
# November 2017


# QC of combiMS models after fixing median strategy to achieve a single model. 
# When optimization did not converge due to excessive solutions within given relative tolerance of best model,
# optimisation was repeated until convergence was achieved 10 times. Then the resulting models were merged using mergeModels.R. 
#
# Here, we analyse and discard artefacts such as model performance being affected by model size, or by number of solutions grouped.
# The script also produces figure S2 on model QC. 
# This script must be run after mergeModels.R -or requires some vars to be loaded-.
# Marti Bernardo-Faura. June 2015. 

# Adapted for directory structure of the GitHub project by Jakob Wirbel in 2017
# 
# 


Results_storage_name = "OptCombiMSclusterFork10_MR__using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__postCellNOptRupdate090318"         # MR inserted
data_used = "phosphos_processed__original_MB_JW"          # MR inserted
#model_used =   "combiMS_rerun_whole_pkn_preprocessed_MR_211117.sif"                      # MR inserted
model_used =   "combiMSplaneCUT.sif"                      # MR inserted




Results_storage_path = paste("../../files/model_merging_MR/", Results_storage_name ,sep="")                        # MR inserted



ifelse(!dir.exists(file.path(Results_storage_path,"Figures")), dir.create(file.path(Results_storage_path,"Figures")), FALSE)
sub_dir__Figures = file.path(Results_storage_path,"Figures")
#sub_dir__Figures = paste("../../files/model_merging_MR/", Results_storage_name, "/Figures" ,sep="")                        # MR inserted




library(ggplot2)
library(ggrepel)
library(corrplot)
library(CellNOptR)

# *************** set working directory for relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../utils/newJaccard.R")       # MR modified
#source("../utils/newJaccard_MR.R")     # MR modified

#source("../utils/calculateScore.R")   # MR modified
source("../utils/calculateScore_MR.R")  # MR modified

# calculateScore.R was benchmarked by editing computeScoreT1multiple.R, 
# where the difference with computeScoreT1.R was to return not just the score but also MSE and num measurements
#source("/Users/marti/Documents/r/combiMS/computeScoreT1multiple.R")
#source("/Users/marti/Documents/r/combiMS/getFitmultipe.R")

#***********************************************************************
# *************** as a reminder, this is the metadata about the models
#***********************************************************************

# load("../../files/model_merging/statsModels.RData")
load(paste("../../files/model_merging_MR/", Results_storage_name ,"/statsModels.RData",sep=""))                        # MR modified
load(paste("../../files/model_merging_MR/", Results_storage_name ,"/statsModels__with_wrong_results.RData",sep=""))    # MR inserted; data frame is called total_patientDF__with_wrong_results after loading





names(total_patientDF)
# [1] "name"               "num_nws"            "avg_score"          "total_nws"          "best_score"         "merged_model_score" "nws_new_reltol"     "convergence"        "best_of_bests"      "run_id"  
total_patientDF[which(total_patientDF$name=="IB030"),]
# name num_nws avg_score total_nws best_score merged_model_score nws_new_reltol convergence best_of_bests run_id
# 671 IB030  281182 0.2041323   4437921  0.1989447          0.1859097           2838           3     0.1859091      1
# 672 IB030  570384 0.2382348   4437921  0.2356783          0.1859097           2838           3     0.1859091      2
# 673 IB030  150570 0.2056265   4437921  0.1859091          0.1859097           2838           3     0.1859091      3
# 674 IB030   38201 0.2296705   4437921  0.1859091          0.1859097           2838           3     0.1859091      4
# 675 IB030  444667 0.2356595   4437921  0.2324375          0.1859097           2838           3     0.1859091      5
# 676 IB030  586940 0.2100592   4437921  0.2070877          0.1859097           2838           3     0.1859091      6
# 677 IB030  610442 0.2378742   4437921  0.2356772          0.1859097           2838           3     0.1859091      7
# 678 IB030  691897 0.2378967   4437921  0.2356767          0.1859097           2838           3     0.1859091      8
# 679 IB030  552197 0.2097010   4437921  0.2070866          0.1859097           2838           3     0.1859091      9
# 680 IB030  511441 0.2101008   4437921  0.2085658          0.1859097           2838           3     0.1859091     10

#and these are the median models in allMedianNetworks loaded by mergeModels.R
# load("../../files/median_models/allMedianModels.RData")
load(paste("../../files/median_models_MR/", Results_storage_name ,"/allMedianModels.RData",sep=""))                        # MR modified

numInteractions=length(allMedianNetworks[1,])
numPatients=169

#***********************************************************************
# *************** calculate num edges in each median model
#***********************************************************************
edgesKept=apply(allMedianNetworks,1, function(x){
  length(which(x==1))
})

# START - MR inserted 

# edgesKept2 = rowSums(allMedianNetworks, na.rm = TRUE, dims = 1)
# difference__edgesKept__VS__edgesKept2 = edgesKept - edgesKept2
# min(difference__edgesKept__VS__edgesKept2)
# # [1] 0
# max(difference__edgesKept__VS__edgesKept2)
# # [1] 0

# END - MR inserted 

#***********************************************************************
# *************** calculate mse to investigate whether it is linked to num edges (mse does not account for num edges, as opposed to score)
#***********************************************************************
data_folder=paste("../../data/", data_used , "/",sep="")   # MR modified
patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
model_path=paste("../../files/model/", model_used, sep="")   # MR modified
fileName=patientData[1]

midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)  
#sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))         # MR modified
sprintf("*********The used model has  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))          # MR modified

# model=preprocessing(midas,model,expansion=FALSE)
# numInteractions=length(model$reacID)
# sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))

#midas MUST be that of right patient
pb = txtProgressBar(min=0, max=numPatients)
medianNetworkMSE=list()

MSEs_MR = NA*vector(mode = "numeric",length = numPatients)                                    # MR modified/inserted
N0_number_of_data_points_at_t0_MR = NA*vector(mode = "numeric",length = numPatients)          # MR modified/inserted



for (index in 1:numPatients){
  fileName=patientData[index]
  midas=CNOlist(paste(data_folder,fileName,sep=""))
  setTxtProgressBar(pb, (pb$getVal()+1))
  # cat("calculating MSE of median network, data is",fileName,"\n")    
  medianNetworkMSE[[index]]=calculateScore_MR(model=model,midas = midas,bString=allMedianNetworks[index,])    # MR modified
  dev.off()
  
  
  # START MR inserted
  
  MSEs_MR[index] = ((medianNetworkMSE[[index]]$deviationPen)*2) / (sum(!is.na(midas@signals[[1]])) + sum(!is.na(midas@signals[[2]])))
  N0_number_of_data_points_at_t0_MR[index] = sum(!is.na(midas@signals[[1]]))
  
  # END MR inserted
  
}

# MSEs = unlist(medianNetworkMSE)
MSEs=unlist(lapply(medianNetworkMSE,function(x){
  (x$deviationPen)/(x$numMeasurements)
}))


# START MR inserted

test_identical__MSEs__Marti__versus__MR = identical(MSEs,MSEs_MR)
test_identical__MSEs__Marti__versus__MR

donor_IDs=sapply(patientData, function(x){
  strsplit(x,"\\.")[[1]][1]
})

N0_matrix__number_of_data_points_at_t0_MR = as.matrix(N0_number_of_data_points_at_t0_MR)
rownames(N0_matrix__number_of_data_points_at_t0_MR) = donor_IDs

# END MR inserted







if(calculate_data_frames_again__instead_of_loading_stored_ones == TRUE){
  
  
  #*************************************************************
  # *************** create simpler metadata structure, rejecting info about all runs, keeping only best
  #***********************************************************************
  patient_names=unique(total_patientDF$name)
  
  mini_DF=data.frame(name=character(),
                     total_nws=vector(),
                     merged_model_score=vector(),
                     nws_new_reltol=vector(),
                     convergence=character(),
                     best_of_bests=vector())
  
  
  
  mini_DF_with_wrong_best_of_bests=data.frame(name=character(),           # MR inserted
                                              total_nws=vector(),
                                              merged_model_score=vector(),
                                              nws_new_reltol=vector(),
                                              convergence=character(),
                                              best_of_bests=vector(),
                                              wrong_best_of_bests=vector())
  
  
  # mini_DF=data.frame()
  # names(mini_DF)=names(total_patientDF[c(1,4,6,7,8,9)])
  for (i in 1:length(patient_names)){
    positions_patient=which(total_patientDF$name==patient_names[i])
    mini_DF=rbind(mini_DF, total_patientDF[positions_patient[1],c("name","total_nws","merged_model_score","nws_new_reltol","convergence","best_of_bests")])
    mini_DF_with_wrong_best_of_bests=rbind(mini_DF_with_wrong_best_of_bests, total_patientDF__with_wrong_results[positions_patient[1],c("name","total_nws","merged_model_score","nws_new_reltol","convergence","best_of_bests","wrong_best_of_bests")])   # MR inserted
  }
  mini_DF$edges=edgesKept
  mini_DF_with_wrong_best_of_bests$edges=edgesKept   # MR inserted
  # mini_DF$MSE=MSEs
  mini_DF$MSE_MR=MSEs_MR              # MR modified
  mini_DF_with_wrong_best_of_bests$MSE_MR=MSEs_MR              # MR modified
  
  # *************** Save intermediate results
  # save(mini_DF,file="../../files/model_merging/mini_DF.RData")
  save(mini_DF,file=paste(Results_storage_path,"/mini_DF.RData",sep=""))   # MR modified
  save(mini_DF_with_wrong_best_of_bests,file=paste(Results_storage_path,"/mini_DF_with_wrong_best_of_bests.RData",sep=""))   # MR modified
  
  
  
  
  
  
  
  # START MR inserted
  
  mini_DF__with_wrongMSE_AND_MSEs_MR = mini_DF
  mini_DF__with_wrongMSE_AND_MSEs_MR$MSEwrong = MSEs
  
  mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR = mini_DF_with_wrong_best_of_bests
  mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$MSEwrong = MSEs
  
  
  save(mini_DF__with_wrongMSE_AND_MSEs_MR,file=paste(Results_storage_path,"/mini_DF__with_wrongMSE_AND_MSEs_MR.RData",sep=""))   # MR modified
  save(mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR,file=paste(Results_storage_path,"/mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR.RData",sep=""))   # MR modified
  
  
  
  
  
}else{
  
  
  
  load(paste(Results_storage_path ,"/mini_DF.RData",sep=""))    # MR inserted; data frame is called mini_DF after loading
  load(paste(Results_storage_path ,"/mini_DF_with_wrong_best_of_bests.RData",sep=""))    # MR inserted; data frame is called mini_DF_with_wrong_best_of_bests after loading
  
  load(paste(Results_storage_path ,"/mini_DF__with_wrongMSE_AND_MSEs_MR.RData",sep=""))    # MR inserted; data frame is called mini_DF__with_wrongMSE_AND_MSEs_MR after loading
  load(paste(Results_storage_path ,"/mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR.RData",sep=""))    # MR inserted; data frame is called mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR after loading
  
  
  
  
}







#***********************************************************************
# *************** Analyzing differences between 
# ****************             (1) avg_score__versus__wrong_avg_score
# ****************             (2) best_score__versus__wrong_best_score
# ****************             (3) best_of_bests__versus__wrong_best_of_bests

# 
# ****************             (4) MSEs_MR__versus__MSEwrong
#
#***********************************************************************



diff__avg_score__minus__wrong_avg_score =  total_patientDF__with_wrong_results$avg_score - total_patientDF__with_wrong_results$wrong_avg_score
min_diff__avg_score__minus__wrong_avg_score = min(diff__avg_score__minus__wrong_avg_score)
max_diff__avg_score__minus__wrong_avg_score = max(diff__avg_score__minus__wrong_avg_score)


FigureWidth__diff = 7
FigureHeight__diff = 3.5

pdf(file.path(sub_dir__Figures,paste("plot__diff__avg_score__minus__wrong_avg_score.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

#plot(total_patientDF__with_wrong_results$avg_score,diff__avg_score__minus__wrong_avg_score,
barplot(diff__avg_score__minus__wrong_avg_score,
        width=0.45,
        xlab="Single optimization runs",
        ylab = "avg_score - wrong_avg_score",
        ylim = c( -(max(abs(diff__avg_score__minus__wrong_avg_score)))*1.1 ,  max(abs(diff__avg_score__minus__wrong_avg_score)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min_diff__avg_score__minus__wrong_avg_score,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max_diff__avg_score__minus__wrong_avg_score,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()




relative_diff__avg_score__minus__wrong_avg_score = diff__avg_score__minus__wrong_avg_score/total_patientDF__with_wrong_results$avg_score

pdf(file.path(sub_dir__Figures,paste("plot__relative_diff__avg_score__minus__wrong_avg_score.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

#plot(total_patientDF__with_wrong_results$avg_score,diff__avg_score__minus__wrong_avg_score,
barplot(relative_diff__avg_score__minus__wrong_avg_score,
        width=0.45,
        xlab="Single optimization runs",
        # ylab = "(avg_score - wrong_avg_score) / avg_score",
        ylab = "Relative difference avg_scores",
        ylim = c( -(max(abs(relative_diff__avg_score__minus__wrong_avg_score)))*1.1 ,  max(abs(relative_diff__avg_score__minus__wrong_avg_score)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min(relative_diff__avg_score__minus__wrong_avg_score),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max(relative_diff__avg_score__minus__wrong_avg_score),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()







diff__best_score__minus__wrong_best_score =  total_patientDF__with_wrong_results$best_score - total_patientDF__with_wrong_results$wrong_best_score
min_diff__best_score__minus__wrong_best_score = min(diff__best_score__minus__wrong_best_score)
max_diff__best_score__minus__wrong_best_score = max(diff__best_score__minus__wrong_best_score)


FigureWidth__diff = 7
FigureHeight__diff = 3.5



pdf(file.path(sub_dir__Figures,paste("plot__diff__best_score__minus__wrong_best_score.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

#plot(total_patientDF__with_wrong_results$best_score,diff__best_score__minus__wrong_best_score,
barplot(diff__best_score__minus__wrong_best_score,
        width=0.45,
        xlab="Single optimization runs",
        ylab = "best_score - wrong_best_score",
        ylim = c( -(max(abs(diff__best_score__minus__wrong_best_score)))*1.1 ,  max(abs(diff__best_score__minus__wrong_best_score)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min_diff__best_score__minus__wrong_best_score,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max_diff__best_score__minus__wrong_best_score,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()




relative_diff__best_score__minus__wrong_best_score = diff__best_score__minus__wrong_best_score/total_patientDF__with_wrong_results$best_score

pdf(file.path(sub_dir__Figures,paste("plot__relative_diff__best_score__minus__wrong_best_score.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

#plot(total_patientDF__with_wrong_results$best_score,diff__best_score__minus__wrong_best_score,
barplot(relative_diff__best_score__minus__wrong_best_score,
        width=0.45,
        xlab="Single optimization runs",
        #ylab = "(best_score - wrong_best_score) / best_score",
        ylab = "Relative difference best_scores",
        ylim = c( -(max(abs(relative_diff__best_score__minus__wrong_best_score)))*1.1 ,  max(abs(relative_diff__best_score__minus__wrong_best_score)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min(relative_diff__best_score__minus__wrong_best_score),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max(relative_diff__best_score__minus__wrong_best_score),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()









diff__best_of_bests__minus__wrong_best_of_bests =  mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$best_of_bests - mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$wrong_best_of_bests
min_diff__best_of_bests__minus__wrong_best_of_bests = min(diff__best_of_bests__minus__wrong_best_of_bests)
max_diff__best_of_bests__minus__wrong_best_of_bests = max(diff__best_of_bests__minus__wrong_best_of_bests)


FigureWidth__diff = 7
FigureHeight__diff = 3.5



pdf(file.path(sub_dir__Figures,paste("plot__diff__best_of_bests__minus__wrong_best_of_bests.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

#plot(mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$best_of_bests,diff__best_of_bests__minus__wrong_best_of_bests,
barplot(diff__best_of_bests__minus__wrong_best_of_bests,
        width=0.45,
        xlab="Donors",
        ylab = "best_of_bests - wrong_best_of_bests",
        ylim = c( -(max(abs(diff__best_of_bests__minus__wrong_best_of_bests)))*1.1 ,  max(abs(diff__best_of_bests__minus__wrong_best_of_bests)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min_diff__best_of_bests__minus__wrong_best_of_bests,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max_diff__best_of_bests__minus__wrong_best_of_bests,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()




relative_diff__best_of_bests__minus__wrong_best_of_bests = diff__best_of_bests__minus__wrong_best_of_bests/mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$best_of_bests

pdf(file.path(sub_dir__Figures,paste("plot__relative_diff__best_of_bests__minus__wrong_best_of_bests.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

#plot(mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$best_of_bests,diff__best_of_bests__minus__wrong_best_of_bests,
barplot(relative_diff__best_of_bests__minus__wrong_best_of_bests,
        width=0.45,
        xlab="Donors",
        #ylab = "(best_of_bests - wrong_best_of_bests) / best_of_bests",
        ylab = "Relative difference best_of_bests",
        ylim = c( -(max(abs(relative_diff__best_of_bests__minus__wrong_best_of_bests)))*1.1 ,  max(abs(relative_diff__best_of_bests__minus__wrong_best_of_bests)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min(relative_diff__best_of_bests__minus__wrong_best_of_bests),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max(relative_diff__best_of_bests__minus__wrong_best_of_bests),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()












diff__MSE_MR__minus__MSEwrong =  mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$MSE_MR - mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$MSEwrong
min_diff__MSE_MR__minus__MSEwrong = min(diff__MSE_MR__minus__MSEwrong)
max_diff__MSE_MR__minus__MSEwrong = max(diff__MSE_MR__minus__MSEwrong)


FigureWidth__diff = 7
FigureHeight__diff = 3.5



pdf(file.path(sub_dir__Figures,paste("plot__diff__MSE_MR__minus__MSEwrong.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

#plot(total_patientDF__with_wrong_results$best_of_bests,diff__MSE_MR__minus__MSEwrong,
barplot(diff__MSE_MR__minus__MSEwrong,
        width=0.45,
        xlab="Donors",
        ylab = "MSE_MR - MSEwrong",
        ylim = c( -(max(abs(diff__MSE_MR__minus__MSEwrong)))*1.1 ,  max(abs(diff__MSE_MR__minus__MSEwrong)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min_diff__MSE_MR__minus__MSEwrong,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max_diff__MSE_MR__minus__MSEwrong,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()




relative_diff__MSE_MR__minus__MSEwrong = diff__MSE_MR__minus__MSEwrong/mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$MSE_MR

pdf(file.path(sub_dir__Figures,paste("plot__relative_diff__MSE_MR__minus__MSEwrong.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

#plot(total_patientDF__with_wrong_results$best_of_bests,diff__MSE_MR__minus__MSEwrong,
barplot(relative_diff__MSE_MR__minus__MSEwrong,
        width=0.45,
        xlab="Donors",
        #ylab = "(MSE_MR - MSEwrong) / MSE_MR",
        ylab = "Relative difference MSEs",
        ylim = c( -(max(abs(relative_diff__MSE_MR__minus__MSEwrong)))*1.1 ,  max(abs(relative_diff__MSE_MR__minus__MSEwrong)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min(relative_diff__MSE_MR__minus__MSEwrong),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max(relative_diff__MSE_MR__minus__MSEwrong),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()






#***********************************************************************
# *************** plot linear: 
# 
# ****************             (1) avg_score__versus__wrong_avg_score
# ****************             (2) best_score__versus__wrong_best_score
# ****************             (3) best_of_bests__versus__wrong_best_of_bests

# 
# ****************             (4) MSEs_MR__versus__MSEwrong
#
#***********************************************************************

FigureWidth_geom_point = 3.5      # MR inserted
FigureHeight_geom_point= 3.5      # MR inserted 


w1=ggplot(data=total_patientDF__with_wrong_results,aes(avg_score,wrong_avg_score))+
  geom_point(shape=1)+
  geom_smooth(method=lm)+                                                               # By default, it is the 95% confidence level interval for predictions from a linear model ("lm"). 
  geom_abline(slope = 1,intercept = 0)   # MR inserted

# and label those patients that are from diagonal
w1 +  
  #geom_text(aes(label=ifelse(abs(best_of_bests-merged_model_score)>0.01,as.character(name),'')),hjust=0,vjust=1) + 
  theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__total_patientDF__with_wrong_results__avg_score__versus__wrong_avg_score.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted







x<-total_patientDF__with_wrong_results$avg_score
ix_or<-order(x)     # order x = all_mean by value 
x<-x[ix_or]
y<-total_patientDF__with_wrong_results$wrong_avg_score[ix_or]

nls__avg_score = nls(y ~ m__avg * x + b__avg)         # nls = Nonlinear Least Squares
# nls__avg_score
# # Nonlinear regression model
# # model: y ~ m__avg * x + b__avg
# # data: parent.frame()
# # m__avg    b__avg 
# # 0.9989067 0.0003914 
# # residual sum-of-squares: 0.002674
# # 
# # Number of iterations to convergence: 1 
# # Achieved convergence tolerance: 5.282e-08

# w1_with_nls=ggplot(data=total_patientDF__with_wrong_results,aes(avg_score,wrong_avg_score))+
#    geom_point(shape=1)+
#    geom_line(color='red',data = predict(nls__avg_score))+
#    #geom_smooth(method=lm)+
#    geom_abline(slope = 1,intercept = 0)   # MR inserted
# 
# # and label those patients that are from diagonal
# w1_with_nls +  
#    #geom_text(aes(label=ifelse(abs(best_of_bests-merged_model_score)>0.01,as.character(name),'')),hjust=0,vjust=1) + 
#    theme_bw()
# ggsave(file.path(sub_dir__Figures,paste("geom_point_with_nls__total_patientDF__with_wrong_results__avg_score__versus__wrong_avg_score.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted
# 







w2=ggplot(data=total_patientDF__with_wrong_results,aes(best_score,wrong_best_score))+
  geom_point(shape=1)+
  geom_smooth(method=lm)+
  geom_abline(slope = 1,intercept = 0)   # MR inserted

# and label those patients that are from diagonal
w2 +  
  #geom_text(aes(label=ifelse(abs(best_of_bests-merged_model_score)>0.01,as.character(name),'')),hjust=0,vjust=1) + 
  theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__total_patientDF__with_wrong_results__best_score__versus__wrong_best_score.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted






w3=ggplot(data=mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR,aes(best_of_bests,wrong_best_of_bests))+
  geom_point(shape=1)+
  geom_smooth(method=lm)+
  geom_abline(slope = 1,intercept = 0)   # MR inserted

# and label those patients that are from diagonal
w3 +  
  #geom_text(aes(label=ifelse(abs(best_of_bests-merged_model_score)>0.01,as.character(name),'')),hjust=0,vjust=1) + 
  theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR__best_of_bests__versus__wrong_best_of_bests.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted













w4=ggplot(data=mini_DF__with_wrongMSE_AND_MSEs_MR,aes(MSEs_MR,MSEwrong))+
  geom_point(shape=1)+
  geom_smooth(method=lm)+
  geom_abline(slope = 1,intercept = 0)   # MR inserted

# and label those patients that are from diagonal
w4 +  
  #geom_text(aes(label=ifelse(abs(best_of_bests-merged_model_score)>0.01,as.character(name),'')),hjust=0,vjust=1) + 
  theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__mini_DF__with_wrongMSE_AND_MSEs_MR__MSEs_MR__versus__MSEwrong.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted











# END MR inserted


#***********************************************************************
# *************** plot distributions
#***********************************************************************
#

FigureWidth_geom_bar = 20      # MR inserted
FigureHeight_geom_bar=3       # MR inserted 

#merged
ggplot(data=mini_DF, aes(reorder(name, merged_model_score), merged_model_score,fill=convergence))+
  geom_bar(stat="identity") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(sub_dir__Figures,paste("geom_bar__mini_DF__merged_model_score.pdf",sep = "")), width = FigureWidth_geom_bar, height = FigureHeight_geom_bar)             # MR inserted


#best
ggplot(data=mini_DF, aes(reorder(name, merged_model_score), best_of_bests,fill=convergence))+
  geom_bar(stat="identity") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(sub_dir__Figures,paste("geom_bar__mini_DF__best_of_bests.pdf",sep = "")), width = FigureWidth_geom_bar, height = FigureHeight_geom_bar)             # MR inserted





#***********************************************************************
# *************** plot linear model best error vs merged model error
#***********************************************************************

FigureWidth_geom_point = 3.5      # MR inserted
FigureHeight_geom_point= 3.5      # MR inserted 


a=ggplot(data=mini_DF,aes(best_of_bests,merged_model_score))+
  geom_point(shape=1)+
  geom_smooth(method=lm)+
  geom_abline(slope = 1,intercept = 0)   # MR inserted

# and label those patients that are from diagonal
a +  geom_text(aes(label=ifelse(abs(best_of_bests-merged_model_score)>0.01,as.character(name),'')),hjust=0,vjust=1) + 
  theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__mini_DF__best_of_bests__versus__merged_model_score.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted






diff__merged_model_score__minus__best_of_bests =  mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$merged_model_score - mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$best_of_bests
min_diff__merged_model_score__minus__best_of_bests = min(diff__merged_model_score__minus__best_of_bests)
max_diff__merged_model_score__minus__best_of_bests = max(diff__merged_model_score__minus__best_of_bests)


FigureWidth__diff = 7
FigureHeight__diff = 3.5



pdf(file.path(sub_dir__Figures,paste("plot__diff__merged_model_score__minus__best_of_bests.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

barplot(diff__merged_model_score__minus__best_of_bests,
        width=0.45,
        xlab="Donors",
        ylab = "merged_model_score - best_of_bests",
        ylim = c( -(max(abs(diff__merged_model_score__minus__best_of_bests)))*1.1 ,  max(abs(diff__merged_model_score__minus__best_of_bests)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min_diff__merged_model_score__minus__best_of_bests,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max_diff__merged_model_score__minus__best_of_bests,b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()




relative_diff__merged_model_score__minus__best_of_bests = diff__merged_model_score__minus__best_of_bests/mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR$merged_model_score

pdf(file.path(sub_dir__Figures,paste("plot__relative_diff__merged_model_score__minus__best_of_bests.pdf",sep = "")),width = FigureWidth__diff, height = FigureHeight__diff)

barplot(relative_diff__merged_model_score__minus__best_of_bests,
        width=0.45,
        xlab="Donors",
        #ylab = "(MSE_MR - MSEwrong) / MSE_MR",
        ylab = "Relative difference merged_model_score minus best_of_bests",
        ylim = c( -(max(abs(relative_diff__merged_model_score__minus__best_of_bests)))*1.1 ,  max(abs(relative_diff__merged_model_score__minus__best_of_bests)))*1.1 )
abline(h = 0,col = "green4",lwd=2)  # add horizontal line
abline(a=min(relative_diff__merged_model_score__minus__best_of_bests),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line
abline(a=max(relative_diff__merged_model_score__minus__best_of_bests),b=0,col = "darkorange2",lwd=2,lty=2)  # add horizontal line

dev.off()












FigureWidth_geom_point = 3.5  +2    # MR inserted
FigureHeight_geom_point= 3.5      # MR inserted 

# scatter plot to color by convergence
b=ggplot(data=mini_DF,aes(best_of_bests,merged_model_score, color=convergence))+
  geom_point(shape=1, size=4)
b+theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__mini_DF__best_of_bests__versus__merged_model_score__colored_by_convergence.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted

#b +scale_colour_gradient(low = "blue")





#***********************************************************************
# *************** for all patients together: plot linear model error vs total num networks
#***********************************************************************

FigureWidth_geom_point = 3.5  +2    # MR inserted
FigureHeight_geom_point= 3.5      # MR inserted

c= ggplot(data=mini_DF,aes(best_of_bests,total_nws))+
  #geom_point(size=4,color="grey")+
  geom_point(size=2,color="grey")+      # MR modified
  geom_smooth(method=lm)

c+theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__mini_DF__best_of_bests__versus__total__nws.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted


# plot distributions 
d= ggplot(data=mini_DF,aes(reorder(name,total_nws),total_nws, fill=convergence))+
  geom_bar(stat="identity")

d  + theme_bw() +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(sub_dir__Figures,paste("geom_bar__mini_DF__best_of_bests__versus__total__nws.pdf",sep = "")), width = FigureWidth_geom_bar, height = FigureHeight_geom_bar)             # MR inserted


#***********************************************************************
# *************** for each patient: plot linear model error vs total num networks
#***********************************************************************


FigureWidth_facet_wrap = 24    # MR inserted
FigureHeight_facet_wrap= 10     # MR inserted



e= ggplot(data=total_patientDF,aes(best_score,num_nws,color=name))+
  #geom_point(size=4)+
  geom_point(size=1)+                                                    # MR modified
  guides(color=FALSE)#+geom_smooth(method=lm)
# e=e+theme_bw()+facet_grid(name~)
# e+theme_bw()
# 
e+facet_wrap(~name,ncol=17)+theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__total_patientDF__best_score__versus__num_nws.pdf",sep = "")), width = FigureWidth_facet_wrap, height = FigureHeight_facet_wrap)             # MR inserted



e_wrong= ggplot(data=total_patientDF__with_wrong_results,aes(wrong_best_score,num_nws,color=name))+
  #geom_point(size=4)+
  geom_point(size=1)+                                                    # MR modified
  guides(color=FALSE)#+geom_smooth(method=lm)
# e=e+theme_bw()+facet_grid(name~)
# e+theme_bw()
# 
e_wrong+facet_wrap(~name,ncol=17)+theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__total_patientDF__with_wrong_results__wrong_best_score__versus__num_nws.pdf",sep = "")), width = FigureWidth_facet_wrap, height = FigureHeight_facet_wrap)             # MR inserted





# plot distributions 
f= ggplot(data=total_patientDF,aes(reorder(run_id,num_nws),num_nws, fill=convergence))+
  geom_bar(stat="identity")

f=f  + theme_bw() +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
f+facet_wrap(~name,ncol=17)
ggsave(file.path(sub_dir__Figures,paste("geom_bar__total_patientDF__num_nws.pdf",sep = "")), width = FigureWidth_facet_wrap, height = FigureHeight_facet_wrap)             # MR inserted






#***********************************************************************
# *************** plot num edges vs MSE
#***********************************************************************
#g= ggplot(data=mini_DF,aes(edges,MSE))+
g= ggplot(data=mini_DF,aes(edges,MSE_MR))+                  # MR modified
  geom_smooth(method=lm)+
  #geom_point(size=4)+
  geom_point(size=2)+                                      # MR modified
  theme_bw()

g+theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_point__mini_DF__edges__versus__MSE_MR.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted






g= ggplot(data=mini_DF,aes(edges))#+geom_bar()+theme_bw()

g+geom_density(fill=NA,aes(y=..count..))+
  geom_histogram(aes(y=..count..))+
  theme_bw()
ggsave(file.path(sub_dir__Figures,paste("geom_density__mini_DF__edges.pdf",sep = "")), width = FigureWidth_geom_point, height = FigureHeight_geom_point)             # MR inserted




#g=ggplot(data=mini_DF, aes(reorder(name, edges), edges,fill=MSE)) #+ scale_fill_manual(values=Greys)
g=ggplot(data=mini_DF, aes(reorder(name, edges), edges,fill=MSE_MR)) #+ scale_fill_manual(values=Greys)                     # MR modified

g=g+geom_bar(stat="identity") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
g + scale_fill_continuous(low='grey',high='grey11')
ggsave(file.path(sub_dir__Figures,paste("geom_bar__mini_DF__edges__filled_by_MSE_MR.pdf",sep = "")), width = FigureWidth_geom_bar, height = FigureHeight_geom_bar)             # MR inserted







#***********************************************************************
# *************** calculate similarity between patients with new median network
#***********************************************************************
#distanceModels=jaccNello(allMedianNetworks[1,],allMedianNetworks[2,])
distanceModels=jaccNello_2(allMedianNetworks[1,],allMedianNetworks[2,])     # MR modified, as the function jaccNello calculates the old / uncorrect Jaccard index



print("Script finished!")   # MR inserted

#***********************************************************************
# solve problem: why are there patients below the diagonal in plot "a" best_of_best vs merged_model_score?
# this are models better than found during optimization!
# answer found: this was a bug in computeScore, now fixed and benchmarked against own version of computeScoreT1.R
#***********************************************************************
# 
# mini_DF[which(mini_DF$best_of_bests-mini_DF$merged_model_score>0.01),'name']
# # "IB064" "IB069" "UZ010" "UZ011" "UZ017"
# # from here, it runs again
# networks_folder='../../files/all_solutions_models/'# "/Users/marti/Documents/R/combiMS/cluster/all/"
# patient="IB064.RData"
# load(paste0(networks_folder,patient)) # this loads PatientResults
# which(PatientResults$AllScores==min(PatientResults$AllScores)) # there are several models with the exaxt same score!
# #372303 372696 372706 372754 373726 373811 374029 
# # get the first model out of the 7 with (same) best score
# 
# best_model=PatientResults$AllNwsPatient[which(PatientResults$AllScores==min(PatientResults$AllScores))[1],]
# #plot simulation and data
# fileName=strsplit(patient,"\\.")[[1]][1]
# midas=CNOlist(paste0(data_folder,fileName,".csv"))
# model=readSIF(model_path)  
# sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
# #model=preprocessing(midas,model,expansion=FALSE)
# #numInteractions=length(model$reacID)
# #sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
# # Here it fails again
# prova=calculateScore(model = model,midas = midas,bString=best_model)
# #plotModel(model,midas,bString=best_model)
# prova=computeScoreT1(CNOlist=midas,model=model,bString=best_model)
# computeScoreT1multiple(CNOlist=midas,model=model,bString=best_model)

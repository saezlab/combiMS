
#loadPatientRuns_MR.R




# For a single patient, merge cluster jobs until 10 successful cluster jobs have been concatenated
# Called by completePatientsFinal.R
# Marti Bernardo-Faura
# June 2015








# with some notes added by 
# Melanie Rinas
# November 2017

# AND the following modification
#    # AllNwsPatient=matrix(ncol=numInteractions)    # original
# AllNwsPatient= c()                                 # MR modified reason: 
# initialization using matrix(ncol=numInteractions) creates a NA line 
# which is later the first line in networksToAdd$AllNwsPatient
# because of AllNwsPatient=rbind(AllNwsPatient,network$Opt$stringsTol)
















# # START - MR inserted 
# 
# 
# # for testing
# 
# patient_name = "CH003"
# runsNeededPatient = 10
# numInteractions = 186
# networks_for_Re_Completion = "/Users/melanie/Documents/GitHub/combiMS-master_MR/files/Cluster_MR/Results/OptCombiMSclusterFork10_MR__using_MIDAS_phosphos_processed__original_MB_JW__AND__optimization_fct__gaBinaryT1__with_set_seed_number/CombiMS_rerun_all_donors/"
# 
# # END - MR inserted 
# 




loadPatientRuns_MR=function(patient_name,runsNeededPatient,numInteractions,networks_for_Re_Completion){
   
   
   
   
   
   
   # generate a list of all runs for given incomplete patient
   #  cat("the evil folder is",networks_for_Re_Completion,"\n")
   modelsForCompletion=list.files(networks_for_Re_Completion,pattern="*.RData",full.names=FALSE)
   # cat("num files total",length(modelsForCompletion),"\n")
   modelsForCompletionAllPatient=modelsForCompletion[grep(patient_name,modelsForCompletion)]
   
   networksToAdd=list()
   cat("fetching all Nws for all runs for patient",sep="\n")
   cat(patient_name,sep="\n")
   
   
   # AllNwsPatient=matrix(ncol=numInteractions)
   AllNwsPatient= c()                                 # MR modified reason: 
                                                      # initialization using matrix(ncol=numInteractions) creates a NA line 
                                                      # which is later the first line in networksToAdd$AllNwsPatient
                                                      # because of AllNwsPatient=rbind(AllNwsPatient,network$Opt$stringsTol)
   
   
   NwsInRuns=vector()
   AllScores=vector()
   cat("runs needed for patient",runsNeededPatient,"\n")
   runsAvailable=length(modelsForCompletionAllPatient)
   cat("runs available for patient",runsAvailable,"\n")
   
   if(runsAvailable<runsNeededPatient){
      cat("**************","\n")
      cat("Ups! Patient lacks ",runsNeededPatient-runsAvailable,"runs","\n")
      cat("**************","\n")
      
   }
   successfulTries=0
   j=1;
   while ( (successfulTries<runsNeededPatient) && (j<=runsAvailable)){
      fileName=modelsForCompletionAllPatient[j]
      
      cat("******* loading data ")
      cat(modelsForCompletionAllPatient[j],sep="\n")
      
      
      
      result=tryCatch(
         {load(paste(networks_for_Re_Completion,fileName,sep="")); successfulTries=successfulTries+1; NwsInRuns=c(NwsInRuns,dim(network$Opt$stringsTol)[1]); AllNwsPatient=rbind(AllNwsPatient,network$Opt$stringsTol); AllScores=c(AllScores,network$Opt$stringsTolScores)},
         error=function(e) {cat("error!!!!!!!!!! Data could not be loaded",sep="\n")})
      
      
      cat("runs tried",j,"\n")
      cat("runs completed",successfulTries,"\n")
      
      
      j=j+1
   }
   networksToAdd$AllNwsPatient=AllNwsPatient
   networksToAdd$AllScores=AllScores
   #here the sequence of number of solutions for each run
   networksToAdd$NwsInRuns=NwsInRuns
   networksToAdd$name=patient_name
   
   return(networksToAdd)
}
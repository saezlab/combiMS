# For a single patient, merge cluster jobs until 10 successful cluster jobs have been concatenated
# Called by completePatientsFinal.R
# Marti Bernardo-Faura
# June 2015




loadPatientRuns=function(patient_name,runsNeededPatient,numInteractions,networks_for_Re_Completion){
  # generate a list of all runs for given incomplete patient
#  cat("the evil folder is",networks_for_Re_Completion,"\n")
  modelsForCompletion=list.files(networks_for_Re_Completion,pattern="*.RData",full.names=FALSE)
 # cat("num files total",length(modelsForCompletion),"\n")
  modelsForCompletionAllPatient=modelsForCompletion[grep(patient_name,modelsForCompletion)]
  
  networksToAdd=list()
  cat("fetching all Nws for all runs for patient",sep="\n")
  cat(patient_name,sep="\n")
  AllNwsPatient=matrix(ncol=numInteractions)
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
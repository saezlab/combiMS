# Function to sanity check if concatenning the networks resulting of optimization for 10 successful cluster
# runs produces a valid merge of all networks within relative tolerance of the best model in the 10 cluster runs
# Called by completePatientsFinal.R
# Marti Bernardo-Faura, June 2015

checkNwsAndScores=function(PatientResults,recalculated_score){
  #This function tests if removing first NA line messed up, and in general if results in structure are consistent
  all_fine=T
  #compare new scores with previously calculated in cluster along with optimisation
  if(! all(PatientResults$AllScores[1:2]==recalculated_score)){
    cat("problem: content scores do not match nws ***********************","\n")
    all_fine=F  
  }
  
  #check if num of scores and num of networks matches
  if(! length(PatientResults$AllScores) == dim(PatientResults$AllNwsPatient)[1]){
    cat("problem: num networks does not match num scores ******","\n")
    all_fine=F
  }
  
  if(! dim(PatientResults$AllNwsPatient)[1] == sum(PatientResults$NwsInRuns)){
    cat("problem: num networks does not matches that in runs *****","\n")
    all_fine=F
  }
  
  if(all_fine==T){
    cat("network dimension and score test: success!","\n")
  }
}
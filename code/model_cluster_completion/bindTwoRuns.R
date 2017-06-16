# For a single patient, merge cluster jobs until 10 successful cluster jobs have been concatenated
# Called by completePatientsFinal.R
# Marti Bernardo-Faura
# June 2015



bindTwoRuns=function(nws_and_scores_first_run,nws_and_scores_second_run){
  #***********************************************************************
  # *************** bind networks, scores, and num networks from all runs
  #***********************************************************************
  completed=list()
  completed$name=strsplit(nws_and_scores_first_run$name,"\\.")[[1]][1]
  
  completed$AllNwsPatient=rbind(nws_and_scores_first_run$AllNwsPatient,nws_and_scores_second_run$AllNwsPatient)
  
  completed$AllScores=c(nws_and_scores_first_run$AllScores,nws_and_scores_second_run$AllScores)
  completed$NwsInRuns=c(nws_and_scores_first_run$NwsInRuns,nws_and_scores_second_run$NwsInRuns)
  
  return(completed)
  
}
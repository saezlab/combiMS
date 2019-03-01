#***********************************************************************
# This scripts loads the results of loadAllModelsPlaneNW2.R,
# which is the concatenation of all networks solutions within 5% reltol for 
# the runs that did not crush out of 10 attempts per patient
#***********************************************************************

# Use relative paths instead of absolut paths for the files
# in Rstudio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('~/Documents/combiMS/code/similarity/')

library(pryr)

source("./SMC_of_2networks.R")

argsJob= commandArgs(trailingOnly=TRUE)

patients = read.csv('../../files/annot169pat_v2.csv')$ID

patientOneOriginal=toString(patients[as.numeric(argsJob[1])])
patientTwoOriginal=toString(patients[as.numeric(argsJob[2])])

#******************************************************************
#*************** define directory, load each patient's median model
#***********************************************************************

networks_folder="../../all_solutions_models/"
similarity_folder="../../files/similarity/"
load(paste0(networks_folder,patientOneOriginal,".RData"))
patientOne=apply(PatientResults$AllNwsPatient, 2, median)
load(paste0(networks_folder,patientTwoOriginal,".RData"))
patientTwo=apply(PatientResults$AllNwsPatient, 2, median)

#***********************************************************************
# *************** calculate similarity
#***********************************************************************

distanceModels=jaccNello(patientOne,patientTwo)

patientPair=paste0(patientOneOriginal,"_",patientTwoOriginal,".RData")
cat("----------saving similarity for",patientPair,"\n")
save(distanceModels,file=paste0(similarity_folder,patientPair))

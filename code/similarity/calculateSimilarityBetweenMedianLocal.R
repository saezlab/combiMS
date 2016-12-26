#***********************************************************************
# This scripts loads the results of loadAllModelsPlaneNW2.R,
# which is the concatenation of all networks solutions within 5% reltol for 
# the runs that did not crush out of 10 attempts per patient
#***********************************************************************
#library(reshape2)
#library(ggplot2)
# library(CellNOptR)
#library(corrplot)
#library(gdata)

#source("/homes/bernardo/combiMS/modelling/jaccNello.R")
source("/Users/marti/Documents/r/combiMS/jaccNello.R")

argsJob= commandArgs(trailingOnly=TRUE)
patientOneOriginal=toString(argsJob[1])
patientTwoOriginal=toString(argsJob[2])

patientOneOriginal="CH003"
patientTwoOriginal="CH004"


#******************************************************************
#*************** define directory, load each patient's median model
#***********************************************************************
#networks_folder="/gpfs/nobackup/saezgrp/bernardo/combiMS/allNetworks"
networks_folder="/Users/marti/Documents/R/combiMS/cluster/singleMedianModels/"
#networks_folder="/gpfs/nobackup/saezgrp/bernardo/combiMS/medianModels/"
#similarity_folder="/homes/bernardo/combiMS/modelling/between_patient_similarity/"
#similarity_folder="/gpfs/nobackup/saezgrp/bernardo/combiMS/median_patient_similarity/"
similarity_folder="/Users/marti/Documents/R/combiMS/cluster/similarity/median/"
load(paste0(networks_folder,patientOneOriginal,".RData"))
patientOne=medianModel
load(paste0(networks_folder,patientTwoOriginal,".RData"))
patientTwo=medianModel
#***********************************************************************
# *************** calculate similarity
#***********************************************************************

distanceModels=jaccNello(patientOne,patientTwo)
#distanceModels=as.vector(distanceModels)
# distance_stats=summary(distanceModels)
# distance_stats=append(distance_stats,sd(distanceModels))
# names(distance_stats)=c("Min","1stQu","Median","Mean","3rdQu","Max","std")

# distanceModels=rbind(TolPatientOne[1:2,],TolPatientTwo[1:2,])
patientPair=paste0(patientOneOriginal,"_",patientTwoOriginal,".RData")
cat("----------saving similarity for",patientPair,"\n")
save(distanceModels,file=paste0(similarity_folder,patientPair))






#hist(lowerTriangle(distanceModels))








#write.table(distanceModels,file=paste0(similarity_folder,patientName),sep=",")

#prune networks to remove spurious interactions
# prova=edgeContribution2(ModelsInTol[i,],model,midas)
# are all networks unique? this seems to indicate too large solution space, as they should duplicate through runs
#duplicated(ModelsInTol)


#subsetModels=ModelsInTol[1:10,]
#distanceModels=vegdist(subsetModels, method="jaccard")

#corrplot.mixed(as.matrix(distanceModels),is.corr=F)
#transform intor string, find out which ones are the same
#stringSolutions=apply(PatientResults$AllNwsPatient, 1, function(x) paste(x, collapse=","))


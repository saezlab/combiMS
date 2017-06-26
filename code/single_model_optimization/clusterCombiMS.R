# Script to parallelize at EBI's computer cluster the optimization for all patients
# Marti Bernardo Faura, final version Juni 2015

library(CellNOptR)
source("/homes/bernardo/combiMS/modelling/OptCombiMScluster.R")

# ****************** to parallelize optimisation in the cluster, this scritp receives as an argument the patient
argsJob= commandArgs(trailingOnly=TRUE)
fileName=toString(argsJob[1])
repIndex=as.numeric(argsJob[2])
# ****************** use data from each patients to optimize one model per patient
#midas=CNOlist("../processedForCNO/normalized/IB015.csv")
#setwd("/homes/bernardo/combiMS/processedForCNO/normalized/")
#filenames=list.files("/homes/bernardo/combiMS/processedForCNO/normalized/",pattern="*.csv",full.names=TRUE)
#   
# networks=list()
# for (i in 1:length(filenames)){
#   
#   networks[[i]]=OptCombiMScluster(filenames[i])
# }


cat("*********** the patient is",fileName,sep="\n")
cat("*********** the repetition is",repIndex,sep="\n")
  
network=OptCombiMScluster(fileName)

#**************

fileName2=unlist(strsplit(fileName,"\\."))[1]
cat("*********** saving patient ",fileName2,sep="\n")
save(network,file=paste("/homes/bernardo/combiMS/modelling/OnePatientForCompletion10hReltol005/results/",repIndex,"_",fileName2,".RData",sep=""))

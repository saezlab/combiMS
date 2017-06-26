# Script to parallelize at EBI's computer cluster the optimization for all patients
# Marti Bernardo Faura, final version Juni 2015

# The script is called with two arguments: the file with the normalised patient data and the index of the optimization, e.g.
#   R clusterCombiMS.R CH003.csv 8

# Minor changes to adjust the paths to the combiMS Github project
# Jakob Wirbel, June 2017

# set working directory
setwd(getSrcDirectory()[1])

library(CellNOptR)
source("./OptCombiMScluster.R")

# ****************** to parallelize optimisation in the cluster, this script receives as an argument the patient
argsJob = commandArgs(trailingOnly=TRUE)
fileName = toString(argsJob[1])
repIndex = as.numeric(argsJob[2])
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
save(network,file=paste("../../files/modeling/cluster/",repIndex,"_",fileName2,".RData",sep=""))

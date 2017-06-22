# ######################################################################################################
# NORMALIZATION PIPELINE FOR COMBIMS
#
# Compiled from OptAllPatients.R
#               processCombiData.R
#               normaliseSimp.R
#
# ######################################################################################################

# Load libraries
library(reshape2)
library(CellNOptR)
library(NMF)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# STILL TO DO
# ****************** FUNCTION FOR COMBINATION OF BOTH TIMEPOINTS
# filenames = list.files('../../data/phosphos_raw_clean', pattern='*.csv', full.names=FALSE)


# ****************** process merged files
filenames=list.files("../../data/phosphos_merged/",pattern="*.csv",full.names=FALSE)
source('./processCombiData.R')

for (i in 1:length(filenames)){
  processCombiData(filenames[i]) 
}

rm(list=ls())
# ****************** normalise all files
filenames=list.files("../../data/phosphos_normalised/",pattern="*.csv",full.names=FALSE)
source('./normaliseSimp.R')

for (i in 1:length(filenames)){
  fileName=filenames[i]
  cat(sprintf("*********The patient is %s\n",strsplit(fileName, split='\\.')[[1]][1]))
  taula=read.csv(paste("../../data/phosphos_normalised/",fileName,sep=""),header=TRUE,dec=".",check.names=FALSE, stringsAsFactors=FALSE)
  taula=normaliseSimp(taula)
  write.table(taula,file=paste0('../../data/phosphos_normalised/', fileName),sep=",",row.names=FALSE,quote=F)
}

rm(list=ls())
# ****************** further processing for matching names between midas and model
filenames=list.files("../../data/phosphos_normalised/",pattern="*.csv",full.names=FALSE)
for (i in 1:length(filenames)){
  fileName=filenames[i]
  cat(sprintf("*********The patient is %s\n",strsplit(fileName, split='\\.')[[1]][1]))
  taula=read.csv(paste("../../data/phosphos_normalised/",fileName,sep=""),header=TRUE,dec=".",check.names=FALSE, stringsAsFactors=FALSE)
  colnames(taula)[2:dim(taula)[2]] = toupper(colnames(taula)[2:dim(taula)[2]])
  write.table(taula,file=paste0('../../data/phosphos_processed/', fileName),sep=",",row.names=FALSE,quote=F)
}

rm(list=ls())
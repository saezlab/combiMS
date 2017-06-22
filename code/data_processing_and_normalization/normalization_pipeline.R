# ######################################################################################################
# NORMALIZATION PIPELINE FOR COMBIMS
#
# Compiled from OptAllPatients.R
#               OptCombiMS.R
#               processCombiData.R
#
#
# ######################################################################################################

# Load libraries
library(reshape2)
library(CellNOptR)
library(NMF)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# STILL TO DO
# # Find data files
# file_list = list.files('../../data/phosphos_raw_clean/')
# 
# ****************** FUNCTION FOR COMBINATION OF BOTH TIMEPOINTS
# combine_time_points = function(file_name){}



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


# STILL TO DO
# ******************FUNCTION FOR FURTHER PROCESSING
# process_normalized = function(file_name){}

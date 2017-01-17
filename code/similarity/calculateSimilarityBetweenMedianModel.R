#***********************************************************************
# This scripts loads the results of loadAllModelsPlaneNW2.R,
# which is the concatenation of all networks solutions within 5% reltol for 
# the runs that did not crush out of 10 attempts per patient
#***********************************************************************

# Use relative paths instead of absolut paths for the files
# in Rstudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd('~/Documents/combiMS/code/similarity/')

source("./newJaccard.R")

# load median models
load('../../files/allMedianModels.RData')

# compute pairwise similarity 
similarityMatrix = sapply(unname(row.names(allMedianNetworks)), function(x){sapply(unname(row.names(allMedianNetworks)), function(y){return(jaccNello(allMedianNetworks[x,], allMedianNetworks[y,]))})})
diag(similarityMatrix) = NA

# save similarity Matrix
save(similarityMatrix, file='../../files/similarityMatrix.RData')

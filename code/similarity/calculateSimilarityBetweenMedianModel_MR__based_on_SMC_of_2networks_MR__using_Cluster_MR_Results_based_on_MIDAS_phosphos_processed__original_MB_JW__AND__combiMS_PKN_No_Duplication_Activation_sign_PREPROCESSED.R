
#***********************************************************************
# This scripts loads the results of loadAllModelsPlaneNW2.R,
# which is the concatenation of all networks solutions within 5% reltol for 
# the runs that did not crush out of 10 attempts per patient
# 
# by Marti Bernardo-Faura.
# 
# 
# 
# adapted by Melanie Rinas, 2019
#***********************************************************************

# Use relative paths instead of absolut paths for the files
# in Rstudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd('~/Documents/combiMS/code/similarity/')


source("../utils/SMC_of_2networks.R")    



files_Similarity_folder = "../../files/similarity"



# load median models

load(paste("../../files/median_models/allMedianModels.RData",sep=""))        # called allMedianNetworks            




# compute pairwise similarity 
similarityMatrix = sapply(unname(row.names(allMedianNetworks)), function(x){sapply(unname(row.names(allMedianNetworks)), function(y){return(SMC_of_2networks(allMedianNetworks[x,], allMedianNetworks[y,]))})})  
diag(similarityMatrix) = NA

# save similarity Matrix
save(similarityMatrix, file=file.path(files_Similarity_folder,"similarityMatrix.RData"))            




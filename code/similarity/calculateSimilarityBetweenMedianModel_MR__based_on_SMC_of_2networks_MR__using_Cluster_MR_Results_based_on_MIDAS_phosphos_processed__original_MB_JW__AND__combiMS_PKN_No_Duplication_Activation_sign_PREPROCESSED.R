

# calculateSimilarityBetweenMedianModel_MR__based_on_SMC_of_2networks_MR__using_Cluster_MR_Results_based_on_MIDAS_phosphos_processed__original_MB_JW__AND__combiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED


Results_storage_name = "Results_based_on_SIF_combiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED"        # MR inserted
Results_storage_name__median_models = "median_models__based_on__SIFcombiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED"    # MR inserted


# Modifications by MR
# 
#  using code SMC_of_2networks_MR.R to calculate
# SMC = simple matching coefficient
# 
# 





# adapted by
# Melanie Rinas, December 2017
# 
#
# 
# 
# from the
# original script called 
# 
# calculateSimilarityBetweenMedianModel.R
# 
# which was established by
# Marti Bernardo-Faura.



#***********************************************************************
# This scripts loads the results of loadAllModelsPlaneNW2.R,
# which is the concatenation of all networks solutions within 5% reltol for 
# the runs that did not crush out of 10 attempts per patient
#***********************************************************************

# Use relative paths instead of absolut paths for the files
# in Rstudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd('~/Documents/combiMS/code/similarity/')

#source("./newJaccard.R")              # MR modified
source("../utils/SMC_of_2networks_MR.R")     # MR modified



files_Similarity_MR_folder = "../../files/similarity"

ifelse(!dir.exists(file.path(files_Similarity_MR_folder,Results_storage_name)), dir.create(file.path(files_Similarity_MR_folder,Results_storage_name)), FALSE)
sub_dir__Results__current_script = file.path(files_Similarity_MR_folder,Results_storage_name)









# load median models
#load('../../files/allMedianModels.RData')   # MR modified

load(paste("../../files/median_models/", Results_storage_name__median_models ,"/allMedianModels.RData",sep=""))  # called allMedianNetworks               # MR modified









# compute pairwise similarity 
similarityMatrix = sapply(unname(row.names(allMedianNetworks)), function(x){sapply(unname(row.names(allMedianNetworks)), function(y){return(SMC_of_2networks_MR(allMedianNetworks[x,], allMedianNetworks[y,]))})})  # MR modified
diag(similarityMatrix) = NA

# save similarity Matrix
#save(similarityMatrix, file='../../files/similarityMatrix.RData')                                           # MR modified
save(similarityMatrix, file=file.path(sub_dir__Results__current_script,"similarityMatrix.RData"))            # MR modified




# Script to iteratively optimise all patients.
# Marti Bernardo-Faura, September 2014

# Minor changes to adjust the paths to the combiMS Github project
# Jakob Wirbel, June 2017

# ****************** Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


source("./OptCombiMS.R")
source("../data_processing_and_normalization/normaliseSimp.R")
source("../data_processing_and_normalization/processCombiData.R")


#filenames=list.files(".",pattern="*.csv",full.names=FALSE)

#filenames=list.files("/Users/marti/Documents/ebi/combiMS/data/midas_clean_merged/processedForCNO/",pattern="*.csv",full.names=FALSE)
filenames=list.files("../../data/phosphos_merged/",pattern="*.csv",full.names=FALSE)
# setwd("/Users/marti/Documents/R/combiMS/modelling/final")



# ****************** process all files to remove controls and be consistent with PKN and midas-like formattted
# do only once :)
#   
# for (i in 1:length(filenames)){
#   processCombiData(filenames[i]) 
# }

# ****************** normalise all files


for (i in 1:length(filenames)){
  fileName=filenames[i]
  #fileName=paste(head(strsplit(filenames[i],"\\_")[[1]],n=1),".csv",sep="") #take only filename
  taula=read.csv(paste("../../data/phosphos_merged/",fileName,sep=""),header=TRUE,dec=".",check.names=FALSE, stringsAsFactors=FALSE)
  taula=normaliseSimp(taula)
  # setwd("/Users/marti/Documents/ebi/combiMS/data/phosphosMergedAbsMax/processed/normalized/")
  write.table(taula,file=paste0('../../data/phosphos_normalised/', fileName),sep=",",row.names=FALSE,quote=F)
}
  
# ****************** use data from all patients to optimize one model per patient
networks=list()
for (i in 1:length(filenames)){
  
  networks[[i]]=OptCombiMS(filenames[i])
}

# for (i in 1:6){
#   
#   networks[[i]]=OptCombiMS(filenames[1])
# }



save(networks,file="networks.Rdata")
#write.table(networks,file="/Users/marti/Documents/R/combiMS/modelling/networksIC5001.csv",sep=",",row.names=FALSE,quote=F)

#networks=read.csv("./networksIC5006.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)

# let's create the annotation labels for a nice aheatmap


# let's get from the optimization results only the average networks for each patient
#averageNetworks=lapply(networks,function(x){x$averageNw}) #not good, creates a list of lists
averageNetworks=do.call(rbind, lapply(networks,function(x){x$averageNw}))


patientAnnot=vector()
for (j in 1:length(filenames)){  
  patientAnnot[j]=strsplit(strsplit(filenames[j],"/")[[1]][10],"_")[[1]][1]
}

annot=read.csv("../../files/annotations/annotations_169_patients.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
# find the indexes of the intersection between the files we have and teo's metadata
# choose only the metadata for the patients that are measured
annot = annot[match(patientAnnot,annot$ID),]

#aheatmap(networks, distfun="euclidean", annRow=annot[,c("Group","Center","Gender")])
#use euclidean distance, although binary works for 0 and 1. Same clustering but smaller distances
#the values are shown in grey scale. The annotations in non continuosus scale so that groups next to each other are distinguishable
aheatmap(averageNetworks, annColors=c("Paired","Dark2","rainbow"),color="grey", Rowv=F, Colv=F,distfun="euclidean", labRow=annot$ID, labCol=model$reacID,annRow=annot[,c("Group","Center","Gender")])
#plot all data




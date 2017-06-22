# Script to iteratively optimise all patients.
# Marti Bernardo-Faura, September 2014

# Minor changes to adjust the paths to the combiMS Github project
# Removed processing and normalization part --> moved to normalization_pipeline.R
# Jakob Wirbel, June 2017

# ****************** Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("./OptCombiMS.R")

filenames=list.files("../../data/phosphos_processed/",pattern="*.csv",full.names=FALSE)

# ****************** use data from all patients to optimize one model per patient
networks=list()
for (i in 1:length(filenames)){
  
  networks[[i]]=OptCombiMS(filenames[i])
}

save(networks,file="../../files/modeling/local/networks.Rdata")

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




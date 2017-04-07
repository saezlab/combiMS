
# Script to create a csv file that can be imported in Cytoscape to visualize the network model for each one of the following groups:
# Healthy','MS','PPMS','RRMS','Gilenya','IFNb','Copaxone','EGCG','Tysabri'

# Marti Bernardo-Faura, November 2015.

# ***transform network file into table of cytoscape properties to annote cytoscape plot
# **************************************************
# ************which phenotype do we want to map on top of model?
# **************************************************
thisMode='mean'
phenotypeNwsFigure_folder='/Users/marti/Documents/R/combiMS/phenotypeNws/'
phenotypes=c('Healthy','MS','PPMS','RRMS','Gilenya','IFNb','Copaxone','EGCG','Tysabri')
i=1
thisPhenotype=phenotypes[i]
interactionWeights=read.table(paste0(phenotypeNwsFigure_folder,thisPhenotype,thisMode,'.csv'),sep=',')



# **************************************************
# ************load anotation to map patients to groups
# **************************************************
data_folder="/Users/marti/Documents/ebi/combiMS/data/phosphosMergedAbsMax/processed/normalized/secondRoundProcessedMidas/"
filenames=list.files(data_folder,pattern="*.csv",full.names=FALSE)
annot=read.csv("/Users/marti/Documents/R/combiMS/modelling/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
numPat=length(filenames)
filenames2=filenames
for (j in 1:numPat){
  filenames2[j]=strsplit(filenames[j],"\\.")[[1]][1]
}
#simplify category
annot2=annot
annot2$Category[which(annot$Category=="dont know")]='PPMS'
annot2$Category[which(annot2$Category=='PPMS slow')]='PPMS'
annot2$Category[which(annot2$Category=='PPMS fast')]='PPMS'

#test if all 169 are correctly labeled
length(which(annot2$Category=='PPMS')) + length(which(annot2$Category=='RRMS')) + length(which(annot2$Category=='healthy'))
length(which(annot2$Category=='RRMS' & annot2$condition=='Untreated')) + length(which(annot2$Category=='PPMS' & annot2$condition=='Untreated')) + length(which(annot2$condition=='Treated')) + length(which(annot2$condition=='Healthy'))
# *************************************************
# ************load model and midas for annotation
# **************************************************
patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
full_model_path='/Users/marti/Documents/R/combiMS/combiMSplane.sif'
model_path='/Users/marti/Documents/R/combiMS/combiMSplaneCUT.sif'
fileName=patientData[1]
midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)
full_model=readSIF(full_model_path)
numInteractions=length(model$reacID)
sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),numInteractions)



# **************************************************
# ************rewrite act and inh in cytoscape format
# **************************************************
#new table to allocate rewriten ints and network properties
interactionWeightsRewriten=interactionWeights

#separate inhibitory and activatory
inhibitoryInts=as.matrix(interactionWeights[grep('!',rownames(interactionWeights)),])
rownames(inhibitoryInts)=rownames(interactionWeights)[grep('!',rownames(interactionWeights))]
colnames(inhibitoryInts)=thisPhenotype
activatoryInts=as.matrix(interactionWeights[-grep('!',rownames(interactionWeights)),])
rownames(activatoryInts)=rownames(interactionWeights)[-grep('!',rownames(interactionWeights))]
colnames(activatoryInts)=thisPhenotype

#in inhibition, replace ! and equal
startNode=numeric(length=dim(inhibitoryInts)[1])
endNode=numeric(length=dim(inhibitoryInts)[1])
inhRewriten=numeric(length=dim(inhibitoryInts)[1])
for (j in 1:dim(inhibitoryInts)[1]){
  startNode[j]=strsplit(rownames(inhibitoryInts)[j],"\\=")[[1]][1]
  startNode[j]=strsplit(startNode[j],"\\!")[[1]][2]
}
for (j in 1:dim(inhibitoryInts)[1]){
  endNode[j]=strsplit(rownames(inhibitoryInts)[j],"\\=")[[1]][2]
}

for(j in 1:dim(inhibitoryInts)[1]){inhRewriten[j]=paste0(startNode[j],"\t","-1","\t",endNode[j])}
#add source and target to table
interactionWeightsRewriten[grep('!',rownames(interactionWeights)),2]=startNode
interactionWeightsRewriten[grep('!',rownames(interactionWeights)),3]=endNode

#in activation, replace equal
startNode=numeric(length=dim(activatoryInts)[1])
endNode=numeric(length=dim(activatoryInts)[1])
actRewriten=numeric(length=dim(activatoryInts)[1])
for (j in 1:dim(activatoryInts)[1]){
  startNode[j]=strsplit(rownames(activatoryInts)[j],"\\=")[[1]][1]
}
for (j in 1:dim(activatoryInts)[1]){
  endNode[j]=strsplit(rownames(activatoryInts)[j],"\\=")[[1]][2]
}

for(j in 1:dim(activatoryInts)[1]){actRewriten[j]=paste0(startNode[j],"\t","1","\t",endNode[j])}
#add source and target to table
interactionWeightsRewriten[-grep('!',rownames(interactionWeights)),2]=startNode
interactionWeightsRewriten[-grep('!',rownames(interactionWeights)),3]=endNode

# **************************************************
#********replace rewriten reactions
# **************************************************
rownames(interactionWeightsRewriten)[grep('!',rownames(interactionWeights))]=inhRewriten
rownames(interactionWeightsRewriten)[-grep('!',rownames(interactionWeights))]=actRewriten


# **************************************************
#********add network properties
# **************************************************
#add interaction id
interactionWeightsRewriten[1:numInteractions,4]=seq(1:numInteractions)
rownames(interactionWeightsRewriten)[grep('!',rownames(interactionWeights))]=inhRewriten
rownames(interactionWeightsRewriten)[-grep('!',rownames(interactionWeights))]=actRewriten
colnames(interactionWeightsRewriten)=c('weight','start','end','id')

#add interaction type
interactionWeightsRewriten[grep('!',rownames(interactionWeights)),5]=-1
interactionWeightsRewriten[-grep('!',rownames(interactionWeights)),5]=1
colnames(interactionWeightsRewriten)=c('weight','start','end','id','type')


#if phenotype is drug, add defective interactions
drugs=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
if(thisPhenotype %in% drugs){
drugScores_folder='/Users/marti/Documents/R/combiMS/drugScores/'
fileName=paste0(drugScores_folder,thisPhenotype,'median',".csv")
drugScores=read.csv(fileName,header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
interactionWeightsRewriten[which(rownames(interactionWeights) %in% rownames(drugScores)),6]='yes'
interactionWeightsRewriten[-which(rownames(interactionWeights) %in% rownames(drugScores)),6]='no'
colnames(interactionWeightsRewriten)=c('weight','start','end','id','type','defective')
}

# **************************************************
# ************are there defective ints that won't appear because of being inactive?
# **************************************************
inactiveInts=which(interactionWeightsRewriten$weight==0)
defectiveInts=which(interactionWeightsRewriten$defective=='yes')
warning(length(which(defectiveInts %in% inactiveInts)),' of ',length(defectiveInts),' defective ints do not appear\n')

# **************************************************
# ************cut model (if we are in mode median)
# **************************************************
if(thisMode=='median'){
interactionWeightsRewriten=interactionWeightsRewriten[-which(interactionWeightsRewriten$weight==0),]
}
# **************************************************
#********save
# **************************************************
write.csv(interactionWeightsRewriten,paste0(phenotypeNwsFigure_folder,thisPhenotype,'Weights','Median','.csv'),quote=F)

# **************************************************
#********create a table of node properties
# **************************************************


# **************************************************
#********create a table of node properties
# **************************************************
# #create a node table specifying if readout is signal or stimulus
# SignalYes=matrix(nrow=length(model$namesSpecies),ncol=2)
# SignalYes[,1]=model$namesSpecies
# SignalYes[which(model$namesSpecies %in% colnames(midas@signals[[2]])),2]='signal'
# SignalYes[-which(model$namesSpecies %in% colnames(midas@signals[[2]])),2]='none'
# SignalYes[which(model$namesSpecies %in% colnames(midas@stimuli)),2]='stimulus'
# 
# write.csv(SignalYes,paste0(phenotypeNwsFigure_folder,'SignalStimulus','.csv'),quote=F)
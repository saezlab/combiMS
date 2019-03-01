# Function called by pathDrugTargetsFinalv3.R to path search and predict experiments coupling a stimulus with a measureable readout
# Specifically, only those paths are selected that include a co-druggable interaction

# Final version April 2016
# Marti Bernardo-Faura

# minor changed for github repository by Jakob Wirbel, July 2017
# 
# minor changed for github repository by Melanie Rinas, 2019

# *************************************************************************************************************************
# ***********determine if a stimulus connects to signals through predicted interactions
# *************************************************************************************************************************



makeMatrixPathStimSignal<-function(graphGroup,link){
  
  source("./is.connected.R")
  
  #graphGroup=subnetwork$graph
  cat("the interaction is",link,"\n")
  
  # ************load model and midas for annotation
  patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
  model_path="../../files/model/combiMS_PKN_No_Duplication_Activation_sign.sif"
  fileName=patientData[1]
  midas=CNOlist(paste(data_folder,fileName,sep=""))
  model=readSIF(model_path)  
  sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
  model=preprocessing(midas,model,expansion=FALSE)
  numInteractions=length(model$reacID)
  sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
  numInteractions=length(model$reacID)
  
  
  # find out nodes & stimulus kept after cutting interactions not found active in group union network
  presentStimuli=colnames(midas@stimuli)[which(colnames(midas@stimuli) %in% graphGroup@nodes)]
  presentReadouts=colnames(midas@signals[[2]])[which(colnames(midas@signals[[2]]) %in% graphGroup@nodes)]
  
  
  # load generic vars to generate matrix of paths
  AllStimulus=colnames(midas@stimuli)
  numStimulus=length(AllStimulus)
  AllReadouts=colnames(midas@signals[[2]])
  numReadouts=length(colnames(midas@signals[[2]]))
  

  pathFound=matrix(nrow=numReadouts,ncol=numStimulus)
  for (i in 1:numReadouts){
    readout=AllReadouts[i]
    cat("the readout is",readout,"\n")
   
    for(j in 1:numStimulus){
      stimulus=AllStimulus[j]
      cat("****stimulus",stimulus,"\n")
      if(stimulus %in% presentStimuli & readout %in% presentReadouts){
        pathYes=is.connected(graph=graphGroup,stimulus=stimulus,readout=readout,link=link)
        pathFound[i,j]=pathYes
      } else{cat("stimulus or readout not present in group graph","\n")}
    }
  }
  colnames(pathFound)=AllStimulus
  rownames(pathFound)=AllReadouts
  pathFound[is.na(pathFound)]=0
    return(pathFound)
}

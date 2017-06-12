# Function to remove signals that were only measured to allow experimental QC and dynamic range calculation (PR; BSA). Also sanity check anc donsistency.
# This is called from command line and OptAllPatients.R
# Marti Bernardo Faura, Final version September 2014

processCombiData=function(fileName){
  
  # instead of directly loading into midas, first clean the headers, which for the TR: have a e.g. :cytokine field added for each column
  #taula=read.csv(paste("/Users/marti/Documents/ebi/combiMS/data/phosphos/midas/",fileName,sep=""),header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
  taula=read.csv(paste("/Users/marti/Documents/ebi/combiMS/data/phosphosMergedAbsMax/",fileName,sep=""),header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
  taula=taula[,-1] # remove the first column, which are row labels
  taula=taula[-2,] # remove the first row, which is repeated (control at 0 and at 5 is the same)
  #pvataula[is.na(taula)]=0
  prova=vector()
  # first remove cytokine labels
  for (i in 1:length(names(taula))){  
    prova[i]=paste(strsplit(names(taula)[i],"\\.")[[1]][1],strsplit(names(taula)[i],"\\.")[[1]][2],sep=":")
  }
  
  #the second position is the cell line and by splitting it we cut CellLine from it, lets restore
  #prova[1]=paste(prova[1],"CellLine",sep=":")
  prova[1]=paste(prova[1],":CellType",sep="")
  #finally save same table with reformated headers
  names(taula)=prova
  #taula=taula[,-which(colnames(taula) == "TR:phophos")]
  taula=taula[,-which(colnames(taula) == "TR:phospho")]
  taula=taula[,-which(colnames(taula) == "TR:medium")]
  #taula=taula[-c(2,3),] #here we remove the two fake columns that repeated the control
  names(taula)[names(taula)=="TR:Teri"]="TR:Teriflunomide"
  taula=taula[,-which(colnames(taula) == "DA:BSA")] #these four lines remove xMAP controls
  taula=taula[,-which(colnames(taula) == "DV:BSA")]
  taula=taula[,-which(colnames(taula) == "DA:PE")]
  taula=taula[,-which(colnames(taula) == "DV:PE")]
  
  # *************************** plot heatmap
  #get treatment names to label heatmap  
  #   namesStimuli=colnames(taula[,grep('^TR',colnames(taula))])
  #   namesStimuli=namesStimuli[-1] # remove cell name
  #   for(i in 1:length(namesStimuli)){
  #     namesStimuli[i]=strsplit(namesStimuli,"\\:")[[i]][2]
  #   }
  #   namesStimuli=append("medium",namesStimuli)
  #   
  #   #remove DV to label heatmap
  #   namesSignals=colnames(taula[,grep('^DV',colnames(taula))])
  #   for(i in 1:length(namesSignals)){
  #     namesSignals[i]=strsplit(namesSignals,"\\:")[[i]][2]
  #   }
  
  #aheatmap(taula[,grep('^DV:', colnames(taula))],labRow=namesStimuli,labCol=namesSignals, Rowv=NA, Colv=NA, scale="column")
  #   aheatmap(taula[,grep('^DV:', colnames(taula))],labRow=namesStimuli,labCol=namesSignals, Rowv=NA, Colv=NA)
  
  fileNameClean=head(strsplit(fileName,"\\_")[[1]],n=1) #take only patient name
  setwd("/Users/marti/Documents/ebi/combiMS/data/phosphosMergedAbsMax/normalized/")
  write.table(taula,file=paste(fileNameClean,"csv",sep="."),sep=",",row.names=FALSE,quote=F)
  
  
}
# Local version of the script Run in cluster for patient optimisation. Called by OptAllPatients.R
# Marti Bernardo Faura, September 2014. 

# Minor changes to adjust the paths to the combiMS Github project
# Jakob Wirbel, June 2017

OptCombiMS=function(fileName){
  library(CellNOptR)
  library(NMF)
  # setwd("/Users/marti/Documents/R/combiMS/modelling/final")
  # source("../data_processing_and_normalization/normaliseSimp.R")
  
  
  # *************************************************
  # ************find patient data, load its midas
  # **************************************************
  data_folder="../../data/phosphos_processed/"
  patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
  i=1; fileName=patientData[i] #give a specific patient if not passed via function
  cat(sprintf("*********The patient is %s",fileName))
  # *************************** load processed midas
  # if already processed (by calling processCombiData.R in OptAllPatients.R), load the processed file
  #taula=read.csv(paste("/Users/marti/Documents/ebi/combiMS/data/midas_clean_merged/processedForCNO/normalized/",fileName,sep=""),header=TRUE,dec=".",check.names=FALSE, stringsAsFactors=FALSE)
  midas=CNOlist(paste(data_folder,fileName,sep=""))
  
  # *************************************************
  # ************load compressed model
  # **************************************************
  model_path='../../files/model/combiMSplaneCUT.sif'
  #full_model_path='/Users/marti/Documents/R/combiMS/combiMSplane.sif'
  model=readSIF(model_path)
  numInteractions=length(model$reacID)
  sprintf("*********After compressing without expanding logic gates, the model has %s nodes and %s reactions",length(model$namesSpecies),numInteractions)
  #plotModel(model,midas,graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))
  
  # *************************************************
  # ************normalise and plot
  # **************************************************
  #taulaHill2=normaliseSimp(taula,HillCoef=2)
  #write.table(taulaHill2,file="/Users/marti/Documents/R/combiMS/modelling/midas.csv",sep=",",row.names=FALSE,quote=F)
  #midas=CNOlist("/Users/marti/Documents/R/combiMS/modelling/midas.csv")
  # plot(midas)
  #midas@signals[[2]][,grep('^DIG1',colnames(midas@signals[[1]]))]

  
  
  # *************************************************
  # ************fit model
  # **************************************************
  initBstring=rep(1,length(model$reacID))
  cutAndPlot(midas,model, list(initBstring)) #simulate
  Opt=gaBinaryT1(CNOlist=midas,model=model,initBstring=initBstring,stallGenMax=600,maxTime=20, maxGens=100000, verbose=FALSE,popSize=100,elitism=2)
#Opt=gaBinaryLargeNw(CNOlist=midas,model=model,initBstring=initBstring,stallGenMax=600,maxTime=20, maxGens=100000, verbose=FALSE,popSize=100,elitism=2)
  
  
#***************************plot result
  #plot fits
  prova=cutAndPlot(model=model,bStrings=list(Opt$bString),CNOlist=midas, plotPDF=FALSE, plotParams = list(maxrow = 25,cex=.8))
  # plot objective function
  plotFit(optRes=Opt)
  #plot network
  #plotModel(model,midas,Opt$bString)
  bestNw=Opt$bString
  
  # *************************** fit model with second timepoint (starting bitString from previous opt)
  #   OptT2=gaBinaryT2(CNOlist=midas,model=model,bStringT1=Opt$bString,verbose=TRUE)
  #   
  #   #plot result
  #   #plot fits
  #   cutAndPlot(model=model,bStrings=list(OptT2$bString),CNOlist=midas, plotPDF=TRUE, plotParams = list(maxrow = 25,cex=.8))
  #   # plot objective function
  #   plotFit(optRes=OptT2)
  #   #plot network
  #   interactionsT1T2=buildBitString(list(Opt$bString,OptT2$bString))
  #   plotModel(model,midas,bString=interactionsT1T2$bs,indexIntegr=interactionsT1T2$bsTimes)
  #   
  # *************************** calculate average times an interaction survives in the all networks within relTol
  #there is a best network expressed as interaction vector and there is the same for all networks within a relTol
  sprintf("*********the number of networks within relTols is %s",dim(Opt$stringsTol)[1])
  countInteractions=numeric(length=dim(Opt$stringsTol)[2])
  for(i in 1:dim(Opt$stringsTol)[2]){
    countInteractions[i]=sum(Opt$stringsTol[,i])
  }
  averageNw=countInteractions/dim(Opt$stringsTol)[1]*100
  
  # *************************** calculate RMSE for each signal across all experiments
  
  #midas@signals contains the value at 0 and 5, hence we take the second
  #then we take the second again, since rows are stimuli and columns signals
  numSignals=dim(midas@signals[[2]])[2]
  numStimuli=dim(midas@signals[[2]])[1]
  rmseAll=numeric(length=numSignals) 
  
  for(i in 1:numSignals){
    rmseAll[i]=sqrt((sum(midas@signals[[2]][,i]-prova$simResults[[1]]$t1[,i])^2)/numStimuli)
  }
  
  
  c=list("averageNw"=averageNw,"rmseAll"=rmseAll,"countInteractions"=countInteractions,"countNws"=dim(Opt$stringsTol)[1],"Opt"=Opt)
  
  return(c)
  
}

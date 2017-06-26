# Final version of the optimization algorithm for a single patient. All parameters applied to final version
# Marti Bernardo-Faura, Juni 2015

# Minor changes to adjust the paths to the combiMS Github project
# Jakob Wirbel, June 2017

OptCombiMScluster=function(fileName){
  
  library(CellNOptR)
  #library(NMF)
  #setwd("/Users/marti/Documents/R/combiMS/modelling/final")
  
  #source("/Users/marti/Documents/R/combiMS/normaliseSimp.R")
  #fileName="/Users/marti/Documents/ebi/combiMS/data/midas_clean_merged/CH003_phosphos_midas.csv"
  cat(sprintf("*********The patient is %s",fileName))
  
  # *************************** load model
  model=readSIF("../../files/model/combiMSplane.sif")
  
  # *************************** load processed midas
  # if already processed (by calling processCombiData.R in OptAllPatients.R), load the processed file
  #taula=read.csv(paste("/Users/marti/Documents/ebi/combiMS/data/midas_clean_merged/processedForCNO/normalized/",fileName,sep=""),header=TRUE,dec=".",check.names=FALSE, stringsAsFactors=FALSE)
  
  
  # *************************** normalise and plot
  #taulaHill2=normaliseSimp(taula,HillCoef=2)
  #write.table(taulaHill2,file="/Users/marti/Documents/R/combiMS/modelling/midas.csv",sep=",",row.names=FALSE,quote=F)
  #midas=CNOlist("/Users/marti/Documents/R/combiMS/modelling/midas.csv")
  midas=CNOlist(paste("../../data/phosphos_processed/",fileName,sep=""))
  # plot(midas)
  #midas@signals[[2]][,grep('^DIG1',colnames(midas@signals[[1]]))]
  
  
  # *************************** format model
  # compress and expand model
  cat(sprintf("*********Before cluster processing, the model has %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID)))
  model=preprocessing(midas,model,expansion=F,compression=TRUE)
  cat(sprintf("*********After cluster compressing and expanding logic gates, the model has %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID)))
  #plotModel(model,midas)
  
  
  # *************************** fit model
  initBstring=rep(1,length(model$reacID))
  #cutAndPlot(midas,model, list(initBstring)) #simulate
  cat("********Optimizing",sep="\n")
  #Opt=gaBinaryLargeNw(CNOlist=midas,model=model,initBstring=initBstring,stallGenMax=600000,maxTime=36000, maxGens=100000, verbose=FALSE,popSize=100,elitism=2,relTol=0.05)
  Opt=gaBinaryT1(CNOlist=midas,
                 model=model,
                 initBstring=initBstring,
                 stallGenMax=600000,
                 maxTime=36000, 
                 maxGens=100000, 
                 verbose=FALSE,
                 popSize=100,
                 elitism=2,
                 relTol=0.05)
  
  
  
  #***************************plot result
  #plot fits
  sim=cutAndPlot(model=model,bStrings=list(Opt$bString),CNOlist=midas, plotPDF=FALSE, plotParams = list(maxrow = 25,cex=.8))
  
  cat("********Simulation completed",sep="\n")
  # plot objective function
  #plotFit(optRes=Opt)
  #plot network
  #plotModel(model,midas,Opt$bString)
  bestNw=Opt$bString
  #print(Opt$stringsTol)
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
  numNetworks=dim(Opt$stringsTol)[1]
  cat(sprintf("*********the number of networks within relTols is %s",numNetworks),sep="\n")
  numEdges=dim(Opt$stringsTol)[2]
  cat(sprintf("*******the number of edges in each network is %s",numEdges),sep="\n")  
  countInteractions=numeric(length=numEdges)
  for(i in 1:numEdges){
    countInteractions[i]=sum(Opt$stringsTol[,i])
  }
  averageNw=countInteractions/numNetworks*100
  
  # *************************** calculate RMSE for each signal across all experiments
  
  #midas@signals contains the value at 0 and 5, hence we take the second
  #then we take the second again, since rows are stimuli and columns signals
  numSignals=dim(midas@signals[[2]])[2]
  numStimuli=dim(midas@signals[[2]])[1]
  rmseAll=numeric(length=numSignals) 
  
  for(i in 1:numSignals){
    rmseAll[i]=sqrt((sum(midas@signals[[2]][,i]-sim$simResults[[1]]$t1[,i])^2)/numStimuli)
  }
  
  
  fileName2=unlist(strsplit(fileName,"\\."))[1]
  c=list("patientName"=fileName2,"averageNw"=averageNw,"rmseAll"=rmseAll,"countInteractions"=countInteractions,"countNws"=dim(Opt$stringsTol)[1],"Opt"=Opt)
  
  return(c)
  
}

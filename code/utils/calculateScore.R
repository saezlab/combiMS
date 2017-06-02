

calculateScore <- function (model,midas,bString){
#   if(!exists("bString")){
#     bString=rep(1,length(model$reacID))
#   }
  sim=cutAndPlot(model=model,bStrings=list(bString),CNOlist=midas, plotPDF=FALSE, plotParams=list(maxrow=50))
  
  #************* calculate deviation
  timeIndex=2
  simResultsT0 = sim$simResults[[1]]$t0
  simResults = sim$simResults[[1]]$t1
  Diff0 <- simResultsT0 - midas@signals[[1]]
  Diff <- simResults - midas@signals[[timeIndex]]
  r0 <- Diff0^2
  r <- Diff^2
  r <- rbind(r0, r)
  deviationPen <- sum(r[!is.na(r)])/2
  
  #************* calculate penalty size
  sizeFac = 1e-04
  #nInTot is the number of inputs prior to optimization after preprocessing
  nInTot <- length(which(model$interMat == -1))
  #nInputs <- length(which(model$interMat == -1))     
  nInputs=length(which(bString==1))
  nDataPts <- dim(midas@signals[[timeIndex]])[1] * dim(midas@signals[[timeIndex]])[2]
  sizePen <- (nDataPts * sizeFac * nInputs)/nInTot
  
  #************* calculate NAPen
  NAFac = 1
  NAPen <- NAFac * length(which(is.na(simResults)))
  
  #************* sum all scores to produce final score
  # achtung: the total after computeScoreT1 is still divideb by numMeasurements in getFit
  numMeasurements=sum(!is.na(midas@signals[[timeIndex]]))
  score = (deviationPen + NAPen + sizePen) / numMeasurements
  scores=list("score"=score,"deviationPen"=deviationPen,"NAPen"=NAPen,"sizePen"=sizePen,"numMeasurements"=numMeasurements)
  return(scores)  
  
}


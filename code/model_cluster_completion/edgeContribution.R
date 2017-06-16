#***********************************************************************
# This scripts sanity checks edge contribution for each model
# Juny 2015
# Marti Bernardo-Faura
#***********************************************************************


edgeContribution <- function(Opt, model, midas){
  

  # ************************************
  # *************************** given Opt, assess edges kept, and how many are ANDs
  # ************************************
  labelsEdges=model$reacID[which(Opt$bString==1)]
  keptEdges=which(Opt$bString==1)
  scoresTest=numeric(length=length(keptEdges))
  cat("----------edges kept:",length(keptEdges),"\n")
  cat("----------number of AND gates in optimized",length(grep('\\+',model$reacID[keptEdges])),"\n")
  
  
  # ************************************
  # *************************** calculate score in full reference model
  # ************************************
  scoreReference=calculateScore(model,midas,Opt$bString)
  
  # ************************************
  # *************************** Remove one edge at a time and calculate simulation score
  # ************************************
  cat("--------interaction removed","\t")
  for (edgeIndex in 1:length(keptEdges)){
    testString=Opt$bString
    testString[keptEdges[edgeIndex]]=0
    cat(edgeIndex,labelsEdges[edgeIndex],"\t")
    scoresTest[edgeIndex]=calculateScore(model,midas,testString)
  }
  
  # ************************************
  # *************************** Remove reference from each score test for better illustration
  # ************************************
  scoresDifference=scoresTest-scoreReference
  
  # ************************************
  # visualize results
  # ************************************
  #create vector ANDgate to visualize if an interaction is with an AND
  ANDgate=rep("no",length(keptEdges))
  ANDgate[grep('\\+',model$reacID[keptEdges])]="yes"
  scoresDF=data.frame(cbind(as.double(scoresDifference),labelsEdges,ANDgate))
  plotErrors=ggplot(scoresDF, aes(labelsEdges, scoresDifference,fill=ANDgate))+geom_bar(stat="identity")+coord_flip()
  
  # ************************************
  # *************************** return edge removed and error of model simulation upon removal
  # ************************************
  c=list("edgesRemoved"=labelsEdges,"errorModels"=scoresDifference,"plot"=plotErrors)
  return(c)  
}
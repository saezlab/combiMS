
# calculateScore_MR
# differs from calculateScore
# by
# (1)  comments for better understanding
# 
# (2)  and modified cutAndPlot surrounding,
#      to prevent that an error occurs as a result of too many open graphical devices:
      
       # graphics.off()                                  # MR modified

       # pdf(file = NULL)                                # MR modified
       # sim=cutAndPlot(model=model,bStrings=list(bString),CNOlist=midas, plotPDF=FALSE, plotParams=list(maxrow=50))  # MR modified
       # dev.off                                         # MR modified



# START MR inserted
# 

# comparison with the settings used for fitting within the script
# OptCombiMSclusterFork10_MR__using_MIDAS_phosphos_processed__original_MB_JW__AND__optimization_fct__gaBinaryT1__with_set_seed_number.R


# sim=cutAndPlot(model=model,bStrings=list(Opt$bString),CNOlist=midas, plotPDF=FALSE, plotParams = list(maxrow = 25,cex=.8))


# for(i in 1:numSignals){
#    rmseAll_MR[i]=sqrt(    (sum( (midas@signals[[2]][,i] - sim$simResults[[1]]$t1[,i] )^2 ) /numStimuli)   )
# }



# gaBinaryT1__with_set_seed_number =  function (CNOlist, model, initBstring = NULL, sizeFac = 1e-04, 
#                                               NAFac = 1, popSize = 50, pMutation = 0.5, maxTime = 60, maxGens = 500, 
#                                               stallGenMax = 100, selPress = 1.2, elitism = 5, relTol = 0.1, 
#                                               verbose = TRUE, priorBitString = NULL, timeIndex = 2,
#                                               seed_number = NULL) 


# Opt=gaBinaryT1__with_set_seed_number(                    # MRmodified
#    CNOlist=midas,                                
#    model=model,
#    initBstring=initBstring,
#    stallGenMax=600000,
#    maxTime=36000, 
#    maxGens=100000, 
#    verbose=FALSE,
#    popSize=100,
#    elitism=2,
#    relTol=0.05,
#    seed_number = seed_number_used)            # MRmodified





# 
# 
# gaBinaryT1__with_set_seed_number =  function (CNOlist, model, initBstring = NULL, sizeFac = 1e-04, 
#                                               NAFac = 1, popSize = 50, pMutation = 0.5, maxTime = 60, maxGens = 500, 
#                                               stallGenMax = 100, selPress = 1.2, elitism = 5, relTol = 0.1, 
#                                               verbose = TRUE, priorBitString = NULL, timeIndex = 2,
#                                               seed_number = NULL) 
# {
#    
#    
#    if (is.null(seed_number) == TRUE) {
#       current_time_as_number = as.numeric(format(Sys.time(), "%s"))
#       # seed_number= current_time_as_number - (generationFile_interval-ID)
#       # seed_number= current_time_as_number 
#       set.seed(current_time_as_number)
#    }else{
#       set.seed(seed_number)
#    }
#    
#    
#    
#    #print(seed_number)
#    cat("*********** current seed number = ",seed_number,sep="\n")
#    
#    
#    if (is.null(initBstring) == TRUE) {
#       initBstring <- rep(1, length(model$reacID))
#    }
#    if ((class(CNOlist) == "CNOlist") == FALSE) {
#       CNOlist = CellNOptR::CNOlist(CNOlist)
#    }
#    checkSignals(CNOlist, model)
#    if (timeIndex < 2) {
#       stop("timeIndex must be >=2")
#    }
#    if (timeIndex > length(CNOlist@timepoints)) {
#       stop(paste("timeIndex must be <= ", length(CNOlist@timepoints), 
#                  sep = " "))
#    }
#    bLength <- length(initBstring)
#    simList = prep4sim(model)
#    indexList = indexFinder(CNOlist, model)
#    Pop <- rbind(initBstring, round(matrix(runif(bLength * (popSize - 
#                                                               1)), nrow = (popSize - 1), ncol = bLength)))
#    Pop <- addPriorKnowledge(Pop, priorBitString)
#    bestbit <- Pop[1, ]
#    bestobj <- Inf
#    stop <- FALSE
#    g <- 0
#    stallGen <- 0
#    res <- rbind(c(g, bestobj, toString(bestbit), stallGen, Inf, 
#                   Inf, toString(bestbit), 0), c(g, bestobj, toString(bestbit), 
#                                                 stallGen, Inf, Inf, toString(bestbit), 0))
#    colnames(res) <- c("Generation", "Best_score", "Best_bitString", 
#                       "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", 
#                       "Best_bit_Gen", "Iter_time")
#    PopTol <- rep(NA, bLength)
#    PopTolScores <- NA
#    library(hash)
#    scores2Hash = hash()
#    getObj <- function(x) {
#       key = toString(.int2dec(x))
#       if (has.key(key, scores2Hash) == TRUE) {
#          return(scores2Hash[[key]])
#       }
#       else {
#          Score = computeScoreT1(CNOlist, model, x, simList, 
#                                 indexList, sizeFac, NAFac, timeIndex)
#          if (length(scores2Hash) < 1000) {
#             scores2Hash[[key]] = Score
#          }
#       }
#       return(Score)
#    }
#    t0 <- Sys.time()
#    t <- t0
#    if (popSize * stallGenMax > 2^bLength) {
#       print("Given your input parameter, an exhaustive search will be faster...")
#       stop = TRUE
#       res = exhaustive(CNOlist, model, relTol = relTol, sizeFac = sizeFac, 
#                        NAFac = NAFac, verbose = verbose)
#       return(res)
#    }
#    if (bLength == 1) {
#       Pop = matrix(c(0, 1), nrow = 2)
#       scores <- apply(Pop, 1, getObj, scoresHash = NULL)
#       rankP <- order(scores, decreasing = TRUE)
#       iBest = rankP[2]
#       return(list(bString = Pop[iBest, ], bScore = scores[iBest], 
#                   results = res, stringsTol = PopTol, stringsTolScores = PopTolScores))
#    }
#    while (!stop) {
#       scores <- apply(Pop, 1, getObj)
#       rankP <- order(scores, decreasing = TRUE)
#       Pop <- Pop[rankP, ]
#       scores <- scores[rankP]
#       fitness <- 2 - selPress + (2 * (selPress - 1) * (c(1:popSize) - 
#                                                           1)/(popSize - 1))
#       wheel1 <- cumsum(fitness/sum(fitness))
#       breaks <- runif(1) * 1/popSize
#       breaks <- c(breaks, breaks + ((1:(popSize - 1)))/popSize)
#       sel <- rep(1, popSize)
#       for (i in 1:length(breaks)) {
#          sel[i] <- which(wheel1 > breaks[i])[1]
#       }
#       Pop2 <- Pop[sel, ]
#       PSize2 <- dim(Pop2)[1]
#       PSize3 <- popSize - elitism
#       mates <- cbind(ceiling(runif(PSize3) * PSize2), ceiling(runif(PSize3) * 
#                                                                  PSize2))
#       InhBit <- matrix(runif((PSize3 * bLength)), nrow = PSize3, 
#                        ncol = bLength)
#       InhBit <- InhBit < 0.5
#       Pop3par1 <- Pop2[mates[, 1], ]
#       Pop3par2 <- Pop2[mates[, 2], ]
#       Pop3 <- Pop3par2
#       Pop3[InhBit] <- Pop3par1[InhBit]
#       MutProba <- matrix(runif((PSize3 * bLength)), nrow = PSize3, 
#                          ncol = bLength)
#       MutProba <- (MutProba < (pMutation/bLength))
#       Pop3[MutProba] <- 1 - Pop3[MutProba]
#       t <- c(t, Sys.time())
#       g <- g + 1
#       thisGenBest <- scores[length(scores)]
#       thisGenBestBit <- Pop[length(scores), ]
#       if (is.na(thisGenBest)) {
#          thisGenBest <- min(scores, na.rm = TRUE)
#          thisGenBestBit <- Pop[which(scores == thisGenBest)[1], 
#                                ]
#       }
#       if (thisGenBest < bestobj) {
#          bestobj <- thisGenBest
#          bestbit <- thisGenBestBit
#          stallGen <- 0
#       }
#       else {
#          stallGen <- stallGen + 1
#       }
#       resThisGen <- c(g, bestobj, toString(bestbit), stallGen, 
#                       (mean(scores, na.rm = TRUE)), thisGenBest, toString(thisGenBestBit), 
#                       as.numeric((t[length(t)] - t[length(t) - 1]), units = "secs"))
#       names(resThisGen) <- c("Generation", "Best_score", "Best_bitString", 
#                              "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", 
#                              "Best_bit_Gen", "Iter_time")
#       if (verbose) {
#          this = resThisGen
#          this[[3]] = substring(this[[3]], 1, 80)
#          this[[7]] = substring(this[[7]], 1, 80)
#          print(this)
#       }
#       res <- rbind(res, resThisGen)
#       Criteria <- c((stallGen > stallGenMax), (as.numeric((t[length(t)] - 
#                                                               t[1]), units = "secs") > maxTime), (g > maxGens))
#       if (any(Criteria)) 
#          stop <- TRUE
#       tolScore <- scores[length(scores)] * relTol
#       TolBs <- which(scores < scores[length(scores)] + tolScore)
#       if (length(TolBs) > 0) {
#          PopTol <- rbind(PopTol, Pop[TolBs, ])
#          PopTolScores <- c(PopTolScores, scores[TolBs])
#       }
#       if (elitism > 0) {
#          Pop <- rbind(Pop3, Pop[(popSize - elitism + 1):popSize, 
#                                 ])
#       }
#       else {
#          Pop <- Pop3
#       }
#       Pop <- addPriorKnowledge(Pop, priorBitString)
#    }
#    PopTol <- as.matrix(PopTol[-1, ])
#    PopTolScores <- PopTolScores[-1]
#    TolBs <- which(PopTolScores < scores[length(scores)] + tolScore)
#    PopTol <- as.matrix(PopTol[TolBs, ])
#    PopTolScores <- PopTolScores[TolBs]
#    PopTolT <- cbind(PopTol, PopTolScores)
#    PopTolT <- unique(PopTolT, MARGIN = 1)
#    if (!is.null(dim(PopTolT))) {
#       PopTol <- PopTolT[, 1:(dim(PopTolT)[2] - 1)]
#       PopTolScores <- PopTolT[, dim(PopTolT)[2]]
#    }
#    else {
#       PopTol <- PopTolT[1:(length(PopTolT) - 1)]
#       PopTolScores <- PopTolT[length(PopTolT)]
#    }
#    res <- res[3:dim(res)[1], ]
#    rownames(res) <- NULL
#    return(list(bString = bestbit, bScore = bestobj, results = res, 
#                stringsTol = PopTol, stringsTolScores = PopTolScores))
# }
# #<environment: namespace:CellNOptR>
# 





# 
# 
# 
# addPriorKnowledge <- function(pop, priorBitString){
#    if (is.null(priorBitString) == TRUE){
#       return(pop)
#    }
#    else{
#       for (i in 1:dim(pop)[1]){
#          pop[i,!is.na(priorBitString)] = priorBitString[!is.na(priorBitString)]
#       }
#    }
#    return(pop)
# }
# 
# 
# 
# .int2dec <- function(x){
#    return(sum(x*2^(rev(seq(x))-1)))
# }

















# 
# 
# computeScoreT1__source_code = function (CNOlist, model, bString, simList = NULL, indexList = NULL, 
#                                         sizeFac = 1e-04, NAFac = 1, timeIndex = 2) 
# {
#    if (is.null(simList) == TRUE) {
#       simList = prep4sim(model)
#    }
#    if (is.null(indexList) == TRUE) {
#       indexList = indexFinder(CNOlist, model)
#    }
#    bs = as.logical(bString)
#    modelCut <- list()
#    modelCut$interMat <- model$interMat[, bs]
#    modelCut$reacID <- model$reacID[bs]
#    modelCut$namesSpecies <- model$namesSpecies
#    simListCut <- cutSimList(simList, bString)
#    nStimuli = length(indexList$stimulated)
#    nInhibitors <- length(indexList$inhibited)
#    nCond <- dim(CNOlist@stimuli)[1]
#    nReacs <- length(modelCut$reacID)
#    nSpecies <- length(model$namesSpecies)
#    nMaxInputs <- dim(simListCut$finalCube)[2]
#    finalCube = as.integer(simListCut$finalCube - 1)
#    ixNeg = as.integer(simListCut$ixNeg)
#    ignoreCube = as.integer(simListCut$ignoreCube)
#    maxIx = as.integer(simListCut$maxIx - 1)
#    indexSignals <- as.integer(indexList$signals - 1)
#    indexStimuli <- as.integer(indexList$stimulated - 1)
#    indexInhibitors <- as.integer(indexList$inhibited - 1)
#    nSignals <- length(indexSignals)
#    valueInhibitors <- as.integer(CNOlist@inhibitors)
#    valueStimuli <- as.integer(CNOlist@stimuli)
#    simResults = .Call("simulatorT1", nStimuli, nInhibitors, 
#                       nCond, nReacs, nSpecies, nSignals, nMaxInputs, finalCube, 
#                       ixNeg, ignoreCube, maxIx, indexSignals, indexStimuli, 
#                       indexInhibitors, valueInhibitors, valueStimuli, as.integer(1))
#    simResultsT0 = .Call("simulatorT1", nStimuli, nInhibitors, 
#                         nCond, nReacs, nSpecies, nSignals, nMaxInputs, finalCube, 
#                         ixNeg, ignoreCube, maxIx, indexSignals, indexStimuli, 
#                         indexInhibitors, valueInhibitors, valueStimuli, as.integer(0))
#    nInTot = length(which(model$interMat == -1))
#    Score <- getFit(simResults = simResults, CNOlist = CNOlist, 
#                    model = modelCut, indexList = indexList, timePoint = timeIndex, 
#                    sizeFac = sizeFac, NAFac = NAFac, 
#                    nInTot = length(which(model$interMat == -1)), simResultsT0 = simResultsT0)
#    if ((class(CNOlist) == "CNOlist") == FALSE) {
#       CNOlist = CellNOptR::CNOlist(CNOlist)
#    }
#    nDataP <- sum(!is.na(CNOlist@signals[[timeIndex]]))
#    Score <- Score/nDataP
#    return(Score)
# }
# <environment: namespace:CellNOptR>
   





# 
# getFit__source_code = function (simResults, CNOlist, model, indexList = NULL, 
#                                 timePoint = c("t1","t2"), sizeFac = 1e-04, NAFac = 1, nInTot, simResultsT0 = NULL) 
# {
#    if ((class(CNOlist) == "CNOlist") == FALSE) {
#       CNOlist = CellNOptR::CNOlist(CNOlist)
#    }
#    if (is.null(indexList) == FALSE) {
#       simResults <- simResults[, indexList$signals]
#       if (is.null(simResultsT0) == FALSE) {
#          simResultsT0 <- simResultsT0[, indexList$signals]
#       }
#    }
#    if (timePoint == "t1") {
#       tPt <- 2
#    }
#    else {
#       if (timePoint == "t2") {
#          tPt <- 3
#       }
#       else {
#          tPt <- timePoint
#       }
#    }
#    if (tPt == 2 && is.null(simResultsT0) == FALSE) {
#       Diff0 <- simResultsT0 - CNOlist@signals[[1]]
#       Diff <- simResults - CNOlist@signals[[tPt]]
#       r0 <- Diff0^2
#       r <- Diff^2
#       r <- rbind(r0, r)
#       deviationPen <- sum(r[!is.na(r)])/2
#    }
#    else {
#       Diff <- simResults - CNOlist@signals[[tPt]]
#       r <- Diff^2
#       deviationPen <- sum(r[!is.na(r)])
#    }
#    NAPen <- NAFac * length(which(is.na(simResults)))
#    nDataPts <- dim(CNOlist@signals[[tPt]])[1] * dim(CNOlist@signals[[tPt]])[2]
#    nInputs <- length(which(model$interMat == -1))
#    sizePen <- (nDataPts * sizeFac * nInputs)/nInTot
#    score <- deviationPen + NAPen + sizePen
#    return(score)
# }
# <environment: namespace:CellNOptR>
   
   
   
   
   
# END MR inserted







calculateScore_MR<- function (model,midas,bString){
   #   if(!exists("bString")){
   #     bString=rep(1,length(model$reacID))
   #   }
   
   graphics.off()                                  # MR modified

   pdf(file = NULL)                                # MR modified
   sim=cutAndPlot(model=model,bStrings=list(bString),CNOlist=midas, plotPDF=FALSE, plotParams=list(maxrow=50))  # MR modified
   dev.off                                         # MR modified
   
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
   # nInTot is the number of inputs prior (!!!) to optimization after preprocessing
   nInTot <- length(which(model$interMat == -1))
   #nInputs <- length(which(model$interMat == -1))


   #   nInputs is the number of inputs after (!!!) to optimization after preprocessing
   #
   # bs = as.logical(bString)
   # modelCut <- list()
   # modelCut$interMat <- model$interMat[, bs]
   # modelCut$reacID <- model$reacID[bs]
   # modelCut$namesSpecies <- model$namesSpecies
   # nInputs <- length(which(modelCut$interMat == -1))
   #
   #
   #
   nInputs=length(which(bString==1))

   nDataPts <- dim(midas@signals[[timeIndex]])[1] * dim(midas@signals[[timeIndex]])[2]
   sizePen <- (nDataPts * sizeFac * nInputs)/nInTot

   #************* calculate NAPen
   NAFac = 1
   NAPen <- NAFac * length(which(is.na(simResults)))

   #************* sum all scores to produce final score
   # achtung: the total after computeScoreT1 is still divideb by numMeasurements in getFit
   #    # achtung: the total after computeScoreT1 is still divided by numMeasurements after getFit in computeScoreT1   # MR modified

   numMeasurements=sum(!is.na(midas@signals[[timeIndex]]))
   score = (deviationPen + NAPen + sizePen) / numMeasurements
   scores=list("score"=score,
               "deviationPen"=deviationPen,
               "NAPen"=NAPen,
               "sizePen"=sizePen,
               "numMeasurements"=numMeasurements)
   
   
   
   return(scores)  
   
}


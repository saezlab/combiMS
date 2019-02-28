# Final version of the optimization algorithm for a single patient. All parameters applied to final version
# Marti Bernardo-Faura, Juni 2015

# Minor changes to adjust the paths to the combiMS Github project
# Jakob Wirbel, June 2017


# Small change to correct the brackets placing concerning the calculation of rmse (rmseAll)
# Modified to use for CombiMS-rerun on the RWTH Aachen cluster
# path adaption
# 
# Melanie Rinas, November 2017









OptCombiMScluster=function(fileName,indexRep){                 
  

  DonorID_file = fileName    # example "CH032.csv"
  indexRep = indexRep        # example 1
  
  DonorID=unlist(strsplit(DonorID_file,"\\."))[1] 
  
  Center = substr(DonorID, start = 1, stop = 2)
  
  if(Center == "CH"){
    CenterID = 1
  }
  if(Center == "IB"){
    CenterID = 2
  }
  if(Center == "KI"){
    CenterID = 3
  }
  if(Center == "UZ"){
    CenterID = 4
  }
  
  DonorCenterNo = as.numeric(substr(DonorID, start = 3, stop = 5))
  
  
  number__merged_by__CenterID_DonorCenterNo_indexRep = as.numeric(paste0(CenterID,DonorCenterNo,indexRep,sep=""))
  
  
  
  
  
  
  # seed_number_used = number__merged_by__CenterID_DonorCenterNo_indexRep
  # set.seed(seed_number_used)
  
  
  
  
  
  
  
  Results_storage_name = "OptCombiMScluster"
  
  
  
  
  
  
  
  
  library(CellNOptR)
  #library(NMF)
  #setwd("/Users/marti/Documents/R/combiMS/modelling/final")
  
  #source("/Users/marti/Documents/R/combiMS/normaliseSimp.R")
  #fileName="/Users/marti/Documents/ebi/combiMS/data/midas_clean_merged/CH003_phosphos_midas.csv"
  cat(sprintf("*********The patient is %s",fileName))
  
  
 
  
  # *************************** define working directory and folder location
  
  #working_directory = "/rwthfs/rz/cluster/home/mr352506/R/CombiMS_rerun"
  working_directory =dirname(rstudioapi::getActiveDocumentContext()$path)    
  
  ifelse(!dir.exists(file.path(working_directory,"Results")), dir.create(file.path(working_directory,"Results")), FALSE)
  sub_dir__Results = file.path(working_directory,"Results")
  
  ifelse(!dir.exists(file.path(sub_dir__Results,Results_storage_name)), dir.create(file.path(sub_dir__Results,Results_storage_name)), FALSE)
  sub_dir__Results__current_analysis = file.path(sub_dir__Results,Results_storage_name)
  
  

  setwd(sub_dir__Results__current_analysis)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  sub_dir__Data = file.path("../../data/phosphos_processed/")
  sub_dir__Data_MIDAS_used = sub_dir__Data
  
  sub_dir__Data_SIF_files = file.path("../../files/model/")
  sub_dir__Data_SIF_used = file.path(sub_dir__Data_SIF_files,"combiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED.sif")


  
  # *************************** load model
  model=readSIF(sub_dir__Data_SIF_used)             

  
  # numInteractions=length(model$reacID)

  # *************************** load processed midas
  # if already processed (by calling processCombiData.R in OptAllPatients.R), load the processed file
  #taula=read.csv(paste("/Users/marti/Documents/ebi/combiMS/data/midas_clean_merged/processedForCNO/normalized/",fileName,sep=""),header=TRUE,dec=".",check.names=FALSE, stringsAsFactors=FALSE)
  
  
  # *************************** normalise and plot
  #taulaHill2=normaliseSimp(taula,HillCoef=2)
  #write.table(taulaHill2,file="/Users/marti/Documents/R/combiMS/modelling/midas.csv",sep=",",row.names=FALSE,quote=F)
  #midas=CNOlist("/Users/marti/Documents/R/combiMS/modelling/midas.csv")
  midas=CNOlist(paste(sub_dir__Data_MIDAS_used,"/",fileName,sep=""))
  # plot(midas)
  #midas@signals[[2]][,grep('^DIG1',colnames(midas@signals[[1]]))]
  
  
  
  
  
  
  
  
  
  
  

  cat(sprintf("*********After cluster compressing and expanding logic gates, the model has %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID)))
  #plotModel(model,midas)
  
  
  
  
  
  
  
  # *************************** fit model
  initBstring=rep(1,length(model$reacID))
  #cutAndPlot(midas,model, list(initBstring)) #simulate
  cat("********Optimizing",sep="\n")
  #Opt=gaBinaryLargeNw(CNOlist=midas,model=model,initBstring=initBstring,stallGenMax=600000,maxTime=36000, maxGens=100000, verbose=FALSE,popSize=100,elitism=2,relTol=0.05)
  
  
  
  
  # Opt=gaBinaryT1(CNOlist=midas,
  #                model=model,
  #                initBstring=initBstring,
  #                stallGenMax=600000,
  #                maxTime=36000, 
  #                maxGens=100000, 
  #                verbose=FALSE,
  #                popSize=100,
  #                elitism=2,
  #                relTol=0.05)
  
  Opt=gaBinaryT1(                    
    CNOlist=midas,                                
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
      rmseAll[i]= sqrt(    (sum( (midas@signals[[2]][,i] - sim$simResults[[1]]$t1[,i] )^2 ) /numStimuli)   )
  }
  
  
  
  
  fileName2=unlist(strsplit(fileName,"\\."))[1]
  # c=list("patientName"=fileName2,"averageNw"=averageNw,"rmseAll"=rmseAll,"countInteractions"=countInteractions,"countNws"=dim(Opt$stringsTol)[1],"Opt"=Opt)
  # 
  # return(c)
  # 
  
  
  Opt_results_list=list("patientName"=fileName2,                         
                        "averageNw"=averageNw,
                        "rmseAll"=rmseAll,                         
                        "countInteractions"=countInteractions,
                        "countNws"=dim(Opt$stringsTol)[1],
                        "Opt"=Opt)   
  
  
  
  

  
  return(Opt_results_list)  
  
}  











# ****************** to parallelize optimisation in the cluster, this script receives as an argument the patient
argsJob = commandArgs(trailingOnly=TRUE)
fileName = toString(argsJob[1])
repIndex = as.numeric(argsJob[2])
# ****************** use data from each patients to optimize one model per patient
#midas=CNOlist("../processedForCNO/normalized/IB015.csv")
#setwd("/homes/bernardo/combiMS/processedForCNO/normalized/")
#filenames=list.files("/homes/bernardo/combiMS/processedForCNO/normalized/",pattern="*.csv",full.names=TRUE)
#   
# networks=list()
# for (i in 1:length(filenames)){
#   
#   networks[[i]]=OptCombiMScluster(filenames[i])
# }


cat("*********** the patient is",fileName,sep="\n")
cat("*********** the repetition is",repIndex,sep="\n")

network=OptCombiMScluster(fileName,repIndex)             



#**************

fileName2=unlist(strsplit(fileName,"\\."))[1]
cat("*********** saving patient ",fileName2,sep="\n")


#save(network,file=paste("../../files/modeling/cluster/",repIndex,"_",fileName2,".RData",sep=""))




#working_directory = "/rwthfs/rz/cluster/home/mr352506/R/CombiMS_rerun"
working_directory =dirname(rstudioapi::getActiveDocumentContext()$path)    

Results_storage_name = "OptCombiMScluster"


ifelse(!dir.exists(file.path(working_directory,"Results")), dir.create(file.path(working_directory,"Results")), FALSE)
sub_dir__Results = file.path(working_directory,"Results")

ifelse(!dir.exists(file.path(sub_dir__Results,Results_storage_name)), dir.create(file.path(sub_dir__Results,Results_storage_name)), FALSE)
sub_dir__Results__current_analysis = file.path(sub_dir__Results,Results_storage_name)

ifelse(!dir.exists(file.path(sub_dir__Results__current_analysis,paste("CombiMS_rerun_all_donors" ,sep=""))), dir.create(file.path(sub_dir__Results__current_analysis,paste("CombiMS_rerun_all_donors" ,sep=""))), FALSE)
sub_dir__Results__all_donors = file.path(file.path(sub_dir__Results__current_analysis,paste("CombiMS_rerun_all_donors" ,sep="")))


save(network,file=paste(sub_dir__Results__all_donors,"/",repIndex,"_",fileName2,".RData",sep=""))

cat("*********** Script finished! ",sep="\n")



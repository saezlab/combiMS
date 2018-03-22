# Final version of the optimization algorithm for a single patient. All parameters applied to final version
# Marti Bernardo-Faura, Juni 2015

# Minor changes to adjust the paths to the combiMS Github project
# Jakob Wirbel, June 2017



# Modified to use for CombiMS-rerun
# on the RWTH Aachen cluster
# Melanie Rinas, November 2017



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# ISSUE !!!!!:
# model=preprocessing(midas,model,expansion=F,compression=TRUE) with model = combiMSplane.sif
# is not identical with
# model = combiMSplaneCUT.sif
#
# Thus taking model = combiMSplaneCUT.sif (without preprocessing)





# ERROR detection by Melanie Rinas  
# 
# Brackets mistake in the calculation of 
# rmseAll
# compare below the versions rmseAll versus rmseAll_MR













#OptCombiMScluster_MR.R




























# OptCombiMScluster=function(fileName){
OptCombiMScluster_MR=function(fileName,indexRep){                 # MRmodified
  
  # MR inserted start
  
  
  
  
  
  
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
  
  
  
  
  
  
  
  Results_storage_name = "OptCombiMScluster_MR"
  
  
  # MR inserted end
  
  
  
  
  
  
  
  
  library(CellNOptR)
  #library(NMF)
  #setwd("/Users/marti/Documents/R/combiMS/modelling/final")
  
  #source("/Users/marti/Documents/R/combiMS/normaliseSimp.R")
  #fileName="/Users/marti/Documents/ebi/combiMS/data/midas_clean_merged/CH003_phosphos_midas.csv"
  cat(sprintf("*********The patient is %s",fileName))
  
  
  # MR inserted start
  
  # *************************** define working directory and folder location
  
  working_directory = "/rwthfs/rz/cluster/home/mr352506/R/CombiMS_rerun"
  
  
  ifelse(!dir.exists(file.path(working_directory,"Results")), dir.create(file.path(working_directory,"Results")), FALSE)
  sub_dir__Results = file.path(working_directory,"Results")
  
  ifelse(!dir.exists(file.path(sub_dir__Results,Results_storage_name)), dir.create(file.path(sub_dir__Results,Results_storage_name)), FALSE)
  sub_dir__Results__current_analysis = file.path(sub_dir__Results,Results_storage_name)
  
  
  
  # ifelse(!dir.exists(file.path(sub_dir__Results__current_analysis,paste("CombiMS_rerun_", DonorID ,sep=""))), dir.create(file.path(sub_dir__Results__current_analysis,paste("CombiMS_rerun_", DonorID ,sep=""))), FALSE)
  # sub_dir__Results__current_donor = file.path(file.path(sub_dir__Results__current_analysis,paste("CombiMS_rerun_", DonorID ,sep="")))
  # 
  # setwd(sub_dir__Results__current_donor)
  
  setwd(sub_dir__Results__current_analysis)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  sub_dir__Data = file.path(working_directory,"Data")
  sub_dir__Data_MIDAS_files = file.path(sub_dir__Data,"MIDAS_files")
  sub_dir__Data_MIDAS_used = file.path(sub_dir__Data_MIDAS_files,"phosphos_processed__original_MB_JW")
  
  sub_dir__Data_SIF_files = file.path(sub_dir__Data,"SIF_files")
  #sub_dir__Data_SIF_used = file.path(sub_dir__Data_SIF_files,"combiMSplane.sif")
  sub_dir__Data_SIF_combiMSplaneCUT = file.path(sub_dir__Data_SIF_files,"combiMSplaneCUT.sif")
  
  # MR inserted end
  
  
  
  # *************************** load model
  # model=readSIF("../../files/model/combiMSplane.sif")
  # model=readSIF(sub_dir__Data_SIF_used)               # MRmodified
  model_combiMSplaneCUT=readSIF(sub_dir__Data_SIF_combiMSplaneCUT)               # MRmodified
  
  
  # numInteractions=length(model$reacID)
  # numNode= length(model$namesSpecies)  # MR inserted
  # number_of_stimuli=dim(midas@signals[[2]])[1]   # MR inserted
  
  
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
  
  
  
  
  
  
  
  
  
  
  
  # ISSUE !!!!!:
  # model=preprocessing(midas,model,expansion=F,compression=TRUE) with model = combiMSplane.sif
  # is not identical with
  # model = combiMSplaneCUT.sif
  
  
  
  # # *************************** format model
  # # compress and expand model
  # cat(sprintf("*********Before cluster processing, the model has %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID)))
  # model=preprocessing(midas,model,expansion=F,compression=TRUE)
  # cat(sprintf("*********After cluster compressing and expanding logic gates, the model has %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID)))
  # #plotModel(model,midas)
  # 
  # 
  # *********The patient is CH003.csv*********Before cluster processing, the model has 168 nodes and 296 reactions[1] "The following species are measured: HSPB1, JUN, PTN11, IKBA, AKT1, WNK1, MK12, TF65, STAT1, GSK3A, STAT3, MP2K1, MKO3, CREB1, FAK1, STAT5, STAT6"
  # [1] "The following species are stimulated: ANTICD3, BDNF, BN201, GILENYA, H2O2, REBIF, IFNG, NACL, IL1A, TERIFLUNOMIDE, IL6, INS, LPS, EGCG, POLYIC, CONA, S1P1, TNFA, VITD3, DMF"
  # [1] "The following species are inhibited: "
  # [1] "The following species are not observable and/or not controllable: BCL10, CABIN, CALCINEURIN, CARD11, CARD11A, CD19, CD28, CD4, CD45, CD8, CDC42, CTLA4, EGF, GAP, ICAM1, IL2, IL4, INTEGRIN, JAK3, KCC4, MALT1, PGFB, X, FOXO1, FOXO3, P21, PSA6, NFAT, SRE, CRE, BCAT, STAT2, MRP2, NQO1, IRF1, PRDM1"
  # *********After cluster compressing and expanding logic gates, the model has 72 nodes and 170 reactions********Optimizing
  # [1] 137
  # 
  
  
  
  model=model_combiMSplaneCUT
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
  #rmseAll=numeric(length=numSignals) 
  rmseAll_MR=numeric(length=numSignals)   # MR inserted
  
  # MR NOTE: WRONG BRACKET PLACING!!!
  
  # for(i in 1:numSignals){
  #    rmseAll[i]=sqrt((sum(midas@signals[[2]][,i]-sim$simResults[[1]]$t1[,i])^2)/numStimuli)
  # }
  
  # MR inserted start ............................................................................................................................................................................................................................................................................................................
  
  for(i in 1:numSignals){
    rmseAll_MR[i]=sqrt(    (sum( (midas@signals[[2]][,i] - sim$simResults[[1]]$t1[,i] )^2 ) /numStimuli)   )
  }
  # MR inserted end ............................................................................................................................................................................................................................................................................................................
  
  
  
  fileName2=unlist(strsplit(fileName,"\\."))[1]
  # c=list("patientName"=fileName2,"averageNw"=averageNw,"rmseAll"=rmseAll,"countInteractions"=countInteractions,"countNws"=dim(Opt$stringsTol)[1],"Opt"=Opt)
  # 
  # return(c)
  # 
  
  
  Opt_results_list=list("patientName"=fileName2,                         # MRmodified
                        "averageNw"=averageNw,
                        "rmseAll_MR"=rmseAll_MR,                         # MRmodified
                        "countInteractions"=countInteractions,
                        "countNws"=dim(Opt$stringsTol)[1],
                        "Opt"=Opt)   
  
  
  
  
  # # MR inserted start
  # 
  # 
  # save(Opt_results_list, file=paste(sub_dir__Results__current_donor , "/Opt_results_list__donor_",DonorID,"_","indexRep_" ,indexRep ,".RData" ,sep=""))
  # 
  # save.image(file=paste(sub_dir__Results__current_donor , "/workspace__fit__donor_",DonorID,"_","indexRep_" ,indexRep ,".RData" ,sep=""))
  # 
  # 
  # # MR inserted end
  
  return(Opt_results_list)  # MRmodified
  
}  # end of function OptCombiMScluster_MR=function(fileName,indexRep){                 # MRmodified











# Copy and paste from Marti's script clusterCombiMS.R  start:


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

network=OptCombiMScluster_MR(fileName,repIndex)              # MRmodified
#OptCombiMScluster_MR=function(fileName,indexRep){                 # MRmodified



#**************

fileName2=unlist(strsplit(fileName,"\\."))[1]
cat("*********** saving patient ",fileName2,sep="\n")


#save(network,file=paste("../../files/modeling/cluster/",repIndex,"_",fileName2,".RData",sep=""))




# Copy and paste from Marti's script clusterCombiMS.R  end



# MR insert start


working_directory = "/rwthfs/rz/cluster/home/mr352506/R/CombiMS_rerun"
Results_storage_name = "OptCombiMScluster_MR"


ifelse(!dir.exists(file.path(working_directory,"Results")), dir.create(file.path(working_directory,"Results")), FALSE)
sub_dir__Results = file.path(working_directory,"Results")

ifelse(!dir.exists(file.path(sub_dir__Results,Results_storage_name)), dir.create(file.path(sub_dir__Results,Results_storage_name)), FALSE)
sub_dir__Results__current_analysis = file.path(sub_dir__Results,Results_storage_name)

ifelse(!dir.exists(file.path(sub_dir__Results__current_analysis,paste("CombiMS_rerun_all_donors" ,sep=""))), dir.create(file.path(sub_dir__Results__current_analysis,paste("CombiMS_rerun_all_donors" ,sep=""))), FALSE)
sub_dir__Results__all_donors = file.path(file.path(sub_dir__Results__current_analysis,paste("CombiMS_rerun_all_donors" ,sep="")))


save(network,file=paste(sub_dir__Results__all_donors,"/",repIndex,"_",fileName2,".RData",sep=""))

cat("*********** Script finished! ",sep="\n")

# MR insert end




#OptCombiMS_MR.R



# Local version of the script Run in cluster for patient optimisation. Called by OptAllPatients.R
# Marti Bernardo Faura, September 2014. 

# Minor changes to adjust the paths to the combiMS Github project
# Jakob Wirbel, June 2017




# Changes to better understand the script and for testing 
# by
# Melanie Rinas. October 2017


# ERROR detection by Melanie Rinas  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 
# Brackets mistake in the calculation of 
# rmseAll
# compare below the versions rmseAll versus rmseAll_MR











OptCombiMS_MR = function(fileName){
  
  
  library(CellNOptR)
  library(NMF)
  # setwd("/Users/marti/Documents/R/combiMS/modelling/final")
  # source("../data_processing_and_normalization/normaliseSimp.R")
  
  # working_directory =dirname(rstudioapi::getActiveDocumentContext()$path)    # MR
  # setwd(working_directory)                                                   # MR
  
  working_directory ="/Users/melanie/Documents/GitHub/combiMS-master_MR/code/single_model_optimization"    # MR
  setwd(working_directory)                                                   # MR
  
  
  
  
  # MR inserted start ............................................................................................................................................................................................................................................................................................................
  
  
  ifelse(!dir.exists(file.path(working_directory,"Results/OptCombiMS_MR")), dir.create(file.path(working_directory,"Results/OptCombiMS_MR")), FALSE)
  sub_dir__Results__current_script = file.path(working_directory,"Results/OptCombiMS_MR")
  
  ifelse(!dir.exists(file.path(sub_dir__Results__current_script,"Figures")), dir.create(file.path(sub_dir__Results__current_script,"Figures")), FALSE)
  sub_dir__Results__current_script_Figures = file.path(sub_dir__Results__current_script,"Figures")
  
  ifelse(!dir.exists(file.path(sub_dir__Results__current_script_Figures,"Figures_plotModel_single_donors")), dir.create(file.path(sub_dir__Results__current_script_Figures,"Figures_plotModel_single_donors")), FALSE)
  sub_dir__Results__current_script_Figures_plotModel = file.path(sub_dir__Results__current_script_Figures,"Figures_plotModel_single_donors")
  
  ifelse(!dir.exists(file.path(sub_dir__Results__current_script_Figures,"Figures_cutAndPlot_single_donors")), dir.create(file.path(sub_dir__Results__current_script_Figures,"Figures_cutAndPlot_single_donors")), FALSE)
  sub_dir__Results__current_script_Figures_cutAndPlot = file.path(sub_dir__Results__current_script_Figures,"Figures_cutAndPlot_single_donors")
  
  ifelse(!dir.exists(file.path(sub_dir__Results__current_script_Figures,"Figures_plotFit_single_donors")), dir.create(file.path(sub_dir__Results__current_script_Figures,"Figures_plotFit_single_donors")), FALSE)
  sub_dir__Results__current_script_Figures_plotModelfittingProcess = file.path(sub_dir__Results__current_script_Figures,"Figures_plotFit_single_donors")
  
  
  
  ifelse(!dir.exists(file.path(sub_dir__Results__current_script,"Results")), dir.create(file.path(sub_dir__Results__current_script,"Results")), FALSE)
  sub_dir__Results__current_script_Results = file.path(sub_dir__Results__current_script,"Results")
  
  
  
  # MR inserted end ............................................................................................................................................................................................................................................................................................................
  
  
  
  
  
  
  
  
  
  
  
  # *************************************************
  # ************find patient data, load its midas
  # **************************************************
  data_folder="../../data/phosphos_processed/"
  patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
  i=1; 
  fileName=patientData[i] #give a specific patient if not passed via function
  
  patient_ID=unlist(strsplit(fileName,"\\."))[1]                      # MR inserted
  
  
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
  numNode= length(model$namesSpecies)  # MR inserted
  number_of_stimuli=dim(midas@signals[[2]])[1]   # MR inserted
  
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
  
  # MR inserted start ............................................................................................................................................................................................................................................................................................................
  
  
  FigureHeight_cutAndPlot = 8
  FigureWidth_cutAndPlot = 20
  
  
  # FigureHeight_cnolist =7
  # FigureWidth_cnolist =14
  
  FigureWidth_plotModel = 16
  FigureHeight_plotModel = 9
  
  
  FigureHeight = 6
  FigureWidth = 10
  
  FontSizeValue=1 
  
  
  pdf(file=file.path(sub_dir__Results__current_script_Figures_cutAndPlot,paste("cutAndPlot_OptCombiMS_simulation_using_initBstring_allReacIDequal1_",patient_ID, ".pdf", sep="")),width=FigureWidth_cutAndPlot, height=FigureHeight_cutAndPlot)
  
  # output_cutAndPlot=cutAndPlotResultsT1(pknmodel, bString=bitstring,CNOlist=CNOlist_data_single_patient, plotPDF=FALSE,plotParams=list(maxrow=25, margin=0.1,cex=FontSizeValue,cex.axis=FontSizeValue,cmap_scale=0.5,width=FigureWidth_cutAndPlot, height=FigureHeight_cutAndPlot))
  output_cutAndPlot_initial=cutAndPlot(CNOlist=midas,model,bString=list(initBstring),plotParams = list(maxrow = (number_of_stimuli+2),cex=.8))
  
  dev.off()
  
  
  # MR inserted end ............................................................................................................................................................................................................................................................................................................
  
  cutAndPlot(midas,model, list(initBstring)) #simulate
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Fitting process -------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  
  indexRep = 3
  
  Opt=gaBinaryT1(CNOlist=midas,model=model,initBstring=initBstring,stallGenMax=600,maxTime=10*60, maxGens=100000, verbose=FALSE,popSize=100,elitism=2)   # MRmodified maxTime                                                                                                              # MRmodified
  
  
                                                                                # MRmodified
  
  
  save(Opt, file=paste(sub_dir__Results__current_script_Results , "/Opt_" ,indexRep ,".RData" ,sep=""))
  
  
  #Opt=gaBinaryT1(CNOlist=midas,model=model,initBstring=initBstring,stallGenMax=600,maxTime=20, maxGens=100000, verbose=FALSE,popSize=100,elitism=2)
  #
  #Opt=gaBinaryLargeNw(CNOlist=midas,model=model,initBstring=initBstring,stallGenMax=600,maxTime=20, maxGens=100000, verbose=FALSE,popSize=100,elitism=2)
  
  
  #***************************plot result
  #plot fits
  
  
  
  # MR inserted start ............................................................................................................................................................................................................................................................................................................
  
  pdf(file=file.path(sub_dir__Results__current_script_Figures_cutAndPlot,paste("cutAndPlot_OptCombiMS_simulation_using_fitted_OptBstring_",patient_ID, ".pdf", sep="")),width=FigureWidth_cutAndPlot, height=FigureHeight_cutAndPlot)
  
  # output_cutAndPlot=cutAndPlotResultsT1(pknmodel, bString=bitstring,CNOlist=CNOlist_data_single_patient, plotPDF=FALSE,plotParams=list(maxrow=25, margin=0.1,cex=FontSizeValue,cex.axis=FontSizeValue,cmap_scale=0.5,width=FigureWidth_cutAndPlot, height=FigureHeight_cutAndPlot))
  output_cutAndPlot_post_fitting=cutAndPlot(CNOlist=midas,model=model,bString=list(Opt$bString),plotPDF=FALSE,plotParams = list(maxrow = (number_of_stimuli+2),cex=.8))
  
  dev.off()
  
  
  # MR inserted end ............................................................................................................................................................................................................................................................................................................
  
  
  
  prova=cutAndPlot(model=model,bStrings=list(Opt$bString),CNOlist=midas, plotPDF=FALSE, plotParams = list(maxrow = 25,cex=.8))
  
  
  
  # MR inserted start ............................................................................................................................................................................................................................................................................................................
  
  pdf(file=file.path(sub_dir__Results__current_script_Figures_plotModelfittingProcess,paste("plotFit_OptCombiMS_simulation_using_fitted_Opt_",patient_ID, ".pdf", sep="")),width=FigureWidth, height=FigureHeight)
  
  # output_cutAndPlot=cutAndPlotResultsT1(pknmodel, bString=bitstring,CNOlist=CNOlist_data_single_patient, plotPDF=FALSE,plotParams=list(maxrow=25, margin=0.1,cex=FontSizeValue,cex.axis=FontSizeValue,cmap_scale=0.5,width=FigureWidth_cutAndPlot, height=FigureHeight_cutAndPlot))
  output_plotFit=plotFit(optRes=Opt)
  
  dev.off()
  
  
  
  
  pdf(file=file.path(sub_dir__Results__current_script_Figures_plotModel,paste("plotModel_OptCombiMS_simulation_using_fitted_OptBstring_",patient_ID, ".pdf", sep="")),width=FigureWidth_plotModel, height=FigureHeight_plotModel)
  
  # output_cutAndPlot=cutAndPlotResultsT1(pknmodel, bString=bitstring,CNOlist=CNOlist_data_single_patient, plotPDF=FALSE,plotParams=list(maxrow=25, margin=0.1,cex=FontSizeValue,cex.axis=FontSizeValue,cmap_scale=0.5,width=FigureWidth_cutAndPlot, height=FigureHeight_cutAndPlot))
  output_plotModel=plotModel(model,midas,Opt$bString)
  
  dev.off()
  
  
  # MR inserted end ............................................................................................................................................................................................................................................................................................................
  
  
  
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
  numNetworks=dim(Opt$stringsTol)[1]        # MR copy and paste from OptCombiMScluster.R
  sprintf("*********the number of networks within relTols is %s",dim(Opt$stringsTol)[1])
  numEdges=dim(Opt$stringsTol)[2]            # MR copy and paste from OptCombiMScluster.R
  cat(sprintf("*******the number of edges in each network is %s",numEdges),sep="\n")      # MR copy and paste from OptCombiMScluster.R 
  
  matrix_stringsTol = Opt$stringsTol    # MR inserted
  
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
  rmseAll_MR=numeric(length=numSignals)   # MR inserted
  
  
  for(i in 1:numSignals){
    rmseAll[i]=sqrt((sum( midas@signals[[2]][,i] - prova$simResults[[1]]$t1[,i] )^2)/numStimuli)
  }
  
  
  
  # MR inserted start ............................................................................................................................................................................................................................................................................................................
  
  for(i in 1:numSignals){
    rmseAll_MR[i]=sqrt(    (sum( (midas@signals[[2]][,i] - prova$simResults[[1]]$t1[,i] )^2 ) /numStimuli)   )
  }
  # MR inserted end ............................................................................................................................................................................................................................................................................................................
  
  
  
  check_rmseAll_versions_of_Marti_and_MR_identical = identical(rmseAll,rmseAll_MR)
  check_rmseAll_versions_of_Marti_and_MR_identical
  
  fileName2=unlist(strsplit(fileName,"\\."))[1]                      # MR copy and paste from OptCombiMScluster.R
  # c=list("patientName"=fileName2,"averageNw"=averageNw,"rmseAll"=rmseAll,"countInteractions"=countInteractions,"countNws"=dim(Opt$stringsTol)[1],"Opt"=Opt)   # MR copy and paste from OptCombiMScluster.R 
  
  #c=list("patientName"=fileName2,"averageNw"=averageNw,"rmseAll_wrong"=rmseAll,"countInteractions"=countInteractions,"countNws"=dim(Opt$stringsTol)[1],"Opt"=Opt)   # MR copy and paste from OptCombiMScluster.R  and modified
  Opt_results_list=list("patientName"=fileName2,"averageNw"=averageNw,"rmseAll_MR"=rmseAll_MR,"countInteractions"=countInteractions,"countNws"=dim(Opt$stringsTol)[1],"Opt"=Opt)   # MR copy and paste from OptCombiMScluster.R  and modified
  
  save(Opt_results_list, file=paste(sub_dir__Results__current_script_Results , "/Opt_results_list_" ,indexRep ,".RData" ,sep=""))
  
  
  #c=list("averageNw"=averageNw,"rmseAll"=rmseAll,"countInteractions"=countInteractions,"countNws"=dim(Opt$stringsTol)[1],"Opt"=Opt)
  
  #return(c)
  return(Opt_results_list)  # MRmodified
  
}























# 
# 
# > load("/Users/melanie/Documents/GitHub/combiMS-master_MR/code/single_model_optimization/Results/OptCombiMS_MR/Results/Opt_3.RData")
# > 
#    > 
#    > Opt_3 = Opt
# > identical(Opt_1$bString,Opt_3$bString)
# [1] TRUE
# >     identical(Opt_1$bScore,Opt_3$bScore)
# [1] TRUE
# >     identical(Opt_1$results,Opt_3$results)
# [1] FALSE
# >         identical(Opt_1$stringsTol,Opt_3$stringsTol)
# [1] TRUE
# >         identical(Opt_1$stringsTolScores,Opt_3$stringsTolScores)
# [1] TRUE
# > Opt_3_results__COLS_1_to_7 = Opt_3_results[,1:7]
# Error: object 'Opt_3_results' not found
# > Opt_3_results = Opt_3$results
# > Opt_3_results__COLS_1_to_7 = Opt_3_results[,1:7]
# >     identical(Opt_1_results__COLS_1_to_7,Opt_3_results__COLS_1_to_7)
# [1] TRUE



# 
# 
# # different  seed_number_used     ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
# 
# 
# > load("/Users/melanie/Documents/GitHub/combiMS-master_MR/code/single_model_optimization/Results/OptCombiMS_MR/Results/Opt_61118.RData")
# > Opt_61118 = Opt
# > identical(Opt_1$bString,Opt_61118)
# [1] FALSE
# > identical(Opt_1$bString,Opt_61118$bString)
# [1] FALSE
# > identical(Opt_1$bScore,Opt_61118$bScore)
# [1] FALSE
# > Opt_61118$bString
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 1 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 1 1 0 1 1
# [85] 1 0 0 0 1 1 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 1 0 0
# [169] 0 0 0 1 0 1 1 1 1 0 0 1 0 1 1 0 1 0
# > Opt_1$bString
# [1] 0 0 0 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 1 0 1 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 1 0 0 1 1
# [85] 1 0 0 0 1 1 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 1 0 0
# [169] 0 0 0 0 0 1 1 1 1 0 1 0 0 1 0 0 0 0
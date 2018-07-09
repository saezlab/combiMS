

# calculateDefective_MR__with_Results_storage_name_Input




# Function to calculate and visualize co-druggable reactions (previously termed defective). Called
# by pathDrugTargetsFinalv3.R
# Marti Bernardo-Faura
# Final version April 2016

# Minor changes for Github repository by Jakob Wirbel, July 2017

#calculateDefective=function(thisMode='mean',drugable="zero"){     # MR modified
calculateDefective_MR__with_Results_storage_name_Input=function(thisMode='mean',drugable="zero",Results_storage_name=format(Sys.time(), "%Y_%m_%d")){     # MR modified
   
   library(ggplot2)
   library(reshape2)
   
   source("./phenotypeNetwork.R")
   
   # phenotypeNws_folder='../../files/group_models/'
   # drugScores_folder='../../files/drugScores/'
   
  
   drugScores_folder='../../files/drugScores'       # MR modified
   phenotypeNws_folder='../../files/group_models'       # MR modified
   
  
   
   ifelse(!dir.exists(file.path(drugScores_folder,Results_storage_name)), dir.create(file.path(drugScores_folder,Results_storage_name)), FALSE)              # MR inserted
   drugScores_folder_of_current_script = file.path(drugScores_folder,Results_storage_name)                                                                              # MR inserted
   drugScores_folder_of_current_script_ = paste(drugScores_folder_of_current_script,"/",sep="")
   
   
   ifelse(!dir.exists(file.path(phenotypeNws_folder,Results_storage_name)), dir.create(file.path(phenotypeNws_folder,Results_storage_name)), FALSE)              # MR inserted
   phenotypeNws_folder_of_current_script = file.path(phenotypeNws_folder,Results_storage_name)                                                                              # MR inserted
   phenotypeNws_folder_of_current_script_ = paste(phenotypeNws_folder_of_current_script,"/",sep="")
   
   
   
   # *************************************************************************************************************************
   # ***********given a patient-model, are we interested in mean or median for the cohort?
   # *************************************************************************************************************************
   warning(paste0(thisMode,' mode ASSUMED to average all patient-models\n'))
   warning(paste0(drugable," scores selected"))
   
   
   
   # *************************************************************************************************************************
   # ***********load anotation
   #*************************************************************************************************************************
   
   # ************load anotation to map patients to groups
   #   data_folder="/Users/marti/Documents/ebi/combiMS/data/phosphosMergedAbsMax/processed/normalized/secondRoundProcessedMidas/"
   #   filenames=list.files(data_folder,pattern="*.csv",full.names=FALSE)
   #   annot=read.csv("/Users/marti/Documents/R/combiMS/modelling/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
   #   numPat=length(filenames)
   #   filenames2=filenames
   #   for (j in 1:numPat){
   #     filenames2[j]=strsplit(filenames[j],"\\.")[[1]][1]
   #   }
   
   # ************load model and midas for annotation
   #   patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
   #   model_path="/Users/marti/Documents/R/combiMS/combiMSplane.sif"
   #   fileName=patientData[1]
   #   midas=CNOlist(paste(data_folder,fileName,sep=""))
   #   model=readSIF(model_path)  
   #   sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
   #   model=preprocessing(midas,model,expansion=FALSE)
   numInteractions=length(model$reacID)
   #   warning(paste0("After compressing: ",length(model$namesSpecies)," nodes, ", numInteractions," reactions"))
   
   
   #writeSIF(model,file="/Users/marti/Documents/R/combiMS/preprocessedModelCombi.sif")
   
   # *************************************************************************************************************************
   # ***********load predicted networks
   # *************************************************************************************************************************
   
   #   load("/Users/marti/Documents/R/combiMS/cluster/analysis/fivePerHundredThousTol2/allMedianModels.RData")
   
   
   
   # *************************************************************************************************************************
   # *********** 1. Calculate the distance between MS and Healthy for each edge
   # *************************************************************************************************************************
   IdxHealthy = annot$Group == 'Healthy'
   IdxMSuntreated = (annot$Treatment == "no" & annot$Disease.Subtype!='')
   
   H = phenotypeNetwork(IdxHealthy, allMedianNetworks,model,mode=thisMode) 
   #phenotypeNetwork<-function(IndexPatients, allMedianNetworks, model, mode='mean', linkActivityThreshold=0.5 ){                
   
   thisPhenotype='Healthy'
   write.table(H$network,file=paste0(phenotypeNws_folder_of_current_script_,thisPhenotype,thisMode,".csv"),sep=",",row.names=T)              # MR modified
   #H_network_matrix = as.matrix(H$network)
   #View(H_network_matrix)
   
   MSuntreatedNw=phenotypeNetwork(IdxMSuntreated, allMedianNetworks,model,mode=thisMode)
   thisPhenotype='MS'
   write.table(MSuntreatedNw$network,file=paste0(phenotypeNws_folder_of_current_script_,thisPhenotype,thisMode,".csv"),sep=",",row.names=T)             # MR modified
   #MSuntreated_network_matrix = as.matrix(MSuntreatedNw$network)
   #View(H_network_matrix)
   
   #distance
   HD=H$network-MSuntreatedNw$network
   #HD_matrix = as.matrix(HD)
   
   
   
   
   
   
   
   # *************************************************************************************************************************
   # *********** 2. Identify edges where difference MS to health is smaller or equal than drug to Health
   # *********** These are drug defective interactions, hence calling for combination therapy that will revert signaling
   # *************************************************************************************************************************
   #Gilenya=Fingolimod
   #Glatiramer=Copaxone
   #Tysabri=Natalizumab
   phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
   numDrugs=5
   
   #prepare dataframe for plot
   phenotypeScores=data.frame(Interactions=colnames(allMedianNetworks),
                              Gilenya=rep(NA,numInteractions),
                              IFNb=rep(NA,numInteractions),
                              Copaxone=rep(NA,numInteractions),
                              EGCG=rep(NA,numInteractions),
                              Tysabri=rep(NA,numInteractions))
   
   phenotypeScores=melt(phenotypeScores,id='Interactions',value.name='Score')
   colnames(phenotypeScores)=c('Interactions','Drug','Score')
   
   numDrugDefective=vector(length=numDrugs)
   names(numDrugDefective)=phenotypes
   
   allDrugNws=list()
   allDefective=list()
   
   for (i in 1:numDrugs){
      thisPhenotype=phenotypes[i]
      patientsThisPhenotype= annot$Treatment==thisPhenotype   #which are the patients in this phenotype?
      
      
      #**************Create phenotype network for Drug i
      phenoNw=phenotypeNetwork(patientsThisPhenotype,allMedianNetworks,model,mode=thisMode) 
      #**************Calculate the distance between i drug-Network and Healthy network
      H2Pheno=H$network-phenoNw$network 
      H2Pheno_matrix = as.matrix(H2Pheno)
      
      #**************Substract the above from Healthy to Ms
      drugScores=abs(HD)-abs(H2Pheno) 
      drugScores_matrix = as.matrix(drugScores)
      
      #phenotypeScores[,i+1]=drugScores
      #**************Fill in dataframe for plotting with scores
      phenotypeScores$Score[which(as.character(phenotypeScores$Drug)==thisPhenotype)]=drugScores
      #barplot(drugScores,ylab='drugScores',xlab='Interactions (original network:186)',main=paste0(thisPhenotype,". Defective:",length(which(drugScores<0))))
      #**************Identify and save defective Ints. Important: there are non-drugable 0 scores,i.e. where situation is ok.
      if(drugable=="negative"){
         defectiveInts=drugScores[which(drugScores<0)]
         allDefective[[i]]=defectiveInts
      } else if(drugable=="zero"){
         #In reactions with 0 score due to situation ok, replace by 1. use names as ids to fulfill two conditions in vectors scores and H2Pheno without messing up
         intsSituationNotOk=names(which(H2Pheno!=0))
         intsSituationOk=names(which(H2Pheno==0))
         ints0AndNotOk=intersect(names(drugScores[which(drugScores==0)]), intsSituationNotOk)
         ints0AndOk=intersect(names(drugScores[which(drugScores==0)]), intsSituationOk)
         intsNeg=which(drugScores<0)
         intsPos=which(drugScores>0)
         warning(paste0("Neg:",length(intsNeg)," Pos:",length(intsPos)," 0AndOk: ",length(ints0AndOk)," 0andNotOk:",length(ints0AndNotOk),"\n"))
         # replace
         drugScores[which(names(drugScores) %in% ints0AndOk)]=1
         defectiveInts=drugScores[which(drugScores<=0)]
         allDefective[[i]]=defectiveInts
         
      }
      
      #**************Calculate number of defective
      numDrugDefective[i]=length(defectiveInts)
      warning(paste0(numDrugDefective[i]," drugable reactions for ",thisPhenotype,"\n"))
      
      #**************Save this network to structure in order to concatenate all networks
      allDrugNws[[i]]=phenotypeNetwork(patientsThisPhenotype,allMedianNetworks,model,mode=thisMode)
      allDrugNws[[i]]$Drug=thisPhenotype
      
      #**************save individual files: for this drug, (i) defective reactions and (ii) network
      write.table(defectiveInts,file=paste0(drugScores_folder_of_current_script_,thisPhenotype,thisMode,drugable,".csv"),sep=",",row.names=T)
      write.table(allDrugNws[[i]]$network,file=paste0(phenotypeNws_folder_of_current_script_,thisPhenotype,thisMode,drugable,".csv"),sep=",",row.names=T)
      
   }
   
   
   #****** which are the interactions that are defective in all drugs?
   alwaysDefective=intersect(intersect(intersect(intersect(names(allDefective[[1]]),names(allDefective[[2]])),names(allDefective[[3]])),names(allDefective[[4]])),names(allDefective[[5]]))
   write.table(alwaysDefective,file=paste0(drugScores_folder_of_current_script_,thisMode,"alwaysDefective.csv"),sep=",",row.names=T)
   
   # *************************************************************************************************************************
   # *********** 3. Visualize defective interactions and phenotype networks
   # *************************************************************************************************************************
   
   #*******plot num of defective for each drug
   numDefectiveDF=data.frame(Drug=character(length=numDrugs),NumDefective=vector(length=numDrugs))
   numDefectiveDF$Drug=phenotypes
   for(i in 1:numDrugs){numDefectiveDF$NumDefective[i]=length(allDefective[[i]])}
   #To prevent ggplot2 from reordering alphabetically the labels according to the factors, which are alphabetical
   #We specify that the factors need to be ordered as they already are
   numDefectiveDF$Drug=factor(numDefectiveDF$Drug,levels=numDefectiveDF$Drug)
   a=ggplot(data=numDefectiveDF,aes(Drug,NumDefective,fill=Drug))+geom_bar(stat='identity',alpha=.5)+theme_bw()
   
   #   #*******plot histogram of num defective by drug
   #   c=ggplot(data=phenotypeScores,aes(Score,fill=Drug))+geom_histogram(position='dodge',alpha=.5)
   #   c=c +theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
   #   
   #   #*******plot density of num defective by drug
   #   c=ggplot(data=phenotypeScores,aes(Score,fill=Drug))+geom_density(alpha=.5)
   #   c=c +theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
   #   
   
   #******plot bars with superposed intearctions to identify common defective
   b=ggplot(data=phenotypeScores,aes(Interactions,Score,fill=Drug))+geom_bar(stat='identity',alpha=0.5)
   #remove background, rotate labels
   b=b + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
   #NOT necessary as this infor is in new plot a. add title including total of defective ints by drug
   # b=b+ggtitle(paste0('Drug-defective Interactions. ',phenotypes[1],':',numDrugDefective[1],"  ",phenotypes[2],':',numDrugDefective[2],"  ",phenotypes[3],':',numDrugDefective[3],"  ",
   #                      phenotypes[4],':',numDrugDefective[4],"  ",phenotypes[5],':',numDrugDefective[5],"  "))
   
   
   defectiveScores=list("defectiveInts"=allDefective,'scorePlot'=b,'numDefectivePlot'=a,'healthyNW'=H,'MSuntreatedNw'=MSuntreatedNw,'allDrugNws'=allDrugNws,'alwaysDefective'=alwaysDefective)
   
   return(defectiveScores)
   #  ModelGraph=list("model"=modelPhenotype,"graph"=graphModel,"network"=phenotypeNw)
   
}

# *************************************************************************************************************************
# *********** 3b. Plot mean networks
# *************************************************************************************************************************
# phenotypeNwsFigure_folder='/Users/marti/Documents/R/combiMS/phenotypeNwsFigure/'
# source("/Users/marti/Documents/R/combiMS/phenotypeNetworkMean.R")
# # thisPhenotype='Healthy'
# # thisPhenotype='MS'
# thisPhenotype=phenotypes[i]
# # IndexPatients=IdxHealthy
# # IndexPatients=IdxMSuntreated
# IndexPatients=patientsThisPhenotype
# networkFigure=phenotypeNetworkMean(IndexPatients,allMedianNetworks,model,mode='mean')
# write.table(networkFigure,file=paste0(phenotypeNwsFigure_folder,thisPhenotype,'Mean',".csv"),sep=",",row.names=T)

# *************************************************************************************************************************
# *********** 4. In case we are interested in the ints where drugs DID work or did not do anything
# *************************************************************************************************************************
# InefectiveInts=drugScores[which(drugScores==0)]
# # reorder is useful for networks averaged with mean rather than median, else order does not apply since they are all 0
# InefectiveInts=InefectiveInts[order(InefectiveInts,decreasing=T)]
# effectSize=0.15
# BestRevertedInts=drugScores[which(drugScores>effectSize)]


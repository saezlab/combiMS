# Function to calculate and visualize co-druggable reactions (previously termed defective). Called
# by pathDrugTargetsFinalv3.R
# Marti Bernardo-Faura
# Final version April 2016

# Minor changes for Github repository by Jakob Wirbel, July 2017
# 
# Changes by Melanie Rinas, August 2018
# 
# 


 
# Copyright information ================================================================================================================================== -->

#  Copyright (c) 2018 - European Molecular Biology Laboratory, European Bioinformatics Institute, UK,
#                       Joint Research Center for Computational Biomedicine (JRC-COMBINE), RWTH-Aachen University, Faculty of Medicine, Aachen, Germany 
# 
#  File author(s): Marti Bernardo-Faura (marti.bernardo.faura@gmail.com), Jakob Wirbel, Melanie Rinas (melrinas@gmail.com) 
# 
#  Distributed under the GPLv3 License. 
#  See accompanying file LICENSE.txt or copy at
#  http://www.gnu.org/licenses/gpl-3.0.html 

#  ======================================================================================================================================================== -->



calculateDefective=function(thisMode='mean',
                            druggable="zero",
                            linkActivityQuantileThreshold_used=0.75,
                            applylinkActivityQuantileThreshold_used ='no'){
  library(ggplot2)
  library(reshape2)
  library(exact2x2)
   
  source("./phenotypeNetwork.R")
  phenotypeNws_folder='../../files/group_models'
  drugScores_folder='../../files/drugScores'
  
  # *************************************************************************************************************************
  # ***********given a patient-model, are we interested in mean or median for the cohort?
  # *************************************************************************************************************************
  warning(paste0(thisMode,' mode ASSUMED to average all patient-models\n'))
  warning(paste0(druggable," scores selected"))
  
  
  applylinkActivityQuantileThreshold__storing_text = ''
  
  
  if(thisMode =='mean' && applylinkActivityQuantileThreshold_used =='no'){     
     
     linkActivityQuantileThreshold_used_text = gsub('\\.', '_', linkActivityQuantileThreshold_used)
     applylinkActivityQuantileThreshold__storing_text = paste('linkActivityQuantileThreshold_',linkActivityQuantileThreshold_used_text,'_NOTroundedRealNumber_',sep="")
     
  } else if (thisMode =='mean' && applylinkActivityQuantileThreshold_used =='yes'){
     
     linkActivityQuantileThreshold_used_text = gsub('\\.', '_', linkActivityQuantileThreshold_used)
     applylinkActivityQuantileThreshold__storing_text = paste('linkActivityQuantileThreshold_',linkActivityQuantileThreshold_used_text,'_rounded_',sep="")
     
  }
  
  
  
  
  drugScores_folder_storage_name = paste("drugScores_based_on__",applylinkActivityQuantileThreshold__storing_text,thisMode,"_",druggable,"_",searchInactiveInts,sep="")        
  
  ifelse(!dir.exists(file.path(drugScores_folder,drugScores_folder_storage_name)), dir.create(file.path(drugScores_folder,drugScores_folder_storage_name)), FALSE)              
  drugScores_folder_of_current_script_settings = file.path(drugScores_folder,drugScores_folder_storage_name)                                                                             
  drugScores_folder_of_current_script_settings_ = paste(drugScores_folder_of_current_script_settings,"/",sep="")
  
  
  phenotypeNws_folder_storage_name = paste("phenotypeNws_based_on__",applylinkActivityQuantileThreshold__storing_text,thisMode,"_",druggable,"_",searchInactiveInts,sep="")        
  
  ifelse(!dir.exists(file.path(phenotypeNws_folder,phenotypeNws_folder_storage_name)), dir.create(file.path(phenotypeNws_folder,phenotypeNws_folder_storage_name)), FALSE)              
  phenotypeNws_folder_of_current_script_settings = file.path(phenotypeNws_folder,phenotypeNws_folder_storage_name)                                                                             
  phenotypeNws_folder_of_current_script_settings_ = paste(phenotypeNws_folder_of_current_script_settings,"/",sep="")
  
  
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
  # *************************************************************************************************************************
  # *********** Identify and visualize drug defective interactions for 
  # 
  #             CASE  I: thisMode=='median'
  #             CASE II: thisMode=='mean'
  #             
  # *************************************************************************************************************************
  # *************************************************************************************************************************
  
  
  # *************************************************************************************************************************
  # *********** CASE  I: thisMode=='median'
  # *************************************************************************************************************************
  
  
  if(thisMode=='median'){
  

     
     # *************************************************************************************************************************
     # *********** I.1. Calculate the distance between MS and Healthy for each edge
     # *************************************************************************************************************************
     IdxHealthy = annot$Group == 'Healthy'
     IdxMSuntreated = (annot$Treatment == "no" & annot$Disease.Subtype!='')
     
     H = phenotypeNetwork(IdxHealthy, allMedianNetworks,model,mode=thisMode,linkActivityQuantileThreshold=linkActivityQuantileThreshold_used,applylinkActivityQuantileThreshold=applylinkActivityQuantileThreshold_used)  
     thisPhenotype='Healthy'
     write.table(H$network,file=paste0(phenotypeNws_folder_of_current_script_settings_,thisPhenotype,"_network.csv"),sep=",",row.names=T)
   
     MSuntreatedNw=phenotypeNetwork(IdxMSuntreated, allMedianNetworks,model,mode=thisMode,linkActivityQuantileThreshold=linkActivityQuantileThreshold_used,applylinkActivityQuantileThreshold=applylinkActivityQuantileThreshold_used) 
     thisPhenotype='MS'
     write.table(MSuntreatedNw$network,file=paste0(phenotypeNws_folder_of_current_script_settings_,thisPhenotype,"_network.csv"),sep=",",row.names=T)
   
     
     #distance
     HD=H$network-MSuntreatedNw$network
     
     absHD = abs(HD)
     thisPhenotype = "AbsDiff_Healthy_MSuntreated"
     write.table(absHD,file=paste0(phenotypeNws_folder_of_current_script_settings_,thisPhenotype,"_network.csv"),sep=",",row.names=T)
     
     
     # *************************************************************************************************************************
     # *********** I.2. Identify drug defective interactions, which have to fullfill the following condition:
     #            
     #            Condition C1: The difference MSuntreated to healthy is smaller or equal than drug to healthy.
     #            
     #
     # *************************************************************************************************************************
     
     
     
     #Gilenya=Fingolimod, Glatiramer=Copaxone
     phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
     numDrugs=5
     phenotypes_colors <- c('IFNb'='#612158', 'Tysabri'='#7A6FAC', 'Gilenya'='#57AB27', 'Copaxone'='#BDCD00', 'EGCG'='#0098A1')
     
     #prepare dataframe for plot
     phenotypeScores=data.frame(Interactions=colnames(allMedianNetworks),Gilenya=rep(NA,numInteractions),
                                IFNb=rep(NA,numInteractions),Copaxone=rep(NA,numInteractions),EGCG=rep(NA,numInteractions),Tysabri=rep(NA,numInteractions))
     phenotypeScores=melt(phenotypeScores,id='Interactions',value.name='Score')
     phenotypeScores$Score_ints0AndOk_fixed_to_1 = phenotypeScores$Score
     colnames(phenotypeScores)=c('Interactions','Drug','Score','Score_ints0AndOk_fixed_to_1')
     
     numDrugDefective=vector(length=numDrugs)
     names(numDrugDefective)=phenotypes
     allDrugNws=list()
     allDefective=list()
     

     
     for (i in 1:numDrugs){
       thisPhenotype=phenotypes[i]
       patientsThisPhenotype= annot$Treatment==thisPhenotype   #which are the patients in this phenotype?
       #**************Create phenotype network for Drug i
       phenoNw=phenotypeNetwork(patientsThisPhenotype,allMedianNetworks,model,mode=thisMode,linkActivityQuantileThreshold=linkActivityQuantileThreshold_used,applylinkActivityQuantileThreshold=applylinkActivityQuantileThreshold_used) 
       #**************Calculate the distance between i drug-Network and Healthy network
       H2Pheno=H$network-phenoNw$network 
       #**************Substract the above from Healthy to Ms
       drugScores=abs(HD)-abs(H2Pheno) 
       #phenotypeScores[,i+1]=drugScores
       drugScores_allInts = drugScores
       write.table(drugScores,file=paste0(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_allInts.csv"),sep=",",row.names=T)
       save(drugScores_allInts,file=paste(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_allInts.RData",sep=""))
       
       
       #**************Fill in dataframe for plotting with scores
       phenotypeScores$Score[which(as.character(phenotypeScores$Drug)==thisPhenotype)]=drugScores
       #barplot(drugScores,ylab='drugScores',xlab='Interactions (original network:186)',main=paste0(thisPhenotype,". Defective:",length(which(drugScores<0))))
      
       
       
       
         
       
       
       
        #**************Identify and save defective Ints. Important: there are non-druggable 0 scores,i.e. where situation is ok.
       if(druggable=="negative"){
         defectiveInts = drugScores[(which(drugScores<0))]  # C1 needs to be fulfilled
         allDefective[[i]]=defectiveInts
         
         
         
         
         #**************Plot DrugScore versus Diff_absHminusDrug
         
         drug_color <- phenotypes_colors[thisPhenotype]
         
         
         df_DrugScore_Diff_absHminusDrug <- data.frame(score=drugScores,
                                                    diff_absHminusDrug=abs(H2Pheno)
         )
         
         df_DrugScore_Diff_absHminusDrug$defectiveInts = (df_DrugScore_Diff_absHminusDrug$score < 0)
         
         c_plot_x_DrugScore_y_Diff_absHminusDrug = ggplot(df_DrugScore_Diff_absHminusDrug, aes(x=score, y=diff_absHminusDrug)) + 
            geom_point(size=2,col=ifelse(df_DrugScore_Diff_absHminusDrug$defectiveInts, drug_color, 'gray50')) + 
            xlab('Score') +
            ylab('abs( H - Drug)') +
            # xlim(xlim_DrugScore__min,xlim_DrugScore__max)+
            theme_bw() + 
            geom_hline(yintercept = 0 , linetype="dashed")+
            geom_vline(xintercept = 0 , linetype="dashed")
         
         ggsave(c_plot_x_DrugScore_y_Diff_absHminusDrug, file=paste(drugScores_folder_of_current_script_settings_,thisPhenotype,'__x_Score_y_Diff_absHminusDrug.pdf',sep=""), width = 5, height = 3)
         
         
         
         phenotypeScores_Neg_NonNeg = phenotypeScores
         phenotypeScores_Neg_NonNeg$Score_Neg_NonNeg = phenotypeScores$Score
         phenotypeScores_Neg_NonNeg$Score_Neg_NonNeg[which(phenotypeScores_Neg_NonNeg$Score_Neg_NonNeg == 0)] = 1


         # b_Neg_NonNeg=ggplot()+
         #    geom_bar(data=phenotypeScores_Neg_NonNeg,aes(Interactions,Score_Neg_NonNeg,fill=Drug),stat='identity',alpha=1)+
         #    scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
         #    #geom_rect(data = rects, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.1) +
         #    #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0, alpha=0.2, fill="blue")+
         #    #annotate("rect", xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf, alpha=0.2, fill="grey") +
         #    #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0), alpha=0.05, fill="blue")+
         #    #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf), alpha=0.05, fill="grey")+
         #    # geom_rect(fill = 'blue', xmin = -Inf, xmax = Inf, ymin =-Inf, ymax = 0, alpha =0.05)+
         #    # geom_rect(fill = 'white', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)+
         #    geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0),
         #              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="orangered", alpha=0.25)+
         #    geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf),
         #              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray50", alpha=0.25)+
         #    geom_bar(data=phenotypeScores_Neg_NonNeg,aes(Interactions,Score_Neg_NonNeg,fill=Drug),stat='identity',alpha=1)+
         #    #geom_bar(stat='identity',alpha=1)+
         #    geom_hline(yintercept = 0)+
         #    scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)=
         #    labs(x = "Model interactions", y = 'Druggable (negative) versus non-druggable (non-negative) \nscores')
         # #geom_rect(fill = 'red', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)
         # #remove background, rotate labels
         # b_Neg_NonNeg=b_Neg_NonNeg + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
         # #print(b_Neg_NonNeg)
         # 
         # ggsave(b_Neg_NonNeg, file=paste(drugScores_folder_of_current_script_settings_,'plot_x_Ints_y_AllDrugs_CoDrugNonCoDrugScores.pdf',sep=""), width = 5, height = 3)
         # 
         # 
         # b_Scores_CoDrug_NonCoDrug = b_Neg_NonNeg
         # 
         
         
         
       } else if(druggable=="zero"){
         #In reactions with 0 score due to situation ok, replace by 1. use names as ids to fulfill two conditions in vectors scores and H2Pheno without messing up
         intsSituationNotOk=names(which(H2Pheno!=0))
         intsSituationOk=names(which(H2Pheno==0))
         ints0AndNotOk=intersect(names(drugScores[which(drugScores==0)]), intsSituationNotOk)
         ints0AndOk=intersect(names(drugScores[which(drugScores==0)]), intsSituationOk)
         intsNeg=which(drugScores<0)
         intsPos=which(drugScores>0)
         warning(paste0("Neg: ",length(intsNeg),"\nPos: ",length(intsPos),"\n0AndOk: ",length(ints0AndOk),"\n0andNotOk: ",length(ints0AndNotOk),"\n"))
         # replace
         drugScores[which(names(drugScores) %in% ints0AndOk)]=1
         defectiveInts = drugScores[(which(drugScores<=0))]  # C1 needs to be fulfilled
         allDefective[[i]]=defectiveInts
       
         
         phenotypeScores$Score_ints0AndOk_fixed_to_1[which(as.character(phenotypeScores$Drug)==thisPhenotype)]=drugScores
         
         
         #**************Plot DrugScore versus Diff_absHminusDrug
         
         drug_color <- phenotypes_colors[thisPhenotype]
         
         
         df_DrugScore_Diff_absHminusDrug <- data.frame(score=drugScores,
                                                    diff_absHminusDrug=abs(H2Pheno)
         )
         
         df_DrugScore_Diff_absHminusDrug$defectiveInts = (df_DrugScore_Diff_absHminusDrug$score <= 0)
         
         c_plot_x_DrugScore_y_Diff_absHminusDrug = ggplot(df_DrugScore_Diff_absHminusDrug, aes(x=score, y=diff_absHminusDrug)) + 
            geom_point(size=2,col=ifelse(df_DrugScore_Diff_absHminusDrug$defectiveInts, drug_color, 'gray50')) + 
            xlab('Score') +
            ylab('abs( H - Drug)') +
            # xlim(xlim_DrugScore__min,xlim_DrugScore__max)+
            theme_bw() + 
            geom_hline(yintercept = 0 , linetype="dashed")+
            geom_vline(xintercept = 0 , linetype="dashed")
         
         ggsave(c_plot_x_DrugScore_y_Diff_absHminusDrug, file=paste(drugScores_folder_of_current_script_settings_,thisPhenotype,'__x_Score_y_Diff_absHminusDrug.pdf',sep=""), width = 5, height = 3)
         
         
         
         
         phenotypeScores_NonPos_Pos = phenotypeScores
         phenotypeScores_NonPos_Pos$Score_NonPos_Pos = phenotypeScores$Score_ints0AndOk_fixed_to_1
         phenotypeScores_NonPos_Pos$Score_NonPos_Pos[which(phenotypeScores_NonPos_Pos$Score_NonPos_Pos == 0)] = -1


         # b_NonPos_Pos=ggplot()+
         #    geom_bar(data=phenotypeScores_NonPos_Pos,aes(Interactions,Score_NonPos_Pos,fill=Drug),stat='identity',alpha=1)+
         #    scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
         #    #geom_rect(data = rects, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.1) +
         #    #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0, alpha=0.2, fill="blue")+
         #    #annotate("rect", xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf, alpha=0.2, fill="grey") +
         #    #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0), alpha=0.05, fill="blue")+
         #    #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf), alpha=0.05, fill="grey")+
         #    # geom_rect(fill = 'blue', xmin = -Inf, xmax = Inf, ymin =-Inf, ymax = 0, alpha =0.05)+
         #    # geom_rect(fill = 'white', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)+
         #    geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0),
         #              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="orangered", alpha=0.5)+
         #    geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf),
         #              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray50", alpha=0.25)+
         #    geom_bar(data=phenotypeScores_NonPos_Pos,aes(Interactions,Score_NonPos_Pos,fill=Drug),stat='identity',alpha=1)+
         #    #geom_bar(stat='identity',alpha=1)+
         #    geom_hline(yintercept = 0)+
         #    scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)
         # #geom_rect(fill = 'red', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)
         # #remove background, rotate labels
         # b_NonPos_Pos=b_NonPos_Pos + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
         # #print(b_NonpositiveScores)
         # 
         # ggsave(b_NonPos_Pos, file=paste(drugScores_folder_of_current_script_settings_,'plot_x_Ints_y_AllDrugs_CoDrugNonCoDrugScores.pdf',sep=""), width = 5, height = 3)
         # 
         # 
         # 
         # 
         # 
         # b_Scores_CoDrug_NonCoDrug = b_NonPos_Pos
         
         
         
       }
       
       #**************Calculate number of defective
       numDrugDefective[i]=length(defectiveInts)
       warning(paste0(numDrugDefective[i]," druggable reactions for ",thisPhenotype,"\n"))
       
       
       

       
       
       #**************Save this network to structure in order to concatenate all networks
       allDrugNws[[i]]=phenotypeNetwork(patientsThisPhenotype,allMedianNetworks,model,mode=thisMode,linkActivityQuantileThreshold=linkActivityQuantileThreshold_used,applylinkActivityQuantileThreshold=applylinkActivityQuantileThreshold_used) 
       allDrugNws[[i]]$Drug=thisPhenotype
       

       
       
       #**************save individual files: for this drug, (i) drugScores, (ii) defective reactions and (iii) network
       
       
       
       write.table(drugScores,file=paste0(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_allInts__modified_ints0AndOk_fixed_to_1.csv"),sep=",",row.names=T)
       
       write.table(defectiveInts,file=paste0(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_defectiveInts.csv"),sep=",",row.names=T)
   
       write.table(allDrugNws[[i]]$network,file=paste0(phenotypeNws_folder_of_current_script_settings_,thisPhenotype,"_network.csv"),sep=",",row.names=T)
   
       
     }
     
   
     #****** which are the interactions that are defective in all drugs?
     alwaysDefective=intersect(intersect(intersect(intersect(names(allDefective[[1]]),names(allDefective[[2]])),names(allDefective[[3]])),names(allDefective[[4]])),names(allDefective[[5]]))
     write.table(alwaysDefective,file=paste0(drugScores_folder_of_current_script_settings_,"alwaysDefective.csv"),sep=",",row.names=T)
   
     
     
     
     
     
     
     
   # *************************************************************************************************************************
     # *********** I.3. Visualize defective interactions and phenotype networks
     # *************************************************************************************************************************
     
     #*******plot num of defective for each drug
     numDefectiveDF=data.frame(Drug=character(length=numDrugs),NumDefective=vector(length=numDrugs))
     numDefectiveDF$Drug=phenotypes
     for(i in 1:numDrugs){numDefectiveDF$NumDefective[i]=length(allDefective[[i]])}
     #To prevent ggplot2 from reordering alphabetically the labels according to the factors, which are alphabetical
     #We specify that the factors need to be ordered as they already are
     numDefectiveDF$Drug=factor(numDefectiveDF$Drug,levels=numDefectiveDF$Drug)

     a=ggplot(data=numDefectiveDF,aes(Drug,NumDefective,fill=Drug))+
        geom_bar(stat='identity',alpha=1)+
        scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE) +
        theme_bw()+
        theme(axis.text=element_text(size=10),axis.title=element_text(size=10))
     
     ggsave(a, file=paste(drugScores_folder_of_current_script_settings_,'plot_NumDefective.pdf',sep=""), width = 4, height = 3)
     
   #   #*******plot histogram of num defective by drug
   #   c=ggplot(data=phenotypeScores,aes(Score,fill=Drug))+geom_histogram(position='dodge',alpha=.5)
   #   c=c +theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
   #   
   #   #*******plot density of num defective by drug
   #   c=ggplot(data=phenotypeScores,aes(Score,fill=Drug))+geom_density(alpha=.5)
   #   c=c +theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
   #   
     
     #******plot bars with superposed intearctions to identify common defective
     b=ggplot(data=phenotypeScores,aes(Interactions,Score,fill=Drug))+
        geom_bar(stat='identity',alpha=1)+
        geom_hline(yintercept = 0)+
        scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)
     #remove background, rotate labels
     b=b + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
     #NOT necessary as this infor is in new plot a. add title including total of defective ints by drug
     # b=b+ggtitle(paste0('Drug-defective Interactions. ',phenotypes[1],':',numDrugDefective[1],"  ",phenotypes[2],':',numDrugDefective[2],"  ",phenotypes[3],':',numDrugDefective[3],"  ",
     #                      phenotypes[4],':',numDrugDefective[4],"  ",phenotypes[5],':',numDrugDefective[5],"  "))
     
     
     ggsave(b, file=paste(drugScores_folder_of_current_script_settings_,'plot_x_Ints_y_AllDrugsScores.pdf',sep=""), width = 20, height = 4.5)
     
     
     
 
     
     if(druggable=="negative"){
        
        
        b_Neg_NonNeg=ggplot()+
           geom_bar(data=phenotypeScores_Neg_NonNeg,aes(Interactions,Score_Neg_NonNeg,fill=Drug),stat='identity',alpha=1)+
           scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
           #geom_rect(data = rects, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.1) +
           #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0, alpha=0.2, fill="blue")+
           #annotate("rect", xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf, alpha=0.2, fill="grey") +
           #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0), alpha=0.05, fill="blue")+
           #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf), alpha=0.05, fill="grey")+
           # geom_rect(fill = 'blue', xmin = -Inf, xmax = Inf, ymin =-Inf, ymax = 0, alpha =0.05)+
           # geom_rect(fill = 'white', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)+
           geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0),
                     aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="orangered", alpha=0.25)+
           geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf),
                     aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray50", alpha=0.25)+
           geom_bar(data=phenotypeScores_Neg_NonNeg,aes(Interactions,Score_Neg_NonNeg,fill=Drug),stat='identity',alpha=1)+
           #geom_bar(stat='identity',alpha=1)+
           geom_hline(yintercept = 0)+
           #scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
              labs(x = "Model interactions", y = 'Druggable (negative) \nversus \nnon-druggable (non-negative) \nscores')
           #geom_rect(fill = 'red', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)
           #remove background, rotate labels
           b_Neg_NonNeg=b_Neg_NonNeg + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
           #print(b_Neg_NonNeg)
           
           ggsave(b_Neg_NonNeg, file=paste(drugScores_folder_of_current_script_settings_,'plot_x_Ints_y_AllDrugs_CoDrugNonCoDrugScores.pdf',sep=""), width = 20, height = 4.5)
           
           
           b_Scores_CoDrug_NonCoDrug = b_Neg_NonNeg
           
     }
     
     if(druggable=="zero"){
        
        b_NonPos_Pos=ggplot()+
           geom_bar(data=phenotypeScores_NonPos_Pos,aes(Interactions,Score_NonPos_Pos,fill=Drug),stat='identity',alpha=1)+
           scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
           #geom_rect(data = rects, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.1) +
           #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0, alpha=0.2, fill="blue")+
           #annotate("rect", xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf, alpha=0.2, fill="grey") +
           #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0), alpha=0.05, fill="blue")+
           #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf), alpha=0.05, fill="grey")+
           # geom_rect(fill = 'blue', xmin = -Inf, xmax = Inf, ymin =-Inf, ymax = 0, alpha =0.05)+
           # geom_rect(fill = 'white', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)+
           geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0),
                     aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="orangered", alpha=0.5)+
           geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf),
                     aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray50", alpha=0.25)+
           geom_bar(data=phenotypeScores_NonPos_Pos,aes(Interactions,Score_NonPos_Pos,fill=Drug),stat='identity',alpha=1)+
           #geom_bar(stat='identity',alpha=1)+
           geom_hline(yintercept = 0)+
           scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
           labs(x = "Model interactions", y = 'Druggable (negative) \nversus \nnon-druggable (non-negative) \nscores')
        #geom_rect(fill = 'red', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)
        #remove background, rotate labels
        b_NonPos_Pos=b_NonPos_Pos + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
        #print(b_NonpositiveScores)
        
        ggsave(b_NonPos_Pos, file=paste(drugScores_folder_of_current_script_settings_,'plot_x_Ints_y_AllDrugs_CoDrugNonCoDrugScores.pdf',sep=""), width = 20, height = 4.5)
        
        
        
        
        
        b_Scores_CoDrug_NonCoDrug = b_NonPos_Pos
        
     }
     
     
     
     
     
     
     defectiveScores=list("defectiveInts"=allDefective,
                          "numDefectiveInts" = numDefectiveDF,
                          'scorePlot'=b,'scorePlot_CoDrug_NonCoDrug'=b_Scores_CoDrug_NonCoDrug,
                          'numDefectivePlot'=a,
                          'healthyNW'=H,'MSuntreatedNw'=MSuntreatedNw,'allDrugNws_allInts'=allDrugNws,
                          'alwaysDefective'=alwaysDefective)

  
  }  # END CASE I: thisMode=='median'
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # *************************************************************************************************************************
  # *********** CASE II: thisMode=='mean'
  # *************************************************************************************************************************
  
  
  if(thisMode=='mean'){
     
     
     
     
     
     
     
     
     # *************************************************************************************************************************
     # *********** I.1. Calculate the distance between MS and Healthy for each edge
     # *************************************************************************************************************************
     IdxHealthy = annot$Group == 'Healthy'
     IdxMSuntreated = (annot$Treatment == "no" & annot$Disease.Subtype!='')
     
     H = phenotypeNetwork(IdxHealthy, allMedianNetworks,model,mode=thisMode,linkActivityQuantileThreshold=linkActivityQuantileThreshold_used,applylinkActivityQuantileThreshold=applylinkActivityQuantileThreshold_used)  
     thisPhenotype='Healthy'
     write.table(H$network,file=paste0(phenotypeNws_folder_of_current_script_settings_,thisPhenotype,"_network.csv"),sep=",",row.names=T)
     
     MSuntreatedNw=phenotypeNetwork(IdxMSuntreated, allMedianNetworks,model,mode=thisMode,linkActivityQuantileThreshold=linkActivityQuantileThreshold_used,applylinkActivityQuantileThreshold=applylinkActivityQuantileThreshold_used) 
     thisPhenotype='MS'
     write.table(MSuntreatedNw$network,file=paste0(phenotypeNws_folder_of_current_script_settings_,thisPhenotype,"_network.csv"),sep=",",row.names=T)
     
     
     #distance
     HD=H$network-MSuntreatedNw$network
     
     absHD = abs(HD)
     thisPhenotype = "AbsDiff_Healthy_MSuntreated"
     write.table(absHD,file=paste0(phenotypeNws_folder_of_current_script_settings_,thisPhenotype,"_network.csv"),sep=",",row.names=T)
     
     
     
     # *************************************************************************************************************************
     # *********** I.2. Identify drug defective interactions, which have to fullfill the following 2 conditions:
     #            
     #            Condition C1: The difference MSuntreated to healthy is smaller or equal than drug to healthy.
     #            Condition C2: Healthy is different from drug.
     #
     # *************************************************************************************************************************
     
     
     
     #Gilenya=Fingolimod, Glatiramer=Copaxone
     phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
     numDrugs=5
     phenotypes_colors <- c('IFNb'='#612158', 'Tysabri'='#7A6FAC', 'Gilenya'='#57AB27', 'Copaxone'='#BDCD00', 'EGCG'='#0098A1')
     
     #prepare dataframe for plot
     phenotypeScores=data.frame(Interactions=colnames(allMedianNetworks),Gilenya=rep(NA,numInteractions),
                                IFNb=rep(NA,numInteractions),Copaxone=rep(NA,numInteractions),EGCG=rep(NA,numInteractions),Tysabri=rep(NA,numInteractions))
     phenotypeScores=melt(phenotypeScores,id='Interactions',value.name='Score')
     phenotypeScores$Score_ints0AndOk_fixed_to_1 = phenotypeScores$Score
     phenotypeScores$CoDruggable = phenotypeScores$Score
     colnames(phenotypeScores)=c('Interactions','Drug','Score','Score_ints0AndOk_fixed_to_1','CoDruggable')
     
     numDrugDefective=vector(length=numDrugs)
     names(numDrugDefective)=phenotypes
     
     allDrugNws_allInts=list()
     allDrugNws_StrongActiveInts=list()
     
     allDefective=list()
     
     

     for (i in 1:numDrugs){
        thisPhenotype=phenotypes[i]
        patientsThisPhenotype= annot$Treatment==thisPhenotype   #which are the patients in this phenotype?
        #**************Create phenotype network for Drug i
        phenoNw=phenotypeNetwork(patientsThisPhenotype,allMedianNetworks,model,mode=thisMode,linkActivityQuantileThreshold=linkActivityQuantileThreshold_used,applylinkActivityQuantileThreshold=applylinkActivityQuantileThreshold_used) 
        #**************Calculate the distance between i drug-Network and Healthy network
        H2Pheno=H$network-phenoNw$network 
        
        
        
        
        # restrict model network to strong active model interactions
        # based on the quantile
        
        
        StrongActiveInts <- union(names(which(H$network > quantile(H$network, probs = linkActivityQuantileThreshold_used))),
                                  union(names(which(MSuntreatedNw$network > quantile(MSuntreatedNw$network, probs = linkActivityQuantileThreshold_used))),
                                        names(which(phenoNw$network > quantile(phenoNw$network, probs = linkActivityQuantileThreshold_used)))))
        
        
        
        HD_StrongActiveInts = HD[StrongActiveInts]
        H2Pheno_StrongActiveInts = H2Pheno[StrongActiveInts]
        H2Pheno_StrongActiveInts_mat = as.matrix(H2Pheno_StrongActiveInts,ncol=1)
        
        H2Pheno_StrongActiveInts_RemovedIntsNA = NA*H2Pheno
        H2Pheno_StrongActiveInts_RemovedIntsNA[StrongActiveInts] = H2Pheno_StrongActiveInts
        
        
        #**************Substract the above from Healthy to MS
        drugScores=abs(HD_StrongActiveInts)-abs(H2Pheno_StrongActiveInts) 
        drugScores_withoutEpsilonRegion = drugScores
        
        # adjust zero-scores for the Epsilon-Region in which the drugScore is in the range (0 - epsilon, 0 + epsilon)
        DrugScoreEpsilonRegion_value <- quantile(abs(MSuntreatedNw$network[StrongActiveInts] - phenoNw$network[StrongActiveInts]), probs=DrugScoreEpsilonRegionQuantileThreshold)
        
        drugScores_withEpsilonRegion = drugScores
        drugScores_withEpsilonRegion[drugScores_withEpsilonRegion > - DrugScoreEpsilonRegion_value & drugScores_withEpsilonRegion < DrugScoreEpsilonRegion_value] <- 0
        
        
        
        
        #phenotypeScores[,i+1]=drugScores
        drugScores_withEpsilonRegion_allInts = NA*HD
        drugScores_withEpsilonRegion_allInts[StrongActiveInts] = drugScores_withEpsilonRegion
        
        
        drugScores_withoutEpsilonRegion_allInts = NA*HD
        drugScores_withoutEpsilonRegion_allInts[StrongActiveInts] = drugScores_withoutEpsilonRegion
        
        drugScores_withEpsilonRegion_allInts_mat = as.matrix(drugScores_withEpsilonRegion_allInts,ncol=1)
        drugScores_withoutEpsilonRegion_allInts_mat = as.matrix(drugScores_withoutEpsilonRegion_allInts,ncol=1)
        
        
        write.table(drugScores_withoutEpsilonRegion,file=paste0(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_withoutEpsilonRegion.csv"),sep=",",row.names=T)
        save(drugScores_withoutEpsilonRegion,file=paste(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_withoutEpsilonRegion.RData",sep=""))
        
        
        
        write.table(drugScores_withEpsilonRegion,file=paste0(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_withEpsilonRegion.csv"),sep=",",row.names=T)
        save(drugScores_withEpsilonRegion,file=paste(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_withEpsilonRegion.RData",sep=""))
        
        
        write.table(drugScores_withEpsilonRegion_allInts,file=paste0(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_withEpsilonRegion_allInts.csv"),sep=",",row.names=T)
        save(drugScores_withEpsilonRegion_allInts,file=paste(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_withEpsilonRegion_allInts.RData",sep=""))
        
        
        
        
        #**************Fill in dataframe for plotting with scores
        phenotypeScores$Score[which(as.character(phenotypeScores$Drug)==thisPhenotype)]=drugScores_withEpsilonRegion_allInts
        #barplot(drugScores,ylab='drugScores',xlab='Interactions (original network:186)',main=paste0(thisPhenotype,". Defective:",length(which(drugScores<0))))
        
        
        
        
        
        #************** Use the quantile function for the condition C2: Healthy is different from drug. 
        
        DiffThreshold_abs_H_minus_D_value <- quantile(abs(H2Pheno_StrongActiveInts), probs= DiffQuantileThreshold_abs_H_minus_D)
        
        
        
        #**************Identify and save defective Ints. Important: there are non-druggable 0 scores,i.e. where situation is ok.
        if(druggable=="negative"){
           
           defectiveInts = drugScores_withEpsilonRegion_allInts[(which(drugScores_withEpsilonRegion_allInts<0 & abs(H2Pheno_StrongActiveInts_RemovedIntsNA) > DiffThreshold_abs_H_minus_D_value))]  # C1 AND C2 need to be fulfilled
           allDefective[[i]]=defectiveInts
           
           phenotypeScores$CoDruggable[which(as.character(phenotypeScores$Drug)==thisPhenotype)] = 0
           phenotypeScores$CoDruggable[which(as.character(phenotypeScores$Drug)==thisPhenotype)][unname((which(drugScores_withEpsilonRegion_allInts<0 & abs(H2Pheno_StrongActiveInts_RemovedIntsNA) > DiffThreshold_abs_H_minus_D_value)))] = 1
           
           
           #**************Plot DrugScore versus Diff_absHminusDrug
           
           drug_color <- phenotypes_colors[thisPhenotype]
           
           
           df_DrugScore_Diff_absHminusDrug <- data.frame(score=drugScores_withEpsilonRegion,
                                                      diff_absHminusDrug=abs(H2Pheno_StrongActiveInts)
           )
           
           df_DrugScore_Diff_absHminusDrug$defectiveInts = (df_DrugScore_Diff_absHminusDrug$score < 0 & df_DrugScore_Diff_absHminusDrug$diff_absHminusDrug > DiffThreshold_abs_H_minus_D_value)
           
           c_plot_x_DrugScore_y_Diff_absHminusDrug = ggplot(df_DrugScore_Diff_absHminusDrug, aes(x=score, y=diff_absHminusDrug)) + 
              geom_point(size=2,col=ifelse(df_DrugScore_Diff_absHminusDrug$defectiveInts, drug_color, 'gray50')) + 
              xlab('Score') +
              ylab('abs( H - Drug)') +
              # xlim(xlim_DrugScore__min,xlim_DrugScore__max)+
              theme_bw() + 
              geom_hline(yintercept = DiffThreshold_abs_H_minus_D_value , linetype="dashed")+
              geom_vline(xintercept = 0 , linetype="dashed")
           
           ggsave(c_plot_x_DrugScore_y_Diff_absHminusDrug, file=paste(drugScores_folder_of_current_script_settings_,thisPhenotype,'__x_Score_y_Diff_absHminusDrug.pdf',sep=""), width = 5, height = 3)
           
           
           
           
           
        } else if(druggable=="zero"){
           
           #In reactions with 0-score or score values within the Epsilon+Region (almost 0) due to situation ok, replace by 1. 
           # use names as ids to fulfill two conditions in vectors scores and H2Pheno without messing up
           
           intsSituationNotOk=names(which(H2Pheno_StrongActiveInts <= - DrugScoreEpsilonRegion_value |  H2Pheno_StrongActiveInts >= DrugScoreEpsilonRegion_value))
           intsSituationOk=names(which(H2Pheno_StrongActiveInts > - DrugScoreEpsilonRegion_value & H2Pheno_StrongActiveInts < DrugScoreEpsilonRegion_value))
           
           ints0AndNotOk=intersect(names(drugScores_withEpsilonRegion[which(drugScores_withEpsilonRegion==0)]), intsSituationNotOk)
           ints0AndOk=intersect(names(drugScores_withEpsilonRegion[which(drugScores_withEpsilonRegion==0)]), intsSituationOk)
           
           intsNeg=which(drugScores_withEpsilonRegion<0)
           intsPos=which(drugScores_withEpsilonRegion>0)
           warning(paste0("Neg: ",length(intsNeg),"\nPos: ",length(intsPos),"\n0AndOk: ",length(ints0AndOk),"\n0andNotOk: ",length(ints0AndNotOk),"\n"))
           # replace
           drugScores_withEpsilonRegion[which(names(drugScores_withEpsilonRegion) %in% ints0AndOk)]=1
           drugScores_withEpsilonRegion_allInts[which(names(drugScores_withEpsilonRegion_allInts) %in% ints0AndOk)]=1
           defectiveInts = drugScores_withEpsilonRegion_allInts[(which(drugScores_withEpsilonRegion_allInts<=0 & abs(H2Pheno_StrongActiveInts_RemovedIntsNA) > DiffThreshold_abs_H_minus_D_value))]  # C1 AND C2 need be fulfilled
           allDefective[[i]]=defectiveInts
           
           
           phenotypeScores$Score_ints0AndOk_fixed_to_1[which(as.character(phenotypeScores$Drug)==thisPhenotype)]=drugScores_withEpsilonRegion_allInts
           phenotypeScores$CoDruggable[which(as.character(phenotypeScores$Drug)==thisPhenotype)] = 0
           phenotypeScores$CoDruggable[which(as.character(phenotypeScores$Drug)==thisPhenotype)][unname((which(drugScores_withEpsilonRegion_allInts<=0 & abs(H2Pheno_StrongActiveInts_RemovedIntsNA) > DiffThreshold_abs_H_minus_D_value)))] = 1
           
           
           
           
           
           #**************Plot DrugScore versus Diff_absHminusDrug
           
           drug_color <- phenotypes_colors[thisPhenotype]
           
           
           df_DrugScore_Diff_absHminusDrug <- data.frame(score=drugScores_withEpsilonRegion,
                                                      diff_absHminusDrug=abs(H2Pheno_StrongActiveInts)
                                                      )
           
           df_DrugScore_Diff_absHminusDrug$defectiveInts = (df_DrugScore_Diff_absHminusDrug$score <= 0 & df_DrugScore_Diff_absHminusDrug$diff_absHminusDrug > DiffThreshold_abs_H_minus_D_value)
           
           c_plot_x_DrugScore_y_Diff_absHminusDrug = ggplot(df_DrugScore_Diff_absHminusDrug, aes(x=score, y=diff_absHminusDrug)) + 
              geom_point(size=2,col=ifelse(df_DrugScore_Diff_absHminusDrug$defectiveInts, drug_color, 'gray50')) + 
              xlab('Score') +
              ylab('abs( H - Drug)') +
              # xlim(xlim_DrugScore__min,xlim_DrugScore__max)+
              theme_bw() + 
              geom_hline(yintercept = DiffThreshold_abs_H_minus_D_value , linetype="dashed")+
              geom_vline(xintercept = 0 , linetype="dashed")
           
           ggsave(c_plot_x_DrugScore_y_Diff_absHminusDrug, file=paste(drugScores_folder_of_current_script_settings_,thisPhenotype,'__x_Score_y_Diff_absHminusDrug.pdf',sep=""), width = 5, height = 3)
           
           
           
           
           
           
           
           
        }
        
        #**************Calculate number of defective
        numDrugDefective[i]=length(defectiveInts)
        warning(paste0(numDrugDefective[i]," druggable reactions for ",thisPhenotype,"\n"))
        
        
        

 
        
        
        
        
        
        
        
        
        
        
        #**************Save this network to structure in order to concatenate all networks
        allDrugNws_allInts[[i]]=phenotypeNetwork(patientsThisPhenotype,allMedianNetworks,model,mode=thisMode,linkActivityQuantileThreshold=linkActivityQuantileThreshold_used,applylinkActivityQuantileThreshold=applylinkActivityQuantileThreshold_used) 
        allDrugNws_allInts[[i]]$Drug=thisPhenotype
        
        
        allDrugNws_StrongActiveInts[[i]]=phenotypeNetwork(patientsThisPhenotype,allMedianNetworks[,StrongActiveInts],model,mode=thisMode,linkActivityQuantileThreshold=linkActivityQuantileThreshold_used,applylinkActivityQuantileThreshold=applylinkActivityQuantileThreshold_used) 
        allDrugNws_StrongActiveInts[[i]]$Drug=thisPhenotype
        

        
        #**************save individual files: for this drug, (i) drugScores_withEpsilonRegion, (ii) defective reactions and (iii) network
        
        write.table(drugScores_withEpsilonRegion_allInts,file=paste0(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_withEpsilonRegion_allInts__modified_ints0AndOk_fixed_to_1.csv"),sep=",",row.names=T)
        write.table(drugScores_withEpsilonRegion,file=paste0(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_withEpsilonRegion__modified_ints0AndOk_fixed_to_1.csv"),sep=",",row.names=T)
        
        write.table(defectiveInts,file=paste0(drugScores_folder_of_current_script_settings_,thisPhenotype,"_drugScores_withEpsilonRegion_defectiveInts.csv"),sep=",",row.names=T)
        
        write.table(allDrugNws_allInts[[i]]$network,file=paste0(phenotypeNws_folder_of_current_script_settings_,thisPhenotype,"_allInts_network.csv"),sep=",",row.names=T)
        write.table(allDrugNws_StrongActiveInts[[i]]$network,file=paste0(phenotypeNws_folder_of_current_script_settings_,thisPhenotype,"_StrongActiveInts_network.csv"),sep=",",row.names=T)
        
        
     }
     
     
     #****** which are the interactions that are defective in all drugs?
     alwaysDefective=intersect(intersect(intersect(intersect(names(allDefective[[1]]),names(allDefective[[2]])),names(allDefective[[3]])),names(allDefective[[4]])),names(allDefective[[5]]))
     write.table(alwaysDefective,file=paste0(drugScores_folder_of_current_script_settings_,"alwaysDefective.csv"),sep=",",row.names=T)
     
     # *************************************************************************************************************************
     # *********** I.3. Visualize defective interactions and phenotype networks
     # *************************************************************************************************************************
     
     #*******plot num of defective for each drug
     numDefectiveDF=data.frame(Drug=character(length=numDrugs),NumDefective=vector(length=numDrugs))
     numDefectiveDF$Drug=phenotypes
     for(i in 1:numDrugs){numDefectiveDF$NumDefective[i]=length(allDefective[[i]])}
     #To prevent ggplot2 from reordering alphabetically the labels according to the factors, which are alphabetical
     #We specify that the factors need to be ordered as they already are
     numDefectiveDF$Drug=factor(numDefectiveDF$Drug,levels=numDefectiveDF$Drug)
     a=ggplot(data=numDefectiveDF,aes(Drug,NumDefective,fill=Drug))+
        geom_bar(stat='identity',alpha=1)+
        scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE) +
        theme_bw()+
        theme(axis.text=element_text(size=10),axis.title=element_text(size=10))
     
     ggsave(a, file=paste(drugScores_folder_of_current_script_settings_,'plot_NumDefective.pdf',sep=""), width = 4, height = 3)
     
     
     #   #*******plot histogram of num defective by drug
     #   c=ggplot(data=phenotypeScores,aes(Score,fill=Drug))+geom_histogram(position='dodge',alpha=.5)
     #   c=c +theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
     #   
     #   #*******plot density of num defective by drug
     #   c=ggplot(data=phenotypeScores,aes(Score,fill=Drug))+geom_density(alpha=.5)
     #   c=c +theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
     #   
     
     #******plot bars with superposed intearctions to identify common defective
     phenotypeScores_count = phenotypeScores
     phenotypeScores_count$ScoreCount = phenotypeScores$Score
     phenotypeScores_count$ScoreCount[which(phenotypeScores_count$Score > 0)] = 1
     phenotypeScores_count$ScoreCount[which(phenotypeScores_count$Score < 0)] = -1
     
     b=ggplot(data=phenotypeScores_count,aes(Interactions,ScoreCount,fill=Drug))+
        geom_bar(stat='identity',alpha=1)+
        geom_hline(yintercept = 0)+
        scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)
     #remove background, rotate labels
     b=b + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
     #NOT necessary as this infor is in new plot a. add title including total of defective ints by drug
     # b=b+ggtitle(paste0('Drug-defective Interactions. ',phenotypes[1],':',numDrugDefective[1],"  ",phenotypes[2],':',numDrugDefective[2],"  ",phenotypes[3],':',numDrugDefective[3],"  ",
     #                      phenotypes[4],':',numDrugDefective[4],"  ",phenotypes[5],':',numDrugDefective[5],"  "))
     
     ggsave(b, file=paste(drugScores_folder_of_current_script_settings_,'plot_x_Ints_y_AllDrugsScores.pdf',sep=""), width = 20, height = 4.5)
     
     
     
     #******plot bars with superposed intearctions to identify common defective only for StrongActiveInts
     phenotypeScores_count_StrongActiveInts = phenotypeScores_count[complete.cases(phenotypeScores_count$Score), ]
     
     
     b_StrongActiveInts=ggplot(data=phenotypeScores_count_StrongActiveInts,aes(Interactions,ScoreCount,fill=Drug))+
        geom_bar(stat='identity',alpha=1)+
        geom_hline(yintercept = 0)+
        scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)
     #remove background, rotate labels
     b_StrongActiveInts=b_StrongActiveInts + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
     
     ggsave(b_StrongActiveInts, file=paste(drugScores_folder_of_current_script_settings_,'plot_x_StrongActiveInts_y_AllDrugsScores.pdf',sep=""), width = 8.5, height = 4)
     
     
     
     
     
     phenotypeScores_CoDrug_NonCoDrug = phenotypeScores
     phenotypeScores_CoDrug_NonCoDrug$Score_CoDrug_NonCoDrug = phenotypeScores$CoDruggable
     phenotypeScores_CoDrug_NonCoDrug$Score_CoDrug_NonCoDrug[which(phenotypeScores$CoDruggable == 1)] = -1
     phenotypeScores_CoDrug_NonCoDrug$Score_CoDrug_NonCoDrug[which(phenotypeScores$CoDruggable == 0)] = 1
     
     b_CoDrug_NonCoDrug=ggplot()+
        geom_bar(data=phenotypeScores_CoDrug_NonCoDrug,aes(Interactions,Score_CoDrug_NonCoDrug,fill=Drug),stat='identity',alpha=1)+
        scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
        geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0),
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="orangered", alpha=0.5)+
        geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf),
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray50", alpha=0.25)+
        geom_bar(data=phenotypeScores_CoDrug_NonCoDrug,aes(Interactions,Score_CoDrug_NonCoDrug,fill=Drug),stat='identity',alpha=1)+
        geom_hline(yintercept = 0)+
        scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
        labs(x = "Model interactions", y = 'Druggable (negative) \nversus \nnon-druggable (non-negative) \nscores')
     #geom_rect(fill = 'red', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)
     #remove background, rotate labels
     b_CoDrug_NonCoDrug=b_CoDrug_NonCoDrug + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
     #print(b_CoDrug_NonCoDrug)
     
     ggsave(b_CoDrug_NonCoDrug, file=paste(drugScores_folder_of_current_script_settings_,'plot_x_Ints_y_AllDrugs_CoDrugNonCoDrugScores.pdf',sep=""), width = 20, height = 4.5)
     
     
     
     b_Scores_CoDrug_NonCoDrug = b_CoDrug_NonCoDrug
     
     
     
     
     
     
     phenotypeScores_CoDrug_NonCoDrug_StrongActiveInts = phenotypeScores_CoDrug_NonCoDrug[complete.cases(phenotypeScores_CoDrug_NonCoDrug$Score), ]
     
     b_CoDrug_NonCoDrug_StrongActiveInts=ggplot()+
        geom_bar(data=phenotypeScores_CoDrug_NonCoDrug_StrongActiveInts,aes(Interactions,Score_CoDrug_NonCoDrug,fill=Drug),stat='identity',alpha=1)+
        scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
        geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=0),
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="orangered", alpha=0.5)+
        geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf),
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray50", alpha=0.25)+
        geom_bar(data=phenotypeScores_CoDrug_NonCoDrug_StrongActiveInts,aes(Interactions,Score_CoDrug_NonCoDrug,fill=Drug),stat='identity',alpha=1)+
        geom_hline(yintercept = 0)+
        scale_fill_manual(values=c("#57AB27","#612158","#BDCD00" ,"#0098A1","#7A6FAC"), guide=FALSE)+
        labs(x = "Model interactions", y = 'Druggable (negative) \nversus \nnon-druggable (non-negative) \nscores')
     #geom_rect(fill = 'red', xmin = -Inf, xmax = Inf, ymin =0, ymax = Inf, alpha =0.05)
     #remove background, rotate labels
     b_CoDrug_NonCoDrug_StrongActiveInts=b_CoDrug_NonCoDrug_StrongActiveInts + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
     #print(b_CoDrug_NonCoDrug_StrongActiveInts)
     
     ggsave(b_CoDrug_NonCoDrug_StrongActiveInts, file=paste(drugScores_folder_of_current_script_settings_,'plot_x_Ints_y_AllDrugs_CoDrugNonCoDrugScores_StrongActiveInts.pdf',sep=""), width = 10.5, height = 3.5)
     
     
     
     b_Scores_CoDrug_NonCoDrug_StrongActiveInts = b_CoDrug_NonCoDrug_StrongActiveInts
     
     
     
     
     
     
     
     
    
     
    
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     defectiveScores=list("defectiveInts"=allDefective,
                          "numDefectiveInts" = numDefectiveDF,
                          'scorePlot'=b,'scorePlot_StrongActiveInts'=b_StrongActiveInts,'scorePlot_CoDrug_NonCoDrug'=b_Scores_CoDrug_NonCoDrug,'scorePlot_CoDrug_NonCoDrug_StrongActiveInts'=b_Scores_CoDrug_NonCoDrug_StrongActiveInts,
                          'numDefectivePlot'=a,
                          'healthyNW'=H,'MSuntreatedNw'=MSuntreatedNw,'allDrugNws_allInts'=allDrugNws_allInts,'allDrugNws_StrongActiveInts'=allDrugNws_StrongActiveInts,
                          'alwaysDefective'=alwaysDefective)
     
     
  }  # END CASE II: thisMode=='mean'
  
  
  
  
  
  
  
  
  
  
  
  
  save(defectiveScores,file=paste(drugScores_folder_of_current_script_settings_,"defectiveScores.RData",sep=""))
  
  
  
  
  
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


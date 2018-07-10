####################################################################################
## PLOT FIGURE PREDICTION AND VALIDATION OF COMBINATION THERAPIES
##
## Panel A: Drug combination scores
## Panel B: FTY network with co-druggable reactions, to be inserted later with Inkscape
## Panel C: Validation FTY-EGCG
## Panel D: Validation FTY-TAK1i
## 
####################################################################################

# ####################################################################################
# User defined input variables
# ####################################################################################

input__median_models_used = "files/median_models/median_models__based_on__SIFcombiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED"
input__Results_based_on_SIF_used = 'Results_based_on_SIF_combiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED'
model_used =   "combiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED.sif"                      # MR inserted
#model_used =   "combiMSplaneCUT.sif"                      # MR inserted






## Load Packages
library(CellNOptR) # Version 1.16
library(reshape2) # Version 1.4.1
library(ggplot2) # Version 2.1.0
# library(cowplot) # Version 0.6.2 # loaded later, otherwise it messes up the facet grid plot
library(gdata) # Version 2.17.0

## Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ####################################################################################
# Panel A # mostly copied from code/combination_therapy_prediction/calculateDefective.R
# ####################################################################################

# Function to calculate and visualize co-druggable reactions (previously termed defective). Called
# by pathDrugTargetsFinalv3.R
# Marti Bernardo-Faura
# Final version April 2016

# Minor changes for Github repository by Jakob Wirbel, July 2017
source("../combination_therapy_prediction/phenotypeNetwork.R")
#load('../../files/median_models/allMedianModels.RData')
file_path__allMedianModels_RData=  paste('../../',input__median_models_used,'/allMedianModels.RData',sep="")
load(file_path__allMedianModels_RData)  # called allMedianNetworks


thisMode='median'
# phenotypeNws_folder='../../files/group_models/'
# drugScores_folder='../../files/drugScores/'
phenotypeNws_folder=paste('../../files/group_models/',input__Results_based_on_SIF_used,sep="")
drugScores_folder=paste('../../files/drugScores/',input__Results_based_on_SIF_used,sep="")


# ************load anotation to map patients to groups
  #   data_folder="/Users/marti/Documents/ebi/combiMS/data/phosphosMergedAbsMax/processed/normalized/secondRoundProcessedMidas/"
  #   filenames=list.files(data_folder,pattern="*.csv",full.names=FALSE)
annot=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
  #   numPat=length(filenames)
  #   filenames2=filenames
  #   for (j in 1:numPat){
  #     filenames2[j]=strsplit(filenames[j],"\\.")[[1]][1]
  #   }
  
# ************load model and midas for annotation
# 


#model_path="../../files/model/combiMSplaneCUT.sif"
model_path=paste("../../files/model/",model_used,sep='')

model=readSIF(model_path)  
numInteractions=length(model$reacID)

# ***********load predicted networks
# *********** 1. Calculate the distance between MS and Healthy for each edge

IdxHealthy = annot$Group == 'Healthy'
IdxMSuntreated = (annot$Treatment == "no" & annot$Disease.Subtype!='')
  
H = phenotypeNetwork(IdxHealthy, allMedianNetworks,model,mode=thisMode)  
  
MSuntreatedNw=phenotypeNetwork(IdxMSuntreated, allMedianNetworks,model,mode=thisMode)

HD=H$network-MSuntreatedNw$network

# *********** 2. Identify edges where difference MS to health is smaller or equal than drug to Health
# *********** These are drug defective interactions, hence calling for combination therapy that will revert signaling

#Gilenya=Fingolimod, Glatiramer=Copaxone
phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
numDrugs=5

#prepare dataframe for plot
phenotypeScores=data.frame(Interactions=colnames(allMedianNetworks),Gilenya=rep(NA,numInteractions),
                           IFNb=rep(NA,numInteractions),Copaxone=rep(NA,numInteractions),EGCG=rep(NA,numInteractions),Tysabri=rep(NA,numInteractions))
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
  #**************Substract the above from Healthy to Ms
  drugScores=abs(HD)-abs(H2Pheno) 
  
  #**************Fill in dataframe for plotting with scores
  phenotypeScores$Score[which(as.character(phenotypeScores$Drug)==thisPhenotype)]=drugScores
    
  #**************Identify and save defective Ints. Important: there are non-drugable 0 scores,i.e. where situation is ok.
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
    
  #**************Calculate number of defective
  numDrugDefective[i]=length(defectiveInts)
  warning(paste0(numDrugDefective[i]," drugable reactions for ",thisPhenotype,"\n"))
    
  #**************Save this network to structure in order to concatenate all networks
  allDrugNws[[i]]=phenotypeNetwork(patientsThisPhenotype,allMedianNetworks,model,mode=thisMode)
  allDrugNws[[i]]$Drug=thisPhenotype
  
}
  
  
#****** which are the interactions that are defective in all drugs?
alwaysDefective=intersect(intersect(intersect(intersect(names(allDefective[[1]]),names(allDefective[[2]])),names(allDefective[[3]])),names(allDefective[[4]])),names(allDefective[[5]]))
  
# *********** 3. Visualize defective interactions and phenotype networks
#*******plot num of defective for each drug
numDefectiveDF=data.frame(Drug=character(length=numDrugs),NumDefective=vector(length=numDrugs))
numDefectiveDF$Drug=phenotypes
for(i in 1:numDrugs){numDefectiveDF$NumDefective[i]=length(allDefective[[i]])}

#To prevent ggplot2 from reordering alphabetically the labels according to the factors, which are alphabetical
#We specify that the factors need to be ordered as they already are
# numDefectiveDF$Drug=factor(numDefectiveDF$Drug,levels=numDefectiveDF$Drug)
# manual re-naming
numDefectiveDF$drug_name = c('FTY', 'IFNb', 'GA', 'EGCG', 'NTZ')
numDefectiveDF$drug_name = factor(numDefectiveDF$drug_name, levels=c('IFNb', 'NTZ', 'FTY', 'GA', 'EGCG'))
# Barplot with the number of defective interaction for each drug, serves as legend for the main part of the figure
a = ggplot(data=numDefectiveDF,
           aes(drug_name,NumDefective,fill=drug_name)) + 
  geom_bar(stat='identity') + 
  theme_classic() + theme(legend.position = 'none') + 
  xlab('') + ylab('Number of \n co-druggable reactions') +
  theme(axis.title.y = element_text(size=8)) +
  scale_fill_manual(values=c('#612158', '#7A6FAC', '#57AB27', '#BDCD00', '#0098A1'), guide=FALSE) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

#******plot bars with superposed intearctions to identify common defective
phenotypeScores$Drug = factor(phenotypeScores$Drug, levels=c('IFNb', 'Tysabri', 'Gilenya', 'Copaxone', 'EGCG'))
b = ggplot(data=phenotypeScores, 
           aes(Interactions,Score,fill=Drug)) + 
  geom_bar(stat='identity') +
  theme_classic() + theme(legend.position = 'none') + 
  xlab('Model Interactions') + ylab('Score') + 
  scale_fill_manual(values=c('#612158', '#7A6FAC', '#57AB27', '#BDCD00', '#0098A1'), guide=FALSE) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

# ####################################################################################
# Panel C
# ####################################################################################

## get data
plot_df_c = read.xls('../../data/validation_experiments/EAE_FTY_EGCG.xls', sheet=1, header=TRUE)
colnames(plot_df_c) = c('days', 'FTY', 'Placebo', 'FTY+\nEGCG', 'EGCG')
plot_df_c = melt(plot_df_c, id='days')
colnames(plot_df_c) = c('days', 'treatment', 'clinical score')

## reorder treatment
plot_df_c$treatment = factor(plot_df_c$treatment, levels=c('Placebo', 'FTY', 'EGCG', 'TAKi', 'FTY+\nEGCG', 'FTY+\nTAKi'))
plot_df_c$comb = 'FTY + EGCG'

# Plot
g_c = ggplot(plot_df_c, aes(x=days, y=`clinical score`, color=treatment)) +
  geom_line() + theme_classic()  + 
  geom_point(shape=15) + theme(strip.background=element_blank(), legend.key.size = unit(1.3, 'lines')) + 
  scale_color_manual(values=c('#9C9E9F', '#57AB27', '#0098A1',  '#006165', '#CC071E', '#A11035'), 
                     labels=levels(plot_df_c$treatment), drop=FALSE)

# ####################################################################################
# Panel C
# ####################################################################################

## get data
plot_df_d = read.xls('../../data/validation_experiments/EAE_FTY_TAK1i.xls', sheet=1, header=TRUE)[,c(1,2,3,4,5)]
colnames(plot_df_d) = c('days', 'FTY+\nTAKi', 'FTY', 'Placebo', 'TAKi')
plot_df_d = melt(plot_df_d, id='days')
colnames(plot_df_d) = c('days', 'treatment', 'clinical score')

## reorder treatment
plot_df_d$treatment = factor(plot_df_d$treatment, levels=c('Placebo', 'FTY', 'TAKi', 'FTY+\nTAKi'))
plot_df_d$comb = 'FTY + TAKi'

# Plot
g_d = ggplot(plot_df_d, aes(x=days, y=`clinical score`, color=treatment)) +
  geom_line() + theme_classic()  + 
  geom_point(shape=15) + theme(strip.background=element_blank(), legend.key.size = unit(1.2, 'lines')) +
  scale_color_manual(values=c('#9C9E9F', '#57AB27', '#006165', '#A11035'), 
                     drop=FALSE, guide=FALSE)

# ####################################################################################
# Combine everything
# ####################################################################################

library(cowplot) # Version 0.6.2 # loaded later, otherwise it messes up the facet grid plot

pdf('../../figures/figure_combination_therapy.pdf', width=7, height=7)
top_row = plot_grid(b, a, labels=c('A', 'B'), nrow=1, rel_widths = c(1, .45))
middle_row = plot_grid(NULL, labels = c('C'))
bottom_row = plot_grid(g_c, g_d, labels=c('D', 'E'), nrow=1, rel_widths = c(1.4, 1))
plot_grid(top_row, middle_row, bottom_row, nrow=3)
dev.off()


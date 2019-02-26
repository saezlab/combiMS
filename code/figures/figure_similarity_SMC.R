
# figure_similarity_SMC



# # ###############################################################################################################################################
# Based on script: analyzePatientSimilarity.R by Marti July 2015
#
# Minor modification for nicer plotting by Jakob in December 2016
#
# Some modification by Melanie in December 2017 and in December 2018
# #################################################################################################################################################



Results_storage_name = "Results_based_on_SIF_combiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED"           
Results_storage_name__median_models = "median_models__based_on__SIFcombiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED"    

#include__similarityMatrixSelf__of_each_donor = include__similarityMatrixSelf__of_each_donor                                 
include__similarityMatrixSelf__of_each_donor = TRUE #FALSE                                 





library(apcluster)
library(gclus)
library(ggplot2)
library(reshape2)
library(NMF)
library(corrplot)
library(cowplot)
library(gdata)
library(ggsignif)




# Use relative paths instead of absolut paths for the files
# in Rstudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))





figures_folder = "../../figures"   


ifelse(!dir.exists(file.path(figures_folder,"figure_similarity_SMC")), dir.create(file.path(figures_folder,"figure_similarity_SMC")), FALSE)              
sub_dir__Figures_folder_of_current_script_class = file.path(figures_folder,"figure_similarity_SMC")                                                                              

ifelse(!dir.exists(file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)), dir.create(file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)), FALSE)              
sub_dir__Figures_folder_of_current_script = file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)                                                                              


ifelse(!dir.exists(file.path(sub_dir__Figures_folder_of_current_script,"Unique_similarityMatrix")), dir.create(file.path(sub_dir__Figures_folder_of_current_script,"Unique_similarityMatrix")), FALSE)              
sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix = file.path(sub_dir__Figures_folder_of_current_script,"Unique_similarityMatrix")     


# ###########################################################################################################
# Load and clean annotation data
# ###########################################################################################################
annot=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)
annot169pat_v2=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)


# ************************************ 
# add some fields to annotation for ggplot




# correct that PPMS patients appear in group instead of by untreated and treatment



#simplify category
annot2=annot169pat_v2
annot2$Category[which(annot$Category=="dont know")]='PPMS'
annot2$Category[which(annot2$Category=='PPMS slow')]='PPMS'
annot2$Category[which(annot2$Category=='PPMS fast')]='PPMS'

#test if all 169 are correctly labeled
length(which(annot2$Category=='PPMS')) + length(which(annot2$Category=='RRMS')) + length(which(annot2$Category=='healthy'))
length(which(annot2$Category=='RRMS' & annot2$condition=='Untreated')) + length(which(annot2$Category=='PPMS' & annot2$condition=='Untreated')) + length(which(annot2$condition=='Treated')) + length(which(annot2$condition=='Healthy'))




Idx_RRMSuntreated= which(annot2$Category=='RRMS' & annot2$condition=='Untreated')
DonorIDs_RRMSuntreated = annot2$ID[Idx_RRMSuntreated]

Idx_PPMSuntreated= which(annot2$Category=='PPMS' & annot2$condition=='Untreated')
DonorIDs_PPMSuntreated = annot2$ID[Idx_PPMSuntreated]






IdxHealthy = annot2$Group == 'Healthy'
DonorIDs_Healthy = annot2$ID[IdxHealthy]

IdxMSuntreated = (annot2$Treatment == "no" & annot2$Disease.Subtype!='')
DonorIDs_MSuntreated = annot2$ID[IdxMSuntreated]





phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')


IdxDrug_Gilenya= annot2$Treatment=='Gilenya' 
DonorIDs_Gilenya = annot2$ID[IdxDrug_Gilenya]

IdxDrug_IFNb= annot2$Treatment=='IFNb' 
DonorIDs_IFNb = annot2$ID[IdxDrug_IFNb]

IdxDrug_Copaxone= annot2$Treatment=='Copaxone' 
DonorIDs_Copaxone = annot2$ID[IdxDrug_Copaxone]

IdxDrug_EGCG= annot2$Treatment=='EGCG' 
DonorIDs_EGCG = annot2$ID[IdxDrug_EGCG]

IdxDrug_Tysabri= annot2$Treatment=='Tysabri' 
DonorIDs_Tysabri = annot2$ID[IdxDrug_Tysabri]


# 
# # Rename groups for consistency with other figures in the paper
# 
# annot2_Group_unique = unique(annot_figure$Group)                   
# annot2_Group_unique
# #[1] "Glatiramer"   "Healthy"      "Untreated"    "PPMS"         "Interferon B" "Fingolimod"   "EGCG"         "Natalizumab" 
# 
# annot2_Disease.Subtype_unique = unique(annot_figure$Disease.Subtype)                   
# annot2_Disease.Subtype_unique
# # [1] "RR"  ""    "CIS" "PP"  "SP"  "PR" 
# 


annot3 = annot2





annot3_Group_unique = unique(annot3$Group)                   
annot3_Group_unique
#[1] "Glatiramer"   "Healthy"      "Untreated"    "PPMS"         "Interferon B" "Fingolimod"   "EGCG"         "Natalizumab" 


annot3[which(annot3$Group=='PPMS' & annot3$Treatment=='no'),'Group']='Untreated'                                                                # 1. PPMS correction
annot3[which(annot3$Group=='PPMS' & annot3$Treatment!='no'),'Group']=annot3[which(annot3$Group=='PPMS' & annot3$Treatment!='no'),'Treatment']   # 2. PPMS correction
annot3[which(annot3$Group=='Gilenya'),'Group']='FTY'                    # Gilenya  is the trade name of Fingolimod=FTY

annot3[which(annot3$Group=='Fingolimod'),'Group']='FTY'
annot3[which(annot3$Group=='Glatiramer'),'Group']='GA'
annot3[which(annot3$Group=='Interferon B'),'Group']='IFNb'
annot3[which(annot3$Group=='Natalizumab'),'Group']='NTZ'

annot3[which(annot3$Group=='EGCG'),'Group']=annot3[which(annot3$Group=='EGCG'),'Treatment']   # EGCG correction

annot3_Group_unique = unique(annot3$Group)                   
annot3_Group_unique
# > annot3_Group_unique
# [1] "GA"                 "Healthy"            "Untreated"          "IFNb"               "FTY"               
# [6] "EGCG"               "IFNb + EGCG"        "EGCG or Placebo"    "NTZ"                "Tysabri or Placebo"




annot3_Disease.Subtype_unique = unique(annot3$Disease.Subtype)                   
annot3_Disease.Subtype_unique
#> annot3_Disease.Subtype_unique
# [1] "RR"  ""    "CIS" "PP"  "SP"  "PR" 

annot3[which(annot3$Disease.Subtype == ''), 'Disease.Subtype'] = 'healthy'


annot3[which(annot3$Disease.Subtype == 'PR' & annot3$Category=='PPMS'), 'Disease.Subtype'] = 'PP'
annot3[which(annot3$Disease.Subtype == 'PR' & annot3$Category=='RRMS'), 'Disease.Subtype'] = 'RR'

annot3[which(annot3$Disease.Subtype == 'SP' & annot3$Category=='PPMS'), 'Disease.Subtype'] = 'PP'
annot3[which(annot3$Disease.Subtype == 'SP' & annot3$Category=='RRMS'), 'Disease.Subtype'] = 'RR'

annot3_Disease.Subtype_unique = unique(annot3$Disease.Subtype)                   
annot3_Disease.Subtype_unique
# > annot3_Disease.Subtype_unique
# [1] "RR"      "Healthy" "CIS"     "PP"
# 
# 
# annot3[which(annot3$Disease.Subtype == 'SP'), 'Disease.Subtype'] = substring(annot3[which(annot3$Disease.Subtype == 'SP'), 'Category'], 1, 2)
# annot3[which(annot3$Disease.Subtype == 'PR'), 'Disease.Subtype'] = substring(annot3[which(annot3$Disease.Subtype == 'PR'), 'Category'], 1, 2)


annot=annot3
rm(annot3)

#annot_figure[which(annot_figure$Group=='PPMS' & annot_figure$Treatment=='no'),'Group']='Untreated'                                                                # 1. PPMS correction
#annot_figure[which(annot_figure$Group=='PPMS' & annot_figure$Treatment!='no'),'Group']=annot_figure[which(annot_figure$Group=='PPMS' & annot_figure$Treatment!='no'),'Treatment']   # 2. PPMS correction
# annot_figure[which(annot_figure$Group=='Gilenya'),'Group']='FTY'                    # Gilenya  is the trade name of Fingolimod=FTY
# 
# annot_figure[which(annot_figure$Group=='Fingolimod'),'Group']='FTY'
# annot_figure[which(annot_figure$Group=='Glatiramer'),'Group']='GA'
# annot_figure[which(annot_figure$Group=='Interferon B'),'Group']='IFNb'
# annot_figure[which(annot_figure$Group=='Natalizumab'),'Group']='NTZ'
# 
# 
# annot_figure[which(annot_figure$Disease.Subtype == ''), 'Disease.Subtype'] = 'healthy'
# annot_figure[which(annot_figure$Disease.Subtype == 'SP'), 'Disease.Subtype'] = substring(annot_figure[which(annot_figure$Disease.Subtype == 'SP'), 'Category'], 1, 2)
# annot_figure[which(annot_figure$Disease.Subtype == 'PR'), 'Disease.Subtype'] = substring(annot_figure[which(annot_figure$Disease.Subtype == 'PR'), 'Category'], 1, 2)
# annot=annot_figure
# rm(annot_figure)




# test consistency of the clinical donor data within the annotation file -------------------------------------------------




Idx_RR_CIS= which(annot$Disease.Subtype=='CIS')
DonorIDs_RR_CIS = annot$ID[Idx_RR_CIS]

Idx_RRuntreated= which(annot$Disease.Subtype=='RR' & annot$condition=='Untreated')
DonorIDs_RRuntreated = annot$ID[Idx_RRuntreated]

Idx_RR_CISuntreated= which(annot$Disease.Subtype=='CIS' & annot$condition=='Untreated')
DonorIDs_RR_CISuntreated = annot$ID[Idx_RR_CISuntreated]

Idx_PPuntreated= which(annot$Disease.Subtype=='PP' & annot$condition=='Untreated')
DonorIDs_PPuntreated = annot$ID[Idx_PPuntreated]


IdxHealthy_subtype = (annot$Disease.Subtype=='healthy')
DonorIDs_Healthy_subtype = annot$ID[IdxHealthy_subtype]

IdxHealthy_condition = (annot$condition=='Healthy')
DonorIDs_Healthy_condition = annot$ID[IdxHealthy_condition]

IdxMSuntreated_condition = (annot$condition=='Untreated')
DonorIDs_MSuntreated_condition = annot$ID[IdxMSuntreated_condition]

IdxMStreated_condition = (annot$condition=='Treated')
DonorIDs_MStreated_condition = annot$ID[IdxMStreated_condition]




phenotypes=c('Gilenya','IFNb','Copaxone','EGCG','Tysabri')
phenotypes_abbreviation=c('FTY','IFNb','GA','EGCG','NTZ')



IdxDrug_FTY= annot$Group=='FTY' 
DonorIDs_FTY = annot$ID[IdxDrug_FTY]

IdxDrug_IFNb= annot$Group=='IFNb' 
DonorIDs_IFNb = annot$ID[IdxDrug_IFNb]

IdxDrug_GA= annot$Group=='GA' 
DonorIDs_GA = annot$ID[IdxDrug_GA]

IdxDrug_EGCG= annot$Group=='EGCG' 
DonorIDs_EGCG = annot$ID[IdxDrug_EGCG]

IdxDrug_NTZ= annot$Group=='NTZ' 
DonorIDs_NTZ = annot$ID[IdxDrug_NTZ]





annot_Group_unique = unique(annot$Group)                   
annot_Group_unique

# > annot_Group_unique   # 10
# [1] "GA"                 "Healthy"            "Untreated"          "IFNb"               "FTY"               
# [6] "EGCG"               "IFNb + EGCG"        "EGCG or Placebo"    "NTZ"                "Tysabri or Placebo"



annot_Disease.Subtype_unique = unique(annot$Disease.Subtype)                   
annot_Disease.Subtype_unique   
# > annot_Disease.Subtype_unique   # 4
# [1] "RR"      "Healthy" "CIS"     "PP"     






# ##########################################################################################################
# PART A OF FIGURE
# ##########################################################################################################

#load('../../files/median_models/allMedianModels.RData')        


files_median_models__folder = "../../files/median_models"                                       
#sub_dir__Results_allMedianModels = file.path(files_median_models__folder,Results_storage_name__median_models)     
sub_dir__Results_allMedianModels = file.path(files_median_models__folder)     

load(file=paste0(sub_dir__Results_allMedianModels,"/allMedianModels.RData"))   





# Merge annotation data (twice, so that it is broader) with model interactions
plot_df = melt(allMedianNetworks)

temp = melt(annot, id.vars = 'ID')
# > 169*9
# [1] 1521 obs. in temp
temp_original = temp                  

temp = temp[which(temp$variable %in% c('Center', 'condition', 'Disease.Subtype')),]
# > 169*3
# [1] 507
# 
 

temp2 = temp
temp2$variable = as.character(temp2$variable)
temp2[which(temp2$variable == 'Center'),'variable'] = 'Center2'
temp2[which(temp2$variable == 'condition'),'variable'] = 'condition2'
temp2[which(temp2$variable == 'Disease.Subtype'),'variable'] = 'Disease.Subtype2'
temp2$variable = factor(temp2$variable)

# Variable for facet grid splitting
plot_df$split= 4
temp$split = 1
temp2$split = 1

temp[which(temp$variable == 'Disease.Subtype'), 'split'] = 3
temp2[which(temp2$variable == 'Disease.Subtype2'), 'split'] = 3

temp[which(temp$variable == 'condition'), 'split'] = 2
temp2[which(temp2$variable == 'condition2'), 'split'] = 2



colnames(plot_df) = c('Patients', 'Model_Interactions', 'value', 'split')
colnames(temp) = c('Patients', 'Model_Interactions', 'value', 'split')
colnames(temp2) = c('Patients', 'Model_Interactions', 'value', 'split')

# strings instead of numeric entries for factor geom_tile
plot_df = rbind(plot_df, temp, temp2)

plot_df[which(plot_df$value == 0),'value'] = 'Inactive'
plot_df[which(plot_df$value == 1),'value'] = 'Active'
plot_df[which(plot_df$value == .5),'value'] = 'Active'

# Reorder factor so that the annotation is broader

plot_df_Model_Interactions_unique = unique(plot_df$Model_Interactions)            
plot_df_Model_Interactions_unique                                                 
# ...
# [166] !LCK=PLCG1          SH3BP2=RAC1         VITD3=STAT1         Center              Disease.Subtype    
# [171] condition           Center2             Disease.Subtype2    condition2         
# 180 Levels: AKT1=CREB1 AKT1=GSK3A !AKT1=RAF1 !AKT1=SLP76 AKT1=STAT1 AKT1=STAT3 !AKT1=ZAP70 ... Disease.Subtype2


NumModelInteractions = as.numeric(ncol(allMedianNetworks))

plot_df$Model_Interactions = factor(plot_df$Model_Interactions, levels = c(levels(plot_df$Model_Interactions)[1:NumModelInteractions], 'Center', 'Center2', 'Disease.Subtype', 'Disease.Subtype2', 'condition', 'condition2'))

# reorder for legend

plot_df_value_unique = unique(plot_df$value)            
plot_df_value_unique                                                 
# > plot_df_value_unique
# [1] "Inactive"  "Active"    "CH"        "IB"        "KI"        "UZ"        "RR"        "Healthy"   "CIS"       "PP"        "Treated"   "Untreated"


# > levels(plot_df$value)
# NULL

plot_df$value = factor(plot_df$value, levels = c('Center', 'CH', 'KI', 'IB', 'UZ', '',
                                                 'Condition', 'Healthy', 'Treated', 'Untreated', ' ',
                                                 'Subtype', 'healthy', 'CIS', 'RR', 'PP', '  ',
                                                 'Reaction', 'Active', 'Inactive'))







# #########################################################################################################################
# Plot as Ggplot heatmap
# #########################################################################################################################

FigureWidth_map =5     
FigureHeight_map=8     


NumAllDonors = as.numeric(length(annot$ID))

heatmap_median_models = ggplot(plot_df, aes(x=Model_Interactions, y=Patients)) + 
   facet_grid(~split, scales = "free_x", space='free', margins=FALSE) + 
   geom_tile(aes(fill=value)) + 
   xlab(paste('Model Reactions (1-',NumModelInteractions,')',sep="")) +
   ylab(paste('Patients (1-', NumAllDonors , ')',sep="")) +
   theme(strip.background=element_blank(), strip.text= element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
         axis.ticks = element_blank(), axis.line.y = element_blank(), legend.text = element_text(size=6), axis.line.x = element_blank(),
         axis.title = element_text(size=9),
         panel.spacing.x = unit(1, 'pt')) + 
   scale_fill_manual(values=c("white", '#2D7F83', '#B65256', '#D0D95C', '#9B91C1', 'white',
                              "white", '#407FB7', '#8DC060', '#FABE50','white',
                              "white", '#407FB7', '#E69679', '#89CCCF','#A8859E', 'white',
                              "white", 'black', 'grey90'),
                     drop=FALSE) +
   guides(fill=guide_legend(ncol=1, keywidth = .6, keyheight = .6, title=''))
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__A__heatmap_median_models.pdf",sep = "")), width = FigureWidth_map, height = FigureHeight_map)             





























# ##########################################################################################################
# PART B OF FIGURE
# ##########################################################################################################

# ##########################################################################################################
# load annotation file and similarity BETWEEN matrix
# ####################################################################################################################



files_Similarity__folder = "../../files/similarity"                                                             
sub_dir__Results_Similarity = file.path(files_Similarity__folder,Results_storage_name)                          

load(paste(sub_dir__Results_Similarity,"/similarityMatrix.RData",sep="")) # called similarityMatrix               



# ************************************ 
# include similarities WITHIN




# load('../../files/similarity/similarityWITHINMeans.RData')
# for (i in 1:dim(similarityMatrix)[1]){
#    similarityMatrix[i,i]=means_within[i]
# }


if(include__similarityMatrixSelf__of_each_donor==TRUE){                                   
   
   #load('/Users/marti/Documents/R/combiMS/cluster/analysis/similarityWITHIN.RData')               
   load(paste(sub_dir__Results_Similarity,"mean_SimilarityBetweenDonorOwnModels.RData",sep="/")) # called mean_SimilarityBetweenDonorOwnModels               
   
   subset_Self_mean = mean(mean_SimilarityBetweenDonorOwnModels,na.rm = TRUE)
   subset_Self_median = median(mean_SimilarityBetweenDonorOwnModels,na.rm = TRUE)
   
   
   
   #colnames(similarityMatrix)==names(similarityVector)
   for (i in 1:dim(similarityMatrix)[1]){
      #similarityMatrix[i,i]=similarityVector[i]
      similarityMatrix[i,i]=mean_SimilarityBetweenDonorOwnModels[i]             
   }
}


# if(include__similarityMatrixSelf__of_each_donor==TRUE){                                   
#   
#   #load('/Users/marti/Documents/R/combiMS/cluster/analysis/similarityWITHIN.RData')               
#   load(paste(sub_dir__Results_Similarity,"mean_SimilarityBetweenDonorOwnModels.RData",sep="/")) # called mean_SimilarityBetweenDonorOwnModels               
#   
#   #colnames(similarityMatrix)==names(similarityVector)
#   for (i in 1:dim(similarityMatrix)[1]){
#     #similarityMatrix[i,i]=similarityVector[i]
#     similarityMatrix[i,i]=mean_SimilarityBetweenDonorOwnModels[i]             
#   }
# }












similarityMatrix_unique = similarityMatrix
similarityMatrix_unique[lower.tri(similarityMatrix_unique,diag=F)]=NA                      


similarityMatrix_unique_NoDiag = similarityMatrix
similarityMatrix_unique_NoDiag[lower.tri(similarityMatrix_unique,diag=T)]=NA                      





# ##########################################################################################################
# create dataframe of similarities to plot
# ##########################################################################################################
sim_df=melt(similarityMatrix_unique)
sim_df$Group=rep(F,dim(similarityMatrix_unique)[1])
sim_df$GroupSim=rep(F,dim(similarityMatrix_unique)[1])
sim_df$Condition=rep(F,dim(similarityMatrix_unique)[1])
sim_df$ConditionSim=rep(F,dim(similarityMatrix_unique)[1])
sim_df$Self=rep(F,dim(similarityMatrix_unique)[1])
sim_df$RR_Untreated=rep(F,dim(similarityMatrix_unique)[1])

colnames(sim_df)=c('SimWith','Patient','Similarity','Group','GroupSim','Condition','ConditionSim','Self','RR_Untreated')  
sim_df$GroupSim=rep(annot$Group,length(annot$ID))
sim_df$ConditionSim=rep(annot$condition,length(annot$ID))
sim_df$all="All"

for (i in 1:length(annot$ID)){
   sim_df$Group[sim_df$Patient==annot[i,'ID']] = annot[i,'Group']
   sim_df$Condition[sim_df$Patient==annot[i,'ID']] = annot[i,'condition'] 
}

pb = txtProgressBar(min=0, max=length(sim_df$Patient))

for (i in 1:length(sim_df$Patient)){
   if (as.character(sim_df$Patient[i])==as.character(sim_df$SimWith[i])){
      sim_df$Self[i]='Self'
   } else {
      sim_df$Self[i]='Another'
   }
   
   if (as.character(sim_df$Patient[i]) %in% annot[which(annot$Disease.Subtype=='RR' & annot$Treatment=='no'),'ID']){
      sim_df$RR_Untreated[i]='RR Untreated'
   } else {
      sim_df$RR_Untreated[i]='N'
   }
   setTxtProgressBar(pb, (pb$getVal()+1))
}

# ##########################################################################################################
# plot distribution of subsets where similarities are subsetting only between patients with same conditions
# and against all
# ############################################################################################################

#sub_same_group_similarity=subset(sim_df,sim_df$Group==sim_df$GroupSim)
#sub_same_group_similarity__WithoutSelf=subset(sim_df,sim_df$Group==sim_df$GroupSim & sim_df$Self == 'Another')
sub_same_group_similarity__WithSelf=subset(sim_df,sim_df$Group==sim_df$GroupSim)
sub_same_group_similarity=subset(sub_same_group_similarity__WithSelf,sub_same_group_similarity__WithSelf$Self == 'Another')


sub_same_group_similarity___Group_unique_names = unique(sub_same_group_similarity$Group)
sub_same_group_similarity___Group_unique_names
# > sub_same_group_similarity___Group_unique_names
# [1] "GA"              "Healthy"         "Untreated"       "IFNb"            "FTY"             "EGCG"            "EGCG or Placebo" "NTZ"  


sub_self_similarity=subset(sub_same_group_similarity__WithSelf,sub_same_group_similarity__WithSelf$Self=='Self')


sub_RRuntreat=subset(sub_same_group_similarity,sub_same_group_similarity$RR_Untreated=='RR Untreated')




# # delete "IFNb + EGCG"   --- not necessary as ONLY in SELF
# sub_same_group_similarity = sub_same_group_similarity[-which(sub_same_group_similarity$Group == 'IFNb + EGCG'),]

# delete "EGCG or Placebo"
sub_same_group_similarity = sub_same_group_similarity[-which(sub_same_group_similarity$Group == "EGCG or Placebo"),]

# # delete Tysabri or Placebo  --- not necessary as ONLY in SELF
# sub_same_group_similarity = sub_same_group_similarity[-which(sub_same_group_similarity$Group == 'Tysabri or Placebo'),]




sub_same_group_similarity___Group_unique_names = unique(sub_same_group_similarity$Group)
sub_same_group_similarity___Group_unique_names
# > sub_same_group_similarity___Group_unique_names
# [1] "GA"        "Healthy"   "Untreated" "IFNb"      "FTY"       "EGCG"      "NTZ"    





# add All similarities, self similarities, and RR_untreated as reference
sub_self_similarity$Group = 'Self'
sub_RRuntreat$Group = 'RR Untreated'
sub_all = sim_df
sub_all$Group = 'All'

# Create additional variable for grouping of boxplots
sub_same_group_similarity$Grid = 'Treatment'
sub_same_group_similarity$Grid[which(sub_same_group_similarity$Group == 'Healthy')] = 'Status'
sub_same_group_similarity$Grid[which(sub_same_group_similarity$Group == 'Untreated')] = 'Status'
sub_self_similarity$Grid = 'All'
sub_RRuntreat$Grid = 'Status'
sub_all$Grid = 'All'







# boxplot by Group
# In a notched box plot, the notches extend 1.58 * IQR / sqrt(n). This gives a roughly 95 interval for comparing medians
# So non-overlapping notches strongly suggest that medians are significantly different
# This is the case for RR untreated vis RR IFNb

plot_df_2 = rbind(sub_same_group_similarity, sub_self_similarity, sub_RRuntreat, sub_all)


plot_df_2_Group__unique = unique(plot_df_2$Group)                 
# plot_df_2_Group__unique
# [1] "GA"           "Healthy"      "Untreated"    "IFNb"         "FTY"          "EGCG"         "NTZ"          "Self"         "RR Untreated" "All"     

## Reorder Groups
sorted_groups = c('All', 'Self', 'Healthy', 'EGCG', 'FTY', 'IFNb', 'GA', 'NTZ', 'RR Untreated', 'Untreated')
plot_df_2$Group = factor(plot_df_2$Group, sorted_groups)

FigureWidth_boxplot = 6   
FigureHeight_boxplot = 4  


if(include__similarityMatrixSelf__of_each_donor==TRUE){
   color_vec_group_fill = c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')
}else{
   color_vec_group_fill = c('#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')
}


# boxplot_similarity = ggplot(plot_df_2, aes(x=Group, y=Similarity, fill=Group)) + 
#    #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
#    geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         
#    theme_classic() + 
#    #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
#    scale_fill_manual(values=color_vec_group_fill) +                        
#    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
#    theme(legend.position = 'none') + 
#    #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
#    #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                                                 
#    #geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
#    geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1,na.rm = TRUE) +                                                 
#    geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.05,na.rm = TRUE) +                                                     
#    ylim(0.55,1.1)               
# ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity__unique.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             


boxplot_similarity = ggplot(plot_df_2, aes(x=Group, y=Similarity, fill=Group)) + 
   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
   geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         
   theme_classic() + 
   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
   scale_fill_manual(values=color_vec_group_fill) +                        
   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
   theme(legend.position = 'none') + 
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                                                 
   #geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
   geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 0.985,na.rm = TRUE) +                                                 
   geom_signif(comparisons=list(c("FTY", "RR Untreated")), map_signif_level = TRUE, y_position = 1.025,na.rm = TRUE) +                                                     
   ylim(0.59,1.03)               
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity__unique.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             


























boxplot_similarity_test_all_drugs__versus__untreated = ggplot(plot_df_2, aes(x=Group, y=Similarity, fill=Group)) + 
   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
   geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         
   theme_classic() + 
   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
   scale_fill_manual(values=color_vec_group_fill) +                        
   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
   theme(legend.position = 'none') + 
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                                                 
   #geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
   geom_signif(comparisons=list(c("EGCG", "Untreated")), map_signif_level = TRUE, y_position = 1.0,na.rm = TRUE) +                                                 
   geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.1,na.rm = TRUE) +                                                     
   geom_signif(comparisons=list(c("IFNb", "Untreated")), map_signif_level = TRUE, y_position =1.2,na.rm = TRUE) +                                                 
   geom_signif(comparisons=list(c("GA", "Untreated")), map_signif_level = TRUE, y_position = 1.3,na.rm = TRUE) +  
   geom_signif(comparisons=list(c("NTZ", "Untreated")), map_signif_level = TRUE, y_position = 1.4,na.rm = TRUE) +                                                     
   
   ylim(0.55,1.1)              
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity_test_all_drugs__versus__untreated__unique.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             

boxplot_similarity_test_all_drugs__versus__RRuntreated = ggplot(plot_df_2, aes(x=Group, y=Similarity, fill=Group)) + 
   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
   geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         
   theme_classic() + 
   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
   scale_fill_manual(values=color_vec_group_fill) +                        
   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
   theme(legend.position = 'none') + 
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                                                 
   #geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
   geom_signif(comparisons=list(c("EGCG", "RR Untreated")), map_signif_level = TRUE, y_position = 1.0,na.rm = TRUE) +                                                 
   geom_signif(comparisons=list(c("FTY", "RR Untreated")), map_signif_level = TRUE, y_position = 1.1,na.rm = TRUE) +                                                     
   geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.2,na.rm = TRUE) +                                                 
   geom_signif(comparisons=list(c("GA", "RR Untreated")), map_signif_level = TRUE, y_position = 1.3,na.rm = TRUE) +  
   geom_signif(comparisons=list(c("NTZ", "RR Untreated")), map_signif_level = TRUE, y_position = 1.4,na.rm = TRUE) +                                                     
   
   ylim(0.55,1.1)        
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity_test_all_drugs__versus__RRuntreated__unique.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             














plot_df_2__subset_without_NA = subset(plot_df_2, !is.na(plot_df_2$Similarity))
plot_df_2__subset_without_NA_Group_unique = unique(plot_df_2__subset_without_NA$Group)
plot_df_2__subset_without_NA_Group_unique
# [1] Untreated    IFNb         EGCG         GA           FTY          Healthy      NTZ          RR Untreated All         
# Levels: All Self Healthy EGCG FTY IFNb GA NTZ RR Untreated Untreated

boxplot_similarity_subset_without_NA = ggplot(plot_df_2__subset_without_NA, aes(x=Group, y=Similarity, fill=Group)) + 
   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
   geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         
   theme_classic() + 
   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
   scale_fill_manual(values=color_vec_group_fill) +                        
   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
   theme(legend.position = 'none') + 
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                                                 
   #geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
   geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05,na.rm = TRUE) +                                                 
   geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15,na.rm = TRUE) +                                                     
   ylim(0.55,1.1)               
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity__unique__using_plot_df_2__subset_without_NA.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             



# ###################################################################################################################################################
# Divert output into pdf
# ###################################################################################################################################################

#pdf("../../figures/figure_similarity.pdf", width=7, height=3.5, onefile = FALSE)            
pdf(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__A__heatmap_median_models__B__boxplot_similarity.pdf",sep = "")), width=7, height=3.5, onefile = FALSE)

plot_grid(heatmap_median_models, boxplot_similarity, labels=c('A', 'B'), label_size = 14)

dev.off()


















sub_plot_df_2_EGCG=subset(plot_df_2,plot_df_2$Group=='EGCG')


boxplot_similarity_EGCG = ggplot(sub_plot_df_2_EGCG, aes(x=Group, y=Similarity, fill=Group)) + 
   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
   geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         
   theme_classic() + 
   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
   scale_fill_manual(values=c('#8DC060')) +                        
   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
   theme(legend.position = 'none') + 
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
   # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         
   # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
   ylim(0.55,1.1)               
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity__unique_EGCG.pdf",sep = "")), width = FigureWidth_boxplot*0.2, height = FigureHeight_boxplot)             



boxplot_similarity_EGCG = ggplot(sub_plot_df_2_EGCG, aes(x=Group, y=Similarity, colour=SimWith)) + 
   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
   geom_point(position = position_jitter(width = 0.5))+
   theme_classic() + 
   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
   #scale_fill_manual(values=c('#8DC060')) +                        
   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
   #theme(legend.position = 'none') + 
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
   # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         
   # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
   ylim(0.55,1.1)               
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity__unique_EGCG__geom_point_colour_SimWith.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             




boxplot_similarity_EGCG_2 = ggplot(sub_plot_df_2_EGCG, aes(x=Group, y=Similarity, colour=Patient)) + 
   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
   geom_point(position = position_jitter(width = 0.5))+
   theme_classic() + 
   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
   #scale_fill_manual(values=c('#8DC060')) +                        
   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
   #theme(legend.position = 'none') + 
   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
   # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         
   # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
   ylim(0.55,1.1)               
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity__unique_EGCG__geom_point_colour_Patient.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             


mat__donor_combinations_EGCG = sub_plot_df_2_EGCG[,c("SimWith","Patient")]

#EGCG_unique_indices=unique(t(apply(mat__donor_combinations_EGCG, 1, sort)))
EGCG_unique_indices=unique(t(apply(mat__donor_combinations_EGCG, 1, sort,decreasing=T)))

#EGCG_unique_indices = mat__donor_combinations_EGCG[!duplicated(t(apply(mat__donor_combinations_EGCG, 1, sort))),]
row.names(EGCG_unique_indices)    
sub_plot_df_2_EGCG_unique = sub_plot_df_2_EGCG[row.names(EGCG_unique_indices), ]




# boxplot_similarity_EGCG_unique = ggplot(sub_plot_df_2_EGCG_unique, aes(x=Group, y=Similarity, fill=Group)) + 
#   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
#   geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         
#   theme_classic() + 
#   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
#   scale_fill_manual(values=c('#8DC060')) +                        
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
#   theme(legend.position = 'none') + 
#   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
#   # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         
#   # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
#   ylim(0.55,1.1)               
# ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity__unique_EGCG_unique.pdf",sep = "")), width = FigureWidth_boxplot*0.2, height = FigureHeight_boxplot)             



# boxplot_similarity_EGCG_unique = ggplot(sub_plot_df_2_EGCG_unique, aes(x=Group, y=Similarity, colour=SimWith)) + 
#   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
#   geom_point(position = position_jitter(width = 0.5))+
#   theme_classic() + 
#   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
#   #scale_fill_manual(values=c('#8DC060')) +                        
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
#   #theme(legend.position = 'none') + 
#   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
#   # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         
#   # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
#   ylim(0.55,1.1)               
# ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity__unique_EGCG_unique__geom_point_colour_SimWith.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             
# 



# boxplot_similarity_EGCG_unique_2 = ggplot(sub_plot_df_2_EGCG_unique, aes(x=Group, y=Similarity, colour=Patient)) + 
#   #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
#   geom_point(position = position_jitter(width = 0.5))+
#   theme_classic() + 
#   #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
#   #scale_fill_manual(values=c('#8DC060')) +                        
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
#   #theme(legend.position = 'none') + 
#   #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
#   # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         
#   # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
#   ylim(0.55,1.1)               
# ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_SMC__B__boxplot_similarity__unique_EGCG_unique__geom_point_colour_Patient.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             
# 
# 
# 
# 


# ##########################################################################################################
# 
# Calculate mean and median for each subgroup
# 
# 
# ############################################################################################################


unique(plot_df_2$Group)
# [1] GA           Healthy      Untreated    IFNb         FTY          EGCG         NTZ          Self         RR Untreated All         
# Levels: All Self Healthy EGCG FTY IFNb GA NTZ RR Untreated Untreated



subset_All=subset(plot_df_2,plot_df_2$Group=='All')
subset_All_mean = mean(subset_All$Similarity,na.rm=TRUE)
subset_All_median = median(subset_All$Similarity,na.rm=TRUE)







subset_SELF=subset(plot_df_2,plot_df_2$Group=='Self')
subset_SELF_mean = mean(subset_SELF$Similarity,na.rm=TRUE)
subset_SELF_median = median(subset_SELF$Similarity,na.rm=TRUE)


test_identical_subset_Self_mean = identical(subset_Self_mean,subset_SELF_mean)
test_identical_subset_Self_mean

test_identical_subset_Self_median = identical(subset_Self_median,subset_SELF_median)
test_identical_subset_Self_median









subset_Healthy=subset(plot_df_2,plot_df_2$Group=='Healthy' & plot_df_2$Self=='Another')
subset_Healthy_mean = mean(subset_Healthy$Similarity,na.rm=TRUE)
subset_Healthy_median = median(subset_Healthy$Similarity,na.rm=TRUE)




# [1] GA           Healthy      Untreated    IFNb         FTY          EGCG         NTZ          Self         RR Untreated All         
# Levels: All Self Healthy EGCG FTY IFNb GA NTZ RR Untreated Untreated
# 
# 
subset_EGCG=subset(plot_df_2,plot_df_2$Group=='EGCG' & plot_df_2$Self=='Another')
subset_EGCG_mean = mean(subset_EGCG$Similarity,na.rm=TRUE)
subset_EGCG_median = median(subset_EGCG$Similarity,na.rm=TRUE)

subset_GA=subset(plot_df_2,plot_df_2$Group=='GA' & plot_df_2$Self=='Another')
subset_GA_mean = mean(subset_GA$Similarity,na.rm=TRUE)
subset_GA_median = median(subset_GA$Similarity,na.rm=TRUE)


subset_IFNb=subset(plot_df_2,plot_df_2$Group=='IFNb' & plot_df_2$Self=='Another')
subset_IFNb_mean = mean(subset_IFNb$Similarity,na.rm=TRUE)
subset_IFNb_median = median(subset_IFNb$Similarity,na.rm=TRUE)


subset_FTY=subset(plot_df_2,plot_df_2$Group=='FTY' & plot_df_2$Self=='Another')
subset_FTY_mean = mean(subset_FTY$Similarity,na.rm=TRUE)
subset_FTY_median = median(subset_FTY$Similarity,na.rm=TRUE)


subset_NTZ=subset(plot_df_2,plot_df_2$Group=='NTZ' & plot_df_2$Self=='Another')
subset_NTZ_mean = mean(subset_NTZ$Similarity,na.rm=TRUE)
subset_NTZ_median = median(subset_NTZ$Similarity,na.rm=TRUE)





#RR Untreated Untreated



subset_RRuntreated=subset(plot_df_2,plot_df_2$Group=='RR Untreated' & plot_df_2$Self=='Another')
subset_RRuntreated_mean = mean(subset_RRuntreated$Similarity,na.rm=TRUE)
subset_RRuntreated_median = median(subset_RRuntreated$Similarity,na.rm=TRUE)





subset_Untreated=subset(plot_df_2,plot_df_2$Group=='Untreated' & plot_df_2$Self=='Another')
subset_Untreated_mean = mean(subset_Untreated$Similarity,na.rm=TRUE)
subset_Untreated_median = median(subset_Untreated$Similarity,na.rm=TRUE)



df_SMC_1mean_2median_values_for_each_subgroup = data.frame(Subgroup = c("All" , "Self" , "Healthy" ,"EGCG" , "FTY" , "IFNb", "GA" , "NTZ" , "RR Untreated" , "Untreated"),
                                                           Mean = c(subset_All_mean,subset_Self_mean,subset_Healthy_mean,
                                                                    subset_EGCG_mean,subset_FTY_mean,subset_IFNb_mean,subset_GA_mean,subset_NTZ_mean,
                                                                    subset_RRuntreated_mean,subset_Untreated_mean),
                                                           Median = c(subset_All_median,subset_Self_median,subset_Healthy_median,
                                                                      subset_EGCG_median,subset_FTY_median,subset_IFNb_median,subset_GA_median,subset_NTZ_median,
                                                                      subset_RRuntreated_median,subset_Untreated_median)
                                                           )

save(df_SMC_1mean_2median_values_for_each_subgroup, file=file.path(paste(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,"/df_SMC_1mean_2median_values_for_each_subgroup.RData",sep="")))            
write.csv(df_SMC_1mean_2median_values_for_each_subgroup, file=file.path(paste(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,"/df_SMC_1mean_2median_values_for_each_subgroup.csv",sep="")))            








# ##########################################################################################################
# 
# Perform wilcox.test
# between some subgroups 
# 
# ############################################################################################################




# (1) EGCG versus RRuntreated ------------------------------------------------------------------------------------------



# subset_EGCG=subset(plot_df_2,plot_df_2$Group=='EGCG' & plot_df_2$Self=='Another')
# subset_EGCG_mean = mean(subset_EGCG$Similarity,na.rm=TRUE)
# subset_EGCG_median = median(subset_EGCG$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_Drug = subset_EGCG$Similarity[!is.na(subset_EGCG$Similarity)]   


# subset_RRuntreated=subset(plot_df_2,plot_df_2$Group=='RR Untreated' & plot_df_2$Self=='Another')
# subset_RRuntreated_mean = mean(subset_RRuntreated$Similarity,na.rm=TRUE)
# subset_RRuntreated_median = median(subset_RRuntreated$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_RRUntreated =  subset_RRuntreated$Similarity[!is.na(subset_RRuntreated$Similarity)]   





WilcoxTest_subgroups_1Drug_2RRUntreated  = wilcox.test(WilcoxTest_subgroup_Drug,WilcoxTest_subgroup_RRUntreated,
                                                       paired=FALSE,
                                                       exact = T)

WilcoxTest_subgroups_1Drug_2RRUntreated__pValue = WilcoxTest_subgroups_1Drug_2RRUntreated$p.value
WilcoxTest_subgroups_1Drug_2RRUntreated__pValue


WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_EGCGversusRRuntreated = WilcoxTest_subgroups_1Drug_2RRUntreated__pValue


WilcoxTest_subgroups_1Drug_2RRUntreated
# Wilcoxon rank sum test with continuity correction
# 
# data:  WilcoxTest_subgroup_Drug and WilcoxTest_subgroup_RRUntreated
# W = 2107, p-value = 0.2098
# alternative hypothesis: true location shift is not equal to 0





# (2) FTY versus RRuntreated ------------------------------------------------------------------------------------------



# subset_FTY=subset(plot_df_2,plot_df_2$Group=='FTY' & plot_df_2$Self=='Another')
# subset_FTY_mean = mean(subset_FTY$Similarity,na.rm=TRUE)
# subset_FTY_median = median(subset_FTY$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_Drug = subset_FTY$Similarity[!is.na(subset_FTY$Similarity)]   
   
   
# subset_RRuntreated=subset(plot_df_2,plot_df_2$Group=='RR Untreated' & plot_df_2$Self=='Another')
# subset_RRuntreated_mean = mean(subset_RRuntreated$Similarity,na.rm=TRUE)
# subset_RRuntreated_median = median(subset_RRuntreated$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_RRUntreated =  subset_RRuntreated$Similarity[!is.na(subset_RRuntreated$Similarity)]   





WilcoxTest_subgroups_1Drug_2RRUntreated  = wilcox.test(WilcoxTest_subgroup_Drug,WilcoxTest_subgroup_RRUntreated,
                                           paired=FALSE,
                                           exact = T)

WilcoxTest_subgroups_1Drug_2RRUntreated__pValue = WilcoxTest_subgroups_1Drug_2RRUntreated$p.value
WilcoxTest_subgroups_1Drug_2RRUntreated__pValue


WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_FTYversusRRuntreated = WilcoxTest_subgroups_1Drug_2RRUntreated__pValue

WilcoxTest_subgroups_1Drug_2RRUntreated
# Wilcoxon rank sum test with continuity correction
# 
# data:  WilcoxTest_subgroup_Drug and WilcoxTest_subgroup_RRUntreated
# W = 46533, p-value = 0.002501
# alternative hypothesis: true location shift is not equal to 0





# (3) IFNb versus RRuntreated ------------------------------------------------------------------------------------------



# subset_IFNb=subset(plot_df_2,plot_df_2$Group=='IFNb' & plot_df_2$Self=='Another')
# subset_IFNb_mean = mean(subset_IFNb$Similarity,na.rm=TRUE)
# subset_IFNb_median = median(subset_IFNb$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_Drug = subset_IFNb$Similarity[!is.na(subset_IFNb$Similarity)]   


# subset_RRuntreated=subset(plot_df_2,plot_df_2$Group=='RR Untreated' & plot_df_2$Self=='Another')
# subset_RRuntreated_mean = mean(subset_RRuntreated$Similarity,na.rm=TRUE)
# subset_RRuntreated_median = median(subset_RRuntreated$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_RRUntreated =  subset_RRuntreated$Similarity[!is.na(subset_RRuntreated$Similarity)]   





WilcoxTest_subgroups_1Drug_2RRUntreated  = wilcox.test(WilcoxTest_subgroup_Drug,WilcoxTest_subgroup_RRUntreated,
                                                       paired=FALSE,
                                                       exact = T)

WilcoxTest_subgroups_1Drug_2RRUntreated__pValue = WilcoxTest_subgroups_1Drug_2RRUntreated$p.value
WilcoxTest_subgroups_1Drug_2RRUntreated__pValue


WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_IFNbversusRRuntreated = WilcoxTest_subgroups_1Drug_2RRUntreated__pValue


WilcoxTest_subgroups_1Drug_2RRUntreated
# Wilcoxon rank sum test with continuity correction
# 
# data:  WilcoxTest_subgroup_Drug and WilcoxTest_subgroup_RRUntreated
# W = 160610, p-value = 7.81e-06
# alternative hypothesis: true location shift is not equal to 0


# (4) GA versus RRuntreated ------------------------------------------------------------------------------------------



# subset_GA=subset(plot_df_2,plot_df_2$Group=='GA' & plot_df_2$Self=='Another')
# subset_GA_mean = mean(subset_GA$Similarity,na.rm=TRUE)
# subset_GA_median = median(subset_GA$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_Drug = subset_GA$Similarity[!is.na(subset_GA$Similarity)]   


# subset_RRuntreated=subset(plot_df_2,plot_df_2$Group=='RR Untreated' & plot_df_2$Self=='Another')
# subset_RRuntreated_mean = mean(subset_RRuntreated$Similarity,na.rm=TRUE)
# subset_RRuntreated_median = median(subset_RRuntreated$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_RRUntreated =  subset_RRuntreated$Similarity[!is.na(subset_RRuntreated$Similarity)]   





WilcoxTest_subgroups_1Drug_2RRUntreated  = wilcox.test(WilcoxTest_subgroup_Drug,WilcoxTest_subgroup_RRUntreated,
                                                       paired=FALSE,
                                                       exact = T)

WilcoxTest_subgroups_1Drug_2RRUntreated__pValue = WilcoxTest_subgroups_1Drug_2RRUntreated$p.value
WilcoxTest_subgroups_1Drug_2RRUntreated__pValue


WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_GAversusRRuntreated = WilcoxTest_subgroups_1Drug_2RRUntreated__pValue


WilcoxTest_subgroups_1Drug_2RRUntreated
# Wilcoxon rank sum test with continuity correction
# 
# data:  WilcoxTest_subgroup_Drug and WilcoxTest_subgroup_RRUntreated
# W = 24460, p-value = 0.2648
# alternative hypothesis: true location shift is not equal to 0





# (5) NTZ versus RRuntreated ------------------------------------------------------------------------------------------



# subset_NTZ=subset(plot_df_2,plot_df_2$Group=='NTZ' & plot_df_2$Self=='Another')
# subset_NTZ_mean = mean(subset_NTZ$Similarity,na.rm=TRUE)
# subset_NTZ_median = median(subset_NTZ$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_Drug = subset_NTZ$Similarity[!is.na(subset_NTZ$Similarity)]   


# subset_RRuntreated=subset(plot_df_2,plot_df_2$Group=='RR Untreated' & plot_df_2$Self=='Another')
# subset_RRuntreated_mean = mean(subset_RRuntreated$Similarity,na.rm=TRUE)
# subset_RRuntreated_median = median(subset_RRuntreated$Similarity,na.rm=TRUE)

WilcoxTest_subgroup_RRUntreated =  subset_RRuntreated$Similarity[!is.na(subset_RRuntreated$Similarity)]   





WilcoxTest_subgroups_1Drug_2RRUntreated  = wilcox.test(WilcoxTest_subgroup_Drug,WilcoxTest_subgroup_RRUntreated,
                                                       paired=FALSE,
                                                       exact = T)

WilcoxTest_subgroups_1Drug_2RRUntreated__pValue = WilcoxTest_subgroups_1Drug_2RRUntreated$p.value
WilcoxTest_subgroups_1Drug_2RRUntreated__pValue


WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_NTZversusRRuntreated = WilcoxTest_subgroups_1Drug_2RRUntreated__pValue


WilcoxTest_subgroups_1Drug_2RRUntreated
# Wilcoxon rank sum test with continuity correction
# 
# data:  WilcoxTest_subgroup_Drug and WilcoxTest_subgroup_RRUntreated
# W = 99129, p-value = 0.0003433
# alternative hypothesis: true location shift is not equal to 0







df_SMC_1mean_2median_3WilcoxTest_W1Drug_W2RRUntreated = data.frame(Subgroup = c("All" , "Self" , "Healthy" ,"EGCG" , "FTY" , "IFNb", "GA" , "NTZ" , "RR Untreated" , "Untreated"),
                                                           Mean = c(subset_All_mean,subset_Self_mean,subset_Healthy_mean,
                                                                    subset_EGCG_mean,subset_FTY_mean,subset_IFNb_mean,subset_GA_mean,subset_NTZ_mean,
                                                                    subset_RRuntreated_mean,subset_Untreated_mean),
                                                           Median = c(subset_All_median,subset_Self_median,subset_Healthy_median,
                                                                      subset_EGCG_median,subset_FTY_median,subset_IFNb_median,subset_GA_median,subset_NTZ_median,
                                                                      subset_RRuntreated_median,subset_Untreated_median),
                                                           WilcoxTest_W1Drug_W2RRUntreated = c(NaN , NaN , NaN ,
                                                                                               WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_EGCGversusRRuntreated,
                                                                                               WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_FTYversusRRuntreated,
                                                                                               WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_IFNbversusRRuntreated,
                                                                                               WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_GAversusRRuntreated,
                                                                                               WilcoxTest_subgroups_1Drug_2RRUntreated__pValue_NTZversusRRuntreated,
                                                                                               NaN , NaN))

save(df_SMC_1mean_2median_3WilcoxTest_W1Drug_W2RRUntreated, file=file.path(paste(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,"/df_SMC_1mean_2median_3WilcoxTest_W1Drug_W2RRUntreated.RData",sep="")))            
write.csv(df_SMC_1mean_2median_3WilcoxTest_W1Drug_W2RRUntreated, file=file.path(paste(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,"/df_SMC_1mean_2median_3WilcoxTest_W1Drug_W2RRUntreated.csv",sep="")))            


print("Script finished!")







# figure_similarity_MR

#figure_similarity_MR__using_Cluster_MR_Results_based_on_using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__gaBinaryT1_postCellNOptRupdate__OKpostCellNOptRupdate



# Some notes added
# 
# AND
# 
# paths adapted for CombiMS rerun
# 
# AND
# 
# some modifications with no influence on the results
# 
# AND
# 



# ERROR detected
# boxplot considers all values in the symmetric similarityMatrix  and 
# not only the unique ones of the similarityMatrix (diagonal entries and upper part (or lower) triangular part)
#
# Suggestion to fix this error:
# 
# (1) insert the following command:  
#     similarityMatrix[lower.tri(similarityMatrix,diag=F)]=NA                      # MR inserted
#     
# (2) use command na.rm within the functions
#     geom_boxplot and
#     geom_signif



# 
# by 
# Melanie Rinas
# December 2017




Results_storage_name = "OptCombiMSclusterFork10_MR__using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__gaBinaryT1_postCellNOptRupdate"        # MR inserted

include__similarityMatrixSelf__of_each_donor = FALSE                                 # MR inserted



# ###############################################################################################################################################
# Based on script: analyzePatientSimilarity.R by Marti July 2015
#
# Minor modification for nicer plotting by Jakob in December 2016
# #################################################################################################################################################

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





figures_MR_folder = "../../figures_MR"   


ifelse(!dir.exists(file.path(figures_MR_folder,"figure_similarity_MR")), dir.create(file.path(figures_MR_folder,"figure_similarity_MR")), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script_class = file.path(figures_MR_folder,"figure_similarity_MR")                                                                              # MR inserted

ifelse(!dir.exists(file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)), dir.create(file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script = file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)                                                                              # MR inserted

ifelse(!dir.exists(file.path(sub_dir__Figures_folder_of_current_script,"Whole_similarityMatrix")), dir.create(file.path(sub_dir__Figures_folder_of_current_script,"Whole_similarityMatrix")), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix = file.path(sub_dir__Figures_folder_of_current_script,"Whole_similarityMatrix")     

ifelse(!dir.exists(file.path(sub_dir__Figures_folder_of_current_script,"Unique_similarityMatrix")), dir.create(file.path(sub_dir__Figures_folder_of_current_script,"Unique_similarityMatrix")), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix = file.path(sub_dir__Figures_folder_of_current_script,"Unique_similarityMatrix")     


# ###########################################################################################################
# Load and clean annotation data
# ###########################################################################################################
annot=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)

# ************************************ 
# add some fields to annotation for ggplot
annot2=annot
# correct that PPMS patients appear in group instead of by untreated and treatment





# Rename groups for consistency with other figures in the paper

annot2_Group_unique = unique(annot2$Group)                   # MR inserted
annot2_Group_unique
# [1] "Glatiramer"   "Healthy"      "Untreated"    "PPMS"         "Interferon B" "Fingolimod"  
#  "EGCG"         "Natalizumab" 

annot2_Disease.Subtype_unique = unique(annot2$Disease.Subtype)                   # MR inserted
annot2_Disease.Subtype_unique
# [1] "RR"  ""    "CIS" "PP"  "SP"  "PR" 






annot2[which(annot2$Group=='PPMS' & annot2$Treatment=='no'),'Group']='Untreated'                                                                # 1. PPMS correction
annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Group']=annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Treatment']   # 2. PPMS correction
annot2[which(annot2$Group=='Gilenya'),'Group']='FTY'                    # Gilenya  is the trade name of Fingolimod=FTY

annot2[which(annot2$Group=='Fingolimod'),'Group']='FTY'
annot2[which(annot2$Group=='Glatiramer'),'Group']='GA'
annot2[which(annot2$Group=='Interferon B'),'Group']='IFNb'
annot2[which(annot2$Group=='Natalizumab'),'Group']='NTZ'

annot2[which(annot2$Disease.Subtype == ''), 'Disease.Subtype'] = 'healthy'
annot2[which(annot2$Disease.Subtype == 'SP'), 'Disease.Subtype'] = substring(annot2[which(annot2$Disease.Subtype == 'SP'), 'Category'], 1, 2)
annot2[which(annot2$Disease.Subtype == 'PR'), 'Disease.Subtype'] = substring(annot2[which(annot2$Disease.Subtype == 'PR'), 'Category'], 1, 2)
annot=annot2
rm(annot2)





annot_Group_unique = unique(annot$Group)                   # MR inserted
annot_Group_unique
# [1] "GA"                 "Healthy"            "Untreated"          "IFNb"               "FTY"               
# [6] "EGCG"               "NTZ"                "Tysabri or Placebo"


annot_Disease.Subtype_unique = unique(annot$Disease.Subtype)                   # MR inserted
annot_Disease.Subtype_unique
# [1] "RR"      "healthy" "CIS"     "PP"     




# ##########################################################################################################
# PART A OF FIGURE
# ##########################################################################################################

#load('../../files/median_models/allMedianModels.RData')        # MR modified

files_median_models_MR_folder = "../../files/median_models_MR"                                       # MR inserted
sub_dir__Results_allMedianModels = file.path(files_median_models_MR_folder,Results_storage_name)     # MR inserted

load(file=paste0(sub_dir__Results_allMedianModels,"/allMedianModels.RData"))   # MR modified



# Merge annotation data (twice, so that it is broader) with model interactions
plot_df = melt(allMedianNetworks)

temp = melt(annot, id.vars = 'ID')
# > 169*9
# [1] 1521 obs. in temp
temp_original = temp                  # MR modified

temp = temp[which(temp$variable %in% c('Center', 'condition', 'Disease.Subtype')),]

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

plot_df[which(plot_df$value == 0),'value'] = 'inactive'
plot_df[which(plot_df$value == 1),'value'] = 'active'
plot_df[which(plot_df$value == .5),'value'] = 'active'

# Reorder factor so that the annotation is broader

plot_df_Model_Interactions_unique = unique(plot_df$Model_Interactions)            # MR inserted
plot_df_Model_Interactions_unique                                                 # MR inserted
# ...
# [186] !GILENYA=KPCZ       
#        Center              Disease.Subtype     condition            
#        Center2             Disease.Subtype2    condition2   

plot_df$Model_Interactions = factor(plot_df$Model_Interactions, levels = c(levels(plot_df$Model_Interactions)[1:186], 'Center', 'Center2', 'Disease.Subtype', 'Disease.Subtype2', 'condition', 'condition2'))

# reorder for legend

plot_df_value_unique = unique(plot_df$value)            # MR inserted
plot_df_value_unique                                                 # MR inserted
# [1] "inactive"  "active"    
#      "CH"        "IB"        "KI"        "UZ"      
#      "RR"        "healthy"   "CIS"      
# [10] "PP"        
#     "Treated"   "Healthy"   "Untreated"

# > levels(plot_df$value)
# NULL

plot_df$value = factor(plot_df$value, levels = c('Center', 'CH', 'KI', 'IB', 'UZ', '',
                                                 'Condition', 'Healthy', 'Treated', 'Untreated', ' ',
                                                 'Subtype', 'healthy', 'CIS', 'RR', 'PP', '  ',
                                                 'Reaction', 'active', 'inactive'))







# #########################################################################################################################
# Plot as Ggplot heatmap
# #########################################################################################################################

FigureWidth_map =5     # MR inserted
FigureHeight_map=8     # MR inserted

heatmap_median_models = ggplot(plot_df, aes(x=Model_Interactions, y=Patients)) + 
  facet_grid(~split, scales = "free_x", space='free', margins=FALSE) + 
  geom_tile(aes(fill=value)) + 
  xlab('Model Reactions (1-186)') + ylab('Patients (1-169)') +
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
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__A__heatmap_median_models.pdf",sep = "")), width = FigureWidth_map, height = FigureHeight_map)             # MR inserted




# ##########################################################################################################
# PART B OF FIGURE
# ##########################################################################################################

# ##########################################################################################################
# load annotation file and similarity BETWEEN matrix
# ####################################################################################################################

#load("../../files/similarity/similarityMatrix.RData")    # MR commented


files_Similarity_MR_folder = "../../files/Similarity_MR"                                                             # MR inserted
sub_dir__Results_Similarity_MR = file.path(files_Similarity_MR_folder,Results_storage_name)                          # MR inserted

load(paste(sub_dir__Results_Similarity_MR,"/similarityMatrix.RData",sep="")) # called similarityMatrix               # MR inserted

# wrong_similarityMatrix
load(paste(sub_dir__Results_Similarity_MR,"/wrong_similarityMatrix.RData",sep="")) # called wrong_similarityMatrix               # MR inserted


# ************************************ 
# include similarities WITHIN




# START MR modified

# load('../../files/similarity/similarityWITHINMeans.RData')
# for (i in 1:dim(similarityMatrix)[1]){
#    similarityMatrix[i,i]=means_within[i]
# }



if(include__similarityMatrixSelf__of_each_donor==TRUE){                                   # MR inserted
  
  #load('/Users/marti/Documents/R/combiMS/cluster/analysis/similarityWITHIN.RData')               # MR modified
  load(paste(sub_dir__Results_Similarity_MR,"mean_SimilarityBetweenDonorOwnModels.RData",sep="/")) # called mean_SimilarityBetweenDonorOwnModels               # MR modified
  
  #colnames(similarityMatrix)==names(similarityVector)
  for (i in 1:dim(similarityMatrix)[1]){
    #similarityMatrix[i,i]=similarityVector[i]
    similarityMatrix[i,i]=mean_SimilarityBetweenDonorOwnModels[i]             # MR modified
  }
}

# END MR modified





# ##########################################################################################################
# create dataframe of similarities to plot
# ##########################################################################################################
sim_df=melt(similarityMatrix)
sim_df$Group=rep(F,dim(similarityMatrix)[1])
sim_df$GroupSim=rep(F,dim(similarityMatrix)[1])
sim_df$Condition=rep(F,dim(similarityMatrix)[1])
sim_df$ConditionSim=rep(F,dim(similarityMatrix)[1])
sim_df$Self=rep(F,dim(similarityMatrix)[1])
sim_df$RR_Untreated=rep(F,dim(similarityMatrix)[1])

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

sub_same_group_similarity=subset(sim_df,sim_df$Group==sim_df$GroupSim)
sub_self_similarity=subset(sub_same_group_similarity,sub_same_group_similarity$Self=='Self')
sub_IFN=subset(sub_same_group_similarity,sub_same_group_similarity$RR_Untreated=='RR Untreated')

# delete Tysabri or Placebo, add All similarities, self similarities, and RR_untreated as reference
sub_same_group_similarity = sub_same_group_similarity[-which(sub_same_group_similarity$Group == 'Tysabri or Placebo'),]
sub_self_similarity$Group = 'Self'
sub_IFN$Group = 'RR Untreated'
sub_all = sim_df
sub_all$Group = 'All'

# Create additional variable for grouping of boxplots
sub_same_group_similarity$Grid = 'Treatment'
sub_same_group_similarity$Grid[which(sub_same_group_similarity$Group == 'Healthy')] = 'Status'
sub_same_group_similarity$Grid[which(sub_same_group_similarity$Group == 'Untreated')] = 'Status'
sub_self_similarity$Grid = 'All'
sub_IFN$Grid = 'Status'
sub_all$Grid = 'All'







# boxplot by Group
# In a notched box plot, the notches extend 1.58 * IQR / sqrt(n). This gives a roughly 95 interval for comparing medians
# So non-overlapping notches strongly suggest that medians are significantly different
# This is the case for RR untreated vis RR IFNb

plot_df_2 = rbind(sub_same_group_similarity, sub_self_similarity, sub_IFN, sub_all)


plot_df_2_Group__unique = unique(plot_df_2$Group)                 # MR modified
plot_df_2_Group__unique
#[1] "GA"           "Healthy"      "Untreated"    "IFNb"         "FTY"          "EGCG"         "NTZ"          "Self"         "RR Untreated" "All"         


## Reorder Groups
sorted_groups = c('All', 'Self', 'Healthy', 'EGCG', 'FTY', 'IFNb', 'GA', 'NTZ', 'RR Untreated', 'Untreated')
plot_df_2$Group = factor(plot_df_2$Group, sorted_groups)

FigureWidth_boxplot = 6   # MR inserted
FigureHeight_boxplot = 4  # MR inserted


if(include__similarityMatrixSelf__of_each_donor==TRUE){
  color_vec_group_fill = c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')
}else{
  color_vec_group_fill = c('#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')
}


boxplot_similarity = ggplot(plot_df_2, aes(x=Group, y=Similarity, fill=Group)) + 
  geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  scale_fill_manual(values=color_vec_group_fill) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             # MR inserted



# ###################################################################################################################################################
# Divert output into pdf
# ###################################################################################################################################################

#pdf("../../figures/figure_similarity.pdf", width=7, height=3.5, onefile = FALSE)            # MR modified
pdf(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__A__heatmap_median_models__B__boxplot_similarity.pdf",sep = "")), width=7, height=3.5, onefile = FALSE)

plot_grid(heatmap_median_models, boxplot_similarity, labels=c('A', 'B'), label_size = 14)

dev.off()














# START MR inserted


# ##########################################################################################################
# 
# Testing whether 
# 
# the whole symmetric similarityMatrix is plotted in the boxplots or 
# just the unique part of the similarityMatrix (diagonal entries and upper part (or lower) triangular part)
# 
# by considering the subgroup EGCG (sample size = 6 patients):
# 
# all combinations: 6*6 = 36
# unique combinations: (n(n+1))/2 = 21
# 
# ############################################################################################################



sub_plot_df_2_EGCG=subset(plot_df_2,plot_df_2$Group=='EGCG')


boxplot_similarity_EGCG = ggplot(sub_plot_df_2_EGCG, aes(x=Group, y=Similarity, fill=Group)) + 
  geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity_EGCG.pdf",sep = "")), width = FigureWidth_boxplot*0.2, height = FigureHeight_boxplot)             # MR inserted



boxplot_similarity_EGCG = ggplot(sub_plot_df_2_EGCG, aes(x=Group, y=Similarity, colour=SimWith)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_point(position = position_jitter(width = 0.5))+
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  #scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity_EGCG__geom_point_colour_SimWith.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             # MR inserted




boxplot_similarity_EGCG_2 = ggplot(sub_plot_df_2_EGCG, aes(x=Group, y=Similarity, colour=Patient)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_point(position = position_jitter(width = 0.5))+
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  #scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity_EGCG__geom_point_colour_Patient.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             # MR inserted


mat__donor_combinations_EGCG = sub_plot_df_2_EGCG[,c("SimWith","Patient")]


#EGCG_unique_indices=unique(t(apply(mat__donor_combinations_EGCG, 1, sort)))
EGCG_unique_indices=unique(t(apply(mat__donor_combinations_EGCG, 1, sort,decreasing=T)))
#EGCG_unique_indices = mat__donor_combinations_EGCG[!duplicated(t(apply(mat__donor_combinations_EGCG, 1, sort))),]
row.names(EGCG_unique_indices)    
sub_plot_df_2_EGCG_unique = sub_plot_df_2_EGCG[row.names(EGCG_unique_indices), ]




boxplot_similarity_EGCG_unique = ggplot(sub_plot_df_2_EGCG_unique, aes(x=Group, y=Similarity, fill=Group)) + 
  geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity_EGCG_unique.pdf",sep = "")), width = FigureWidth_boxplot*0.2, height = FigureHeight_boxplot)             # MR inserted



boxplot_similarity_EGCG_unique = ggplot(sub_plot_df_2_EGCG_unique, aes(x=Group, y=Similarity, colour=SimWith)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_point(position = position_jitter(width = 0.5))+
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  #scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity_EGCG_unique__geom_point_colour_SimWith.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             # MR inserted




boxplot_similarity_EGCG_unique_2 = ggplot(sub_plot_df_2_EGCG_unique, aes(x=Group, y=Similarity, colour=Patient)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_point(position = position_jitter(width = 0.5))+
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  #scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity_EGCG_unique__geom_point_colour_Patient.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             # MR inserted























# ###################################################################################################################################################
# 
# Reduce plot_df_2 to unique similarity combinations
# 
# ###################################################################################################################################################


# The following try between the lines TTTTTT does not work



# TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

# 
# mat__donor_combinations = plot_df_2[,c("SimWith","Patient")]
# 
# unique_indices = mat__donor_combinations[!duplicated(t(apply(mat__donor_combinations, 1, sort))),]
# #row.names(unique_indices)    
# plot_df_2_unique = plot_df_2[row.names(unique_indices), ]
# 
# 
# 
# 
# 
# 
# 
# boxplot_similarity = ggplot(plot_df_2_unique, aes(x=Group, y=Similarity, fill=Group)) + 
#    geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
#    theme_classic() + 
#    #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
#    scale_fill_manual(values=color_vec_group_fill) +                        # MR modified
#    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
#    theme(legend.position = 'none') + 
#    #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
#    geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
#    geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
#    ylim(-0.05,1.2)               # MR modified
# ggsave(file.path(sub_dir__Figures_folder_of_current_script,paste("figure_similarity_MR__B__boxplot_similarity_unique.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             # MR inserted
# 
# 
# 
# # ###################################################################################################################################################
# # Divert output into pdf
# # ###################################################################################################################################################
# 
# #pdf("../../figures/figure_similarity.pdf", width=7, height=3.5, onefile = FALSE)            # MR modified
# pdf(file.path(sub_dir__Figures_folder_of_current_script,paste("figure_similarity_MR__A__heatmap_median_models__B__boxplot_similarity_unique.pdf",sep = "")), width=7, height=3.5, onefile = FALSE)
# 
# plot_grid(heatmap_median_models, boxplot_similarity, labels=c('A', 'B'), label_size = 14)
# 
# dev.off()


# TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT



rm(list = ls())





Results_storage_name = "OptCombiMSclusterFork10_MR__using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__gaBinaryT1_postCellNOptRupdate"        # MR inserted

include__similarityMatrixSelf__of_each_donor = FALSE                                 # MR inserted



# ###############################################################################################################################################
# Based on script: analyzePatientSimilarity.R by Marti July 2015
#
# Minor modification for nicer plotting by Jakob in December 2016
# #################################################################################################################################################

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





figures_MR_folder = "../../figures_MR"   


ifelse(!dir.exists(file.path(figures_MR_folder,"figure_similarity_MR")), dir.create(file.path(figures_MR_folder,"figure_similarity_MR")), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script_class = file.path(figures_MR_folder,"figure_similarity_MR")                                                                              # MR inserted

ifelse(!dir.exists(file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)), dir.create(file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script = file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)                                                                              # MR inserted

ifelse(!dir.exists(file.path(sub_dir__Figures_folder_of_current_script,"Whole_similarityMatrix")), dir.create(file.path(sub_dir__Figures_folder_of_current_script,"Whole_similarityMatrix")), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix = file.path(sub_dir__Figures_folder_of_current_script,"Whole_similarityMatrix")     

ifelse(!dir.exists(file.path(sub_dir__Figures_folder_of_current_script,"Unique_similarityMatrix")), dir.create(file.path(sub_dir__Figures_folder_of_current_script,"Unique_similarityMatrix")), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix = file.path(sub_dir__Figures_folder_of_current_script,"Unique_similarityMatrix")     


# ###########################################################################################################
# Load and clean annotation data
# ###########################################################################################################
annot=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)

# ************************************ 
# add some fields to annotation for ggplot
annot2=annot
# correct that PPMS patients appear in group instead of by untreated and treatment





# Rename groups for consistency with other figures in the paper

annot2_Group_unique = unique(annot2$Group)                   # MR inserted
annot2_Group_unique
# [1] "Glatiramer"   "Healthy"      "Untreated"    "PPMS"         "Interferon B" "Fingolimod"  
#  "EGCG"         "Natalizumab" 

annot2_Disease.Subtype_unique = unique(annot2$Disease.Subtype)                   # MR inserted
annot2_Disease.Subtype_unique
# [1] "RR"  ""    "CIS" "PP"  "SP"  "PR" 






annot2[which(annot2$Group=='PPMS' & annot2$Treatment=='no'),'Group']='Untreated'                                                                # 1. PPMS correction
annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Group']=annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Treatment']   # 2. PPMS correction
annot2[which(annot2$Group=='Gilenya'),'Group']='FTY'                    # Gilenya  is the trade name of Fingolimod=FTY

annot2[which(annot2$Group=='Fingolimod'),'Group']='FTY'
annot2[which(annot2$Group=='Glatiramer'),'Group']='GA'
annot2[which(annot2$Group=='Interferon B'),'Group']='IFNb'
annot2[which(annot2$Group=='Natalizumab'),'Group']='NTZ'

annot2[which(annot2$Disease.Subtype == ''), 'Disease.Subtype'] = 'healthy'
annot2[which(annot2$Disease.Subtype == 'SP'), 'Disease.Subtype'] = substring(annot2[which(annot2$Disease.Subtype == 'SP'), 'Category'], 1, 2)
annot2[which(annot2$Disease.Subtype == 'PR'), 'Disease.Subtype'] = substring(annot2[which(annot2$Disease.Subtype == 'PR'), 'Category'], 1, 2)
annot=annot2
rm(annot2)





annot_Group_unique = unique(annot$Group)                   # MR inserted
annot_Group_unique
# [1] "GA"                 "Healthy"            "Untreated"          "IFNb"               "FTY"               
# [6] "EGCG"               "NTZ"                "Tysabri or Placebo"


annot_Disease.Subtype_unique = unique(annot$Disease.Subtype)                   # MR inserted
annot_Disease.Subtype_unique
# [1] "RR"      "healthy" "CIS"     "PP"     




# ##########################################################################################################
# PART A OF FIGURE
# ##########################################################################################################

#load('../../files/median_models/allMedianModels.RData')        # MR modified

files_median_models_MR_folder = "../../files/median_models_MR"                                       # MR inserted
sub_dir__Results_allMedianModels = file.path(files_median_models_MR_folder,Results_storage_name)     # MR inserted

load(file=paste0(sub_dir__Results_allMedianModels,"/allMedianModels.RData"))   # MR modified



# Merge annotation data (twice, so that it is broader) with model interactions
plot_df = melt(allMedianNetworks)

temp = melt(annot, id.vars = 'ID')
# > 169*9
# [1] 1521 obs. in temp
temp_original = temp                  # MR modified

temp = temp[which(temp$variable %in% c('Center', 'condition', 'Disease.Subtype')),]

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

plot_df[which(plot_df$value == 0),'value'] = 'inactive'
plot_df[which(plot_df$value == 1),'value'] = 'active'
plot_df[which(plot_df$value == .5),'value'] = 'active'

# Reorder factor so that the annotation is broader

plot_df_Model_Interactions_unique = unique(plot_df$Model_Interactions)            # MR inserted
plot_df_Model_Interactions_unique                                                 # MR inserted
# ...
# [186] !GILENYA=KPCZ       
#        Center              Disease.Subtype     condition            
#        Center2             Disease.Subtype2    condition2   

plot_df$Model_Interactions = factor(plot_df$Model_Interactions, levels = c(levels(plot_df$Model_Interactions)[1:186], 'Center', 'Center2', 'Disease.Subtype', 'Disease.Subtype2', 'condition', 'condition2'))

# reorder for legend

plot_df_value_unique = unique(plot_df$value)            # MR inserted
plot_df_value_unique                                                 # MR inserted
# [1] "inactive"  "active"    
#      "CH"        "IB"        "KI"        "UZ"      
#      "RR"        "healthy"   "CIS"      
# [10] "PP"        
#     "Treated"   "Healthy"   "Untreated"

# > levels(plot_df$value)
# NULL

plot_df$value = factor(plot_df$value, levels = c('Center', 'CH', 'KI', 'IB', 'UZ', '',
                                                 'Condition', 'Healthy', 'Treated', 'Untreated', ' ',
                                                 'Subtype', 'healthy', 'CIS', 'RR', 'PP', '  ',
                                                 'Reaction', 'active', 'inactive'))







# #########################################################################################################################
# Plot as Ggplot heatmap
# #########################################################################################################################

FigureWidth_map =5     # MR inserted
FigureHeight_map=8     # MR inserted

heatmap_median_models = ggplot(plot_df, aes(x=Model_Interactions, y=Patients)) + 
  facet_grid(~split, scales = "free_x", space='free', margins=FALSE) + 
  geom_tile(aes(fill=value)) + 
  xlab('Model Reactions (1-186)') + ylab('Patients (1-169)') +
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
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Whole_similarityMatrix,paste("figure_similarity_MR__A__heatmap_median_models.pdf",sep = "")), width = FigureWidth_map, height = FigureHeight_map)             # MR inserted




# ##########################################################################################################
# PART B OF FIGURE
# ##########################################################################################################

# ##########################################################################################################
# load annotation file and similarity BETWEEN matrix
# ####################################################################################################################

#load("../../files/similarity/similarityMatrix.RData")    # MR commented


files_Similarity_MR_folder = "../../files/Similarity_MR"                                                             # MR inserted
sub_dir__Results_Similarity_MR = file.path(files_Similarity_MR_folder,Results_storage_name)                          # MR inserted

load(paste(sub_dir__Results_Similarity_MR,"/similarityMatrix.RData",sep="")) # called similarityMatrix               # MR inserted

# wrong_similarityMatrix
load(paste(sub_dir__Results_Similarity_MR,"/wrong_similarityMatrix.RData",sep="")) # called wrong_similarityMatrix               # MR inserted


# ************************************ 
# include similarities WITHIN




# START MR modified

# load('../../files/similarity/similarityWITHINMeans.RData')
# for (i in 1:dim(similarityMatrix)[1]){
#    similarityMatrix[i,i]=means_within[i]
# }



if(include__similarityMatrixSelf__of_each_donor==TRUE){                                   # MR inserted
  
  #load('/Users/marti/Documents/R/combiMS/cluster/analysis/similarityWITHIN.RData')               # MR modified
  load(paste(sub_dir__Results_Similarity_MR,"mean_SimilarityBetweenDonorOwnModels.RData",sep="/")) # called mean_SimilarityBetweenDonorOwnModels               # MR modified
  
  #colnames(similarityMatrix)==names(similarityVector)
  for (i in 1:dim(similarityMatrix)[1]){
    #similarityMatrix[i,i]=similarityVector[i]
    similarityMatrix[i,i]=mean_SimilarityBetweenDonorOwnModels[i]             # MR modified
  }
}

# END MR modified











similarityMatrix_unique = similarityMatrix

similarityMatrix_unique[lower.tri(similarityMatrix_unique,diag=F)]=NA                      # MR inserted






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

sub_same_group_similarity=subset(sim_df,sim_df$Group==sim_df$GroupSim)
sub_self_similarity=subset(sub_same_group_similarity,sub_same_group_similarity$Self=='Self')
sub_IFN=subset(sub_same_group_similarity,sub_same_group_similarity$RR_Untreated=='RR Untreated')

# delete Tysabri or Placebo, add All similarities, self similarities, and RR_untreated as reference
sub_same_group_similarity = sub_same_group_similarity[-which(sub_same_group_similarity$Group == 'Tysabri or Placebo'),]
sub_self_similarity$Group = 'Self'
sub_IFN$Group = 'RR Untreated'
sub_all = sim_df
sub_all$Group = 'All'

# Create additional variable for grouping of boxplots
sub_same_group_similarity$Grid = 'Treatment'
sub_same_group_similarity$Grid[which(sub_same_group_similarity$Group == 'Healthy')] = 'Status'
sub_same_group_similarity$Grid[which(sub_same_group_similarity$Group == 'Untreated')] = 'Status'
sub_self_similarity$Grid = 'All'
sub_IFN$Grid = 'Status'
sub_all$Grid = 'All'







# boxplot by Group
# In a notched box plot, the notches extend 1.58 * IQR / sqrt(n). This gives a roughly 95 interval for comparing medians
# So non-overlapping notches strongly suggest that medians are significantly different
# This is the case for RR untreated vis RR IFNb

plot_df_2 = rbind(sub_same_group_similarity, sub_self_similarity, sub_IFN, sub_all)


plot_df_2_Group__unique = unique(plot_df_2$Group)                 # MR modified
plot_df_2_Group__unique
#[1] "GA"           "Healthy"      "Untreated"    "IFNb"         "FTY"          "EGCG"         "NTZ"          "Self"         "RR Untreated" "All"         


## Reorder Groups
sorted_groups = c('All', 'Self', 'Healthy', 'EGCG', 'FTY', 'IFNb', 'GA', 'NTZ', 'RR Untreated', 'Untreated')
plot_df_2$Group = factor(plot_df_2$Group, sorted_groups)

FigureWidth_boxplot = 6   # MR inserted
FigureHeight_boxplot = 4  # MR inserted


if(include__similarityMatrixSelf__of_each_donor==TRUE){
  color_vec_group_fill = c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')
}else{
  color_vec_group_fill = c('#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')
}


boxplot_similarity = ggplot(plot_df_2, aes(x=Group, y=Similarity, fill=Group)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         # MR modified
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  scale_fill_manual(values=color_vec_group_fill) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                                                 # MR modified
  #geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05,na.rm = TRUE) +                                                 # MR modified
  geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15,na.rm = TRUE) +                                                     # MR modified
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity__unique.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             # MR inserted








boxplot_similarity_test_all_drugs__versus__untreated = ggplot(plot_df_2, aes(x=Group, y=Similarity, fill=Group)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         # MR modified
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  scale_fill_manual(values=color_vec_group_fill) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                                                 # MR modified
  #geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  geom_signif(comparisons=list(c("EGCG", "Untreated")), map_signif_level = TRUE, y_position = 0.7,na.rm = TRUE) +                                                 # MR modified
  geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 0.8,na.rm = TRUE) +                                                     # MR modified
  geom_signif(comparisons=list(c("IFNb", "Untreated")), map_signif_level = TRUE, y_position = 0.9,na.rm = TRUE) +                                                 # MR modified
  geom_signif(comparisons=list(c("GA", "Untreated")), map_signif_level = TRUE, y_position = 1.0,na.rm = TRUE) +  
  geom_signif(comparisons=list(c("NTZ", "Untreated")), map_signif_level = TRUE, y_position = 1.1,na.rm = TRUE) +                                                     # MR modified
  # MR modified
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity_test_all_drugs__versus__untreated__unique.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             # MR inserted

boxplot_similarity_test_all_drugs__versus__RRuntreated = ggplot(plot_df_2, aes(x=Group, y=Similarity, fill=Group)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         # MR modified
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  scale_fill_manual(values=color_vec_group_fill) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                                                 # MR modified
  #geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  geom_signif(comparisons=list(c("EGCG", "RR Untreated")), map_signif_level = TRUE, y_position = 0.7,na.rm = TRUE) +                                                 # MR modified
  geom_signif(comparisons=list(c("FTY", "RR Untreated")), map_signif_level = TRUE, y_position = 0.8,na.rm = TRUE) +                                                     # MR modified
  geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 0.9,na.rm = TRUE) +                                                 # MR modified
  geom_signif(comparisons=list(c("GA", "RR Untreated")), map_signif_level = TRUE, y_position = 1.0,na.rm = TRUE) +  
  geom_signif(comparisons=list(c("NTZ", "RR Untreated")), map_signif_level = TRUE, y_position = 1.1,na.rm = TRUE) +                                                     # MR modified
  # MR modified
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity_test_all_drugs__versus__RRuntreated__unique.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             # MR inserted














plot_df_2__subset_without_NA = subset(plot_df_2, !is.na(plot_df_2$Similarity))
plot_df_2__subset_without_NA_Group_unique = unique(plot_df_2__subset_without_NA$Group)
plot_df_2__subset_without_NA_Group_unique
# [1] Untreated    IFNb         EGCG         GA           FTY          Healthy      NTZ          RR Untreated All         
# Levels: All Self Healthy EGCG FTY IFNb GA NTZ RR Untreated Untreated

boxplot_similarity_subset_without_NA = ggplot(plot_df_2__subset_without_NA, aes(x=Group, y=Similarity, fill=Group)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         # MR modified
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  scale_fill_manual(values=color_vec_group_fill) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                                                 # MR modified
  #geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05,na.rm = TRUE) +                                                 # MR modified
  geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15,na.rm = TRUE) +                                                     # MR modified
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity__unique__using_plot_df_2__subset_without_NA.pdf",sep = "")), width = FigureWidth_boxplot, height = FigureHeight_boxplot)             # MR inserted



# ###################################################################################################################################################
# Divert output into pdf
# ###################################################################################################################################################

#pdf("../../figures/figure_similarity.pdf", width=7, height=3.5, onefile = FALSE)            # MR modified
pdf(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__A__heatmap_median_models__B__boxplot_similarity.pdf",sep = "")), width=7, height=3.5, onefile = FALSE)

plot_grid(heatmap_median_models, boxplot_similarity, labels=c('A', 'B'), label_size = 14)

dev.off()














# START MR inserted


# ##########################################################################################################
# 
# Testing whether 
# 
# the whole symmetric similarityMatrix_unique is plotted in the boxplots or 
# only the unique part of the similarityMatrix_unique (diagonal entries and upper part (or lower) triangular part)
# 
# by considering the subgroup EGCG (sample size = 6 patients):
# 
# all combinations: 6*6 = 36
# unique combinations: (n(n+1))/2 = 21
# 
# ############################################################################################################



sub_plot_df_2_EGCG=subset(plot_df_2,plot_df_2$Group=='EGCG')


boxplot_similarity_EGCG = ggplot(sub_plot_df_2_EGCG, aes(x=Group, y=Similarity, fill=Group)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         # MR modified
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity__unique_EGCG.pdf",sep = "")), width = FigureWidth_boxplot*0.2, height = FigureHeight_boxplot)             # MR inserted



boxplot_similarity_EGCG = ggplot(sub_plot_df_2_EGCG, aes(x=Group, y=Similarity, colour=SimWith)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_point(position = position_jitter(width = 0.5))+
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  #scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity__unique_EGCG__geom_point_colour_SimWith.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             # MR inserted




boxplot_similarity_EGCG_2 = ggplot(sub_plot_df_2_EGCG, aes(x=Group, y=Similarity, colour=Patient)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_point(position = position_jitter(width = 0.5))+
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  #scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity__unique_EGCG__geom_point_colour_Patient.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             # MR inserted


mat__donor_combinations_EGCG = sub_plot_df_2_EGCG[,c("SimWith","Patient")]

#EGCG_unique_indices=unique(t(apply(mat__donor_combinations_EGCG, 1, sort)))
EGCG_unique_indices=unique(t(apply(mat__donor_combinations_EGCG, 1, sort,decreasing=T)))

#EGCG_unique_indices = mat__donor_combinations_EGCG[!duplicated(t(apply(mat__donor_combinations_EGCG, 1, sort))),]
row.names(EGCG_unique_indices)    
sub_plot_df_2_EGCG_unique = sub_plot_df_2_EGCG[row.names(EGCG_unique_indices), ]




boxplot_similarity_EGCG_unique = ggplot(sub_plot_df_2_EGCG_unique, aes(x=Group, y=Similarity, fill=Group)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0,na.rm = TRUE) +                                                                         # MR modified
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity__unique_EGCG_unique.pdf",sep = "")), width = FigureWidth_boxplot*0.2, height = FigureHeight_boxplot)             # MR inserted



boxplot_similarity_EGCG_unique = ggplot(sub_plot_df_2_EGCG_unique, aes(x=Group, y=Similarity, colour=SimWith)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_point(position = position_jitter(width = 0.5))+
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  #scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity__unique_EGCG_unique__geom_point_colour_SimWith.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             # MR inserted




boxplot_similarity_EGCG_unique_2 = ggplot(sub_plot_df_2_EGCG_unique, aes(x=Group, y=Similarity, colour=Patient)) + 
  #geom_boxplot(notch=T, outlier.alpha = 0.2, outlier.stroke = 0) +
  geom_point(position = position_jitter(width = 0.5))+
  theme_classic() + 
  #scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) +         
  #scale_fill_manual(values=c('#8DC060')) +                        # MR modified
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #theme(legend.position = 'none') + 
  #geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE) +
  # geom_signif(comparisons=list(c("IFNb", "RR Untreated")), map_signif_level = TRUE, y_position = 1.05) +                         # MR modified
  # geom_signif(comparisons=list(c("FTY", "Untreated")), map_signif_level = TRUE, y_position = 1.15) +
  ylim(-0.05,1.2)               # MR modified
ggsave(file.path(sub_dir__Figures_folder_of_current_script__Unique_similarityMatrix,paste("figure_similarity_MR__B__boxplot_similarity__unique_EGCG_unique__geom_point_colour_Patient.pdf",sep = "")), width = FigureWidth_boxplot*1.5, height = FigureHeight_boxplot)             # MR inserted

















print("Script finished!")


# END MR inserted



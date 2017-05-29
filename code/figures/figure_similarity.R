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
  
# Use relative paths instead of absolut paths for the files
# in Rstudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ###########################################################################################################
# Load and clean annotation data
# ###########################################################################################################
annot=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)

# ************************************ 
# add some fields to annotation for ggplot
annot2=annot
# correct that PPMS patients appear in group instead of by untreated and treatment
# Rename groups for consistency with other figures in the paper
annot2[which(annot2$Group=='PPMS' & annot2$Treatment=='no'),'Group']='Untreated'
annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Group']=annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Treatment']
annot2[which(annot2$Group=='Gilenya'),'Group']='FTY'
annot2[which(annot2$Group=='Fingolimod'),'Group']='FTY'
annot2[which(annot2$Group=='Glatiramer'),'Group']='GA'
annot2[which(annot2$Group=='Interferon B'),'Group']='IFNb'
annot2[which(annot2$Group=='Natalizumab'),'Group']='NTZ'
annot2[which(annot2$Disease.Subtype == ''), 'Disease.Subtype'] = 'healthy'
annot2[which(annot2$Disease.Subtype == 'SP'), 'Disease.Subtype'] = substring(annot2[which(annot2$Disease.Subtype == 'SP'), 'Category'], 1, 2)
annot2[which(annot2$Disease.Subtype == 'PR'), 'Disease.Subtype'] = substring(annot2[which(annot2$Disease.Subtype == 'PR'), 'Category'], 1, 2)
annot=annot2
rm(annot2)

# ##########################################################################################################
# PART A OF FIGURE
# ##########################################################################################################
load('../../files/median_models/allMedianModels.RData')

# Merge annotation data (twice, so that it is broader) with model interactions
plot_df = melt(allMedianNetworks)
temp = melt(annot, id.vars = 'ID')
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
plot_df$Model_Interactions = factor(plot_df$Model_Interactions, levels = c(levels(plot_df$Model_Interactions)[1:186], 'Center', 'Center2', 'Disease.Subtype', 'Disease.Subtype2', 'condition', 'condition2'))

# reorder for legend
plot_df$value = factor(plot_df$value, levels = c('Center', 'CH', 'KI', 'IB', 'UZ', '',
                                                 'Condition', 'Healthy', 'Treated', 'Untreated', ' ',
                                                 'Subtype', 'healthy', 'CIS', 'RR', 'PP', '  ',
                                                 'Reaction', 'active', 'inactive'))

# #########################################################################################################################
# Plot as Ggplot heatmap
# #########################################################################################################################

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

# ##########################################################################################################
# PART B OF FIGURE
# ##########################################################################################################

# ##########################################################################################################
# load annotation file and similarity BETWEEN matrix
# ####################################################################################################################
load("../../files/similarity/similarityMatrix.RData")

# ************************************ 
# include similarities WITHIN
load('../../files/similarity/similarityWITHINMeans.RData')
for (i in 1:dim(similarityMatrix)[1]){
  similarityMatrix[i,i]=means_within[i]
}
 

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

## Reorder Groups
sorted_groups = c('All', 'Self', 'Healthy', 'EGCG', 'FTY', 'IFNb', 'GA', 'NTZ', 'RR Untreated', 'Untreated')
plot_df_2$Group = factor(plot_df_2$Group, sorted_groups)

boxplot_similarity = ggplot(plot_df_2) + 
  geom_boxplot(notch=T,aes(x=Group, y=Similarity, fill=Group), outlier.alpha = 0.2, outlier.stroke = 0) +
  theme_classic() + 
  scale_fill_manual(values=c('#ECEDED', '#ECEDED', '#407FB7', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#8DC060', '#FDD48F', '#FABE50')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  ylim(c(0.5, 1))

# ###################################################################################################################################################
# Divert output into pdf
# ###################################################################################################################################################
pdf("../../figures/figure_similarity.pdf", width=7, height=3.5, onefile = FALSE)

plot_grid(heatmap_median_models, boxplot_similarity, labels=c('A', 'B'), label_size = 14)

dev.off()
